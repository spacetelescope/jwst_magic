#!/usr/bin/env python

"""
segment_guiding.py

Segment Guiding Tool (SGT), optionally using tkinter GUI

Created by Colin Cox on 2017-05-01.
Modified by Lauren Chambers January 2018
"""

# Standard Library Imports
import os
import logging

# Third Party Imports
# Work out backend business
import matplotlib
backend = matplotlib.get_backend()
matplotlib.use("Qt5Agg")
import matplotlib.pyplot as plt
import pysiaf
from pysiaf.utils import rotations
import numpy as np
from astropy.io import ascii as asc
from astropy.table import Table

# Local Imports
from .. import coordinate_transforms, utils
from ..segment_guiding import SegmentGuidingGUI

# Establish segment guiding files directory
LOCAL_PATH = os.path.dirname(os.path.realpath(__file__))
PACKAGE_PATH = os.path.split(LOCAL_PATH)[0]
OUT_PATH = os.path.split(PACKAGE_PATH)[0]  # Location of out/ and logs/ directory

# Open the SIAF with pysiaf
FGS_SIAF = pysiaf.Siaf('FGS')

# Start logger
LOGGER = logging.getLogger(__name__)


class SegmentGuidingCalculator:
    def __init__(self, segment_infile, program_id, observation_num, visit_num,
                 root=None, GUI=False, GS_params_dict=None, refonly=False,
                 selected_segs=None, vss_infile=None, out_dir=None,
                 ct_uncert_fctr=0.9, countrate_factor=None):

        self.root = root
        self.out_dir = utils.make_out_dir(out_dir, OUT_PATH, self.root)

        utils.ensure_dir_exists(self.out_dir)

        self.program_id = program_id
        self.observation_num = observation_num
        self.visit_num = visit_num

        # Will the tool be run through the GUI?
        self.GUI = GUI

        # Implement "refonly" label for reference stars?
        self.refonly = refonly

        # Set countrate uncertainty factor and multiplicative factor
        self.ct_uncert_fctr = ct_uncert_fctr
        self.countrate_factor = countrate_factor

        # Parse the input file type (ALLpsfs.txt, regfile.txt, and VSS infile)
        self.get_gs_params(vss_infile, GS_params_dict)
        self.parse_infile(segment_infile)
        if selected_segs is not None:
            self.get_selected_segs(selected_segs)

        # Get aperture parameters from FGS SIAF
        self.FGSsetup()

    def ChosenSeg(self, *args):
        '''
        1) Check that the user-provided segment ID number is valid (0 to 18)
        2) Calculate the central V2/V3 point of either the provided segment or
        the center of the segment array
        3) Update GUI with V2/V3 coordinates and Ideal angle coordinatess
        '''
        # Try to convert provided segment ID number to an integer
        try:
            segN = int(self.segNum)
        except ValueError:
            raise ValueError('Unrecognized segment number: {}'.format(segN))
            return

        # If no segment ID number is provided, do nothing
        if segN == '':
            return

        # Refresh FGS data
        self.FGSsetup()

        # Ensure the provided segment ID is valid
        segMax = len(self.V2SegArray)
        if (segN < 0) or (segN > segMax):
            msg = 'Segment number {} out of range (0, {})'.format(segN, segMax)
            raise ValueError(msg)
            return

        # Determine the central V2/V3 point from the given segment ID

        # If a specific segment was provided, set the V2/V3 ref point to be
        # that segment's location
        if segN > 0:
            self.V2SegN = self.V2SegArray[segN - 1]
            self.V3SegN = self.V3SegArray[segN - 1]
        # Otherwise, if the input segment ID was 0, set the V2/V3 ref point to
        # be the mean of all segments' locations
        else:
            self.V2SegN = self.V2SegArray.mean()
            self.V3SegN = self.V3SegArray.mean()

        # Calculate aim
        self.V2Ref = float(self.fgsV2)   # obtained from FGS setup
        self.V3Ref = float(self.fgsV3)
        dV2Aim = self.V2SegN
        dV3Aim = self.V3SegN
        self.V2Aim = self.V2Ref + dV2Aim
        self.V3Aim = self.V3Ref + dV3Aim

        # Convert to Ideal coordinates
        self.xIdl, self.yIdl = self.fgs_siaf_aperture.tel_to_idl(self.V2Aim, self.V3Aim)

    def FGSsetup(self, *args):
        '''Taking the current guider number (per the radio buttons on the GUI),
        extracts V2Ref, V3Ref, V3IdlYAngle, and VIdlParity from the respective
        aperture in the FGS SIAF.
        '''

        # Read chosen guider number
        # Ensure the guider number is valid
        if str(self.fgsNum) not in ['1', '2']:
            raise ValueError('Invalid guider number: "{}"'.format(self.fgsNum))

        det = 'FGS' + str(self.fgsNum) + '_FULL_OSS'

        # Open SIAF aperture for appropriate guider with pysiaf
        self.fgs_siaf_aperture = FGS_SIAF[det]
        V2Ref = self.fgs_siaf_aperture.V2Ref
        V3Ref = self.fgs_siaf_aperture.V3Ref
        V3IdlYAngle = self.fgs_siaf_aperture.V3IdlYAngle
        VIdlParity = self.fgs_siaf_aperture.VIdlParity

        self.fgsV2 = '%10.4f' % V2Ref
        self.fgsV3 = '%10.4f' % V3Ref
        self.fgsAngle = '%10.4f' % V3IdlYAngle
        self.fgsParity = '%3d' % VIdlParity


    def Calculate(self):
        nseg = len(self.SegIDArray)

        # Convert V2/V3 coordinates to ideal coordinates
        self.xIdlSegs, self.yIdlSegs = self.fgs_siaf_aperture.tel_to_idl(self.V2SegArray + self.V2Ref,
                                                                         self.V3SegArray + self.V3Ref)

        # Verify all guidestar parameters are valid
        self.check_guidestar_params()

        # Get the attitude matrix
        A = rotations.attitude(self.V2Aim + self.V2Boff, self.V3Aim + self.V3Boff,
                               self.RA, self.Dec, float(self.PA))

        # Get RA and Dec for each segment.
        self.SegRA = np.zeros(nseg)
        self.SegDec = np.zeros(nseg)
        for i in range(nseg):
            V2 = self.V2Ref + self.V2SegArray[i]
            V3 = self.V3Ref + self.V3SegArray[i]
            (self.SegRA[i], self.SegDec[i]) = rotations.pointing(A, V2, V3, positive_ra=True)

        # Convert segment coordinates to detector frame
        self.xDet, self.yDet = self.fgs_siaf_aperture.idl_to_det(self.xIdlSegs, self.yIdlSegs)

        # Check to make sure no segments are off the detector
        for x, y, i_seg in zip(self.xDet, self.yDet, self.SegIDArray):
            if x < 0.5 or x > 2048.5:
                LOGGER.warn('Segment Guiding: ' + '%8s off detector in X direction' % i_seg)
            if y < 0.5 or y > 2048.5:
                LOGGER.warn('Segment Guiding: ' + '%8s off detector in Y direction' % i_seg)

        # Check to make sure that RA is between 0 and 360 and Dec is between -90 and 90
        self.check_coords()

        # Print and save final output
        self.write_override_file(nseg)

        # Save .pngs of plots
        self.plot_segments()

    def write_override_file(self, nseg, verbose=True):
        # Define path and name of output override file
        out_file = 'gs-override-{}_{}_{}.txt'.format(self.program_id, self.observation_num,
                                                     self.visit_num)
        out_file = os.path.join(self.out_dir, out_file)

        # Print summary of input data (guide star RA, Dec, and PA, etc...)
        if verbose:
            # Get the guide star and boresight parameters
            V2B, V3B, gsRA, gsDec, gsPA = self.check_guidestar_params()

            summary_output = """Guide Star Parameters
                Aperture FGS: {0}
                V2/V3 Refs: ({1}, {2}) arc-sec
                Guiding segment number: {3}
                V2/V3 Boresight offset: ({4}, {5}) arc-sec
                Guide star RA & Dec: ({6}, {7}) degrees
                Position angle: {8} degrees""".format(self.fgsNum, self.fgsV2, self.fgsV3, self.segNum, V2B, V3B, gsRA, gsDec, gsPA)
            LOGGER.info('Segment Guiding: ' + summary_output)

            all_segments = 'All Segment Locations'
            all_segments += '\n                Segment     dV2    dV3    xIdl   yIdl     RA         Dec         xDet     yDet'
            for p in range(nseg):
                all_segments += ('\n                %5s    %6.2f %6.2f  %6.2f %6.2f  %10.6f %10.6f  %8.2f %8.2f' \
                                  % (self.SegIDArray[p], self.V2SegArray[p], self.V3SegArray[p],
                                     self.xIdlSegs[p], self.yIdlSegs[p], self.SegRA[p],
                                     self.SegDec[p], self.xDet[p], self.yDet[p]))
            LOGGER.info('Segment Guiding: ' + all_segments)

        # Write out override file with RA/Decs of selected segments
        with open(out_file, 'w') as f:
            # Determine whether to include a multiplicative countrate factor
            countrate_qualifier = \
                ' -count_rate_factor=({:.3f})'.format(self.countrate_factor) if self.countrate_factor else ''
            out_string = 'sts -gs_select{} {:4d}:{}:{}'.format(countrate_qualifier, int(self.program_id),
                                                             self.observation_num,
                                                             self.visit_num)

            # If segments have been selected, only use those
            try:
                orientations = list(self.selected_segment_ids)
                guide_segments = [s[0] for s in orientations]
                all_selected_segs = list(set(np.concatenate(orientations)))

                if self.refonly:
                    n_guide_segments = len(guide_segments)

                for i in range(len(all_selected_segs)):
                    # Reorder all possible segments
                    all_selected_segs.append(all_selected_segs[0])
                    all_selected_segs.remove(all_selected_segs[0])
                    new_seg = list(np.copy(all_selected_segs))

                    # Add to orientation list if not already provided as guide star
                    if new_seg[0] not in guide_segments:
                        orientations.append(new_seg)
                        guide_segments.append(new_seg[0])

            # If not, use all 18 segments as default
            except AttributeError:
                orientations = np.array([np.linspace(1, nseg, nseg).astype(int)])

            # If countrates were included in the input file, use them!
            try:
                rate = self.counts_array
                uncertainty = self.counts_array * self.ct_uncert_fctr
            except AttributeError:
                rate = [0.0] * nseg
                uncertainty = [0.0] * nseg

            # Write the commands for each orientation
            for i_o, orientation in enumerate(orientations):
                guide_seg_id = orientation[0]

                label = 'star'
                seg = i_o + 1
                # If implementing "ref_only" labels, determine if star is guide
                # or reference star, and alter label and ID accordingly
                if self.refonly:
                    if i_o >= n_guide_segments:
                        label = 'ref_only'
                        seg = i_o + 1 - n_guide_segments

                # Format segment properties (ID, RA, Dec, countrate, uncertainty)
                star_string = ' -%s%d = %d, %.6f, %.6f, %.1f, %.1f' % (
                    label, seg, guide_seg_id + 1, self.SegRA[guide_seg_id],
                    self.SegDec[guide_seg_id], rate[guide_seg_id],
                    uncertainty[guide_seg_id])

                if not self.refonly or (self.refonly and label == 'star'):
                    # Add list of segment IDs for all reference stars
                    for ref_seg_id in orientation:
                        if ref_seg_id != guide_seg_id:
                            star_string += ', %d' % (ref_seg_id + 1)

                out_string += star_string

            f.write(out_string)

            if verbose:
                LOGGER.info('Segment Guiding: ' + 'Guide Star Override: ' + \
                            out_string.replace('-star', '\n                -star').\
                            replace('-ref_only', '\n                -ref_only'))
            LOGGER.info('Segment Guiding: ' + 'Saved {} segment commands to {}'.format(len(orientations), out_file))

    def check_coords(self):
        ''' Check to make sure that RA is between 0 and 360 and Dec between -90 and 90
        '''
        for i, ra in enumerate(self.SegRA):
            if ra > 360.0:
                LOGGER.warn('Segment Guiding: ' + 'RA = {}'.format(ra))
                self.SegRA -= self.SegRA
            elif ra < 0.0:
                LOGGER.warn('Segment Guiding: ' + 'RA = {}'.format(ra))
                self.SegRA += 360.0
            else:
                continue

        for i, dec in enumerate(self.SegDec):
            if dec > 90.0:
                LOGGER.warn('Segment Guiding: ' + 'Dec = {}'.format(dec))
                self.SegDec -= 180.0
            elif dec < -90.0:
                LOGGER.warn('Segment Guiding: ' + 'Dec = {}'.format(dec))
                self.SegDec += 180.0
            else:
                continue

    def get_gs_params(self, vss_infile, GS_params_dict):
        '''Get GS parameters from dictionary or VSS file
        '''
        if vss_infile and GS_params_dict:
            LOGGER.info('Segment Guiding: ' + 'Reading RA, Dec, and PA from VSS file {}'.format(vss_infile))
            LOGGER.info('Segment Guiding: ' + 'Reading boresight offset and segment number from user-provided dictionary.')
            self.get_guidestar_params_from_visit_file(vss_infile)
            self.V2Boff = GS_params_dict['V2Boff']
            self.V3Boff = GS_params_dict['V3Boff']
            self.segNum = GS_params_dict['segNum']

        elif vss_infile:
            LOGGER.info('Segment Guiding: ' + 'Reading RA, Dec, and PA from VSS file {}'.format(vss_infile))
            LOGGER.info('Segment Guiding: ' + 'Setting boresight offset = 0 and segment number = 0.')
            self.get_guidestar_params_from_visit_file(vss_infile)
            self.V2Boff = 0
            self.V3Boff = 0
            self.segNum = 0

        elif GS_params_dict:
            LOGGER.info('Segment Guiding: ' + 'Reading all GS parameters from user-provided dictionary.')
            # Map GS_params_dict keys to attributes
            for attr_name in GS_params_dict.keys():
                setattr(self, attr_name, float(GS_params_dict[attr_name]))
            self.segNum = int(self.segNum)
            self.fgsNum = int(self.fgsNum)

        else:
            raise ValueError('If running the tool outside of the GUI, must '
                             'supply a dictionary of the required parameters '
                             'to the GS_params_dict argument and/or supply a '
                             'VSS file to the vss_infile argument.')

    def parse_infile(self, segment_infile):
        # If the input file is a .txt file, parse the file
        if segment_infile[-4:] == '.txt':
            read_table = asc.read(segment_infile)
            column_names = read_table.colnames
            n_segs = len(read_table)

            # Are the segment positions already in V2/V3?
            if (any(['V2Seg' == c for c in column_names])) and \
               (any(['V3Seg' == c for c in column_names])):
                segment_coords = read_table

            # If not, convert pixel x/y to V2/V3
            elif (any(['x' == c for c in column_names])) and \
                 (any(['y' == c for c in column_names])):
                segment_coords = Table()
                segment_coords['SegID'] = np.linspace(1, n_segs, n_segs).astype(int)
                v2, v3 = coordinate_transforms.Raw2Tel(read_table['x'], read_table['y'], self.fgsNum)
                segment_coords['V2Seg'], segment_coords['V3Seg'] = v2, v3

            else:
                raise TypeError('Incompatible file type: ', segment_infile)

            # If the countrates are included in the input file, read them!
            if (any(['countrate' == c for c in column_names])):
                self.counts_array = read_table['countrate']

        else:
            raise TypeError('Incompatible file type: ', segment_infile)

        # Define the IDs and coordinates of all segments
        self.SegIDArray = segment_coords['SegID']
        self.V2SegArray = segment_coords['V2Seg']
        self.V3SegArray = segment_coords['V3Seg']
        LOGGER.info('Segment Guiding: ' + '{} segment coordinates read from {}'.format(len(self.SegIDArray), segment_infile))

    def get_selected_segs(self, selected_segs):
        '''If a list of selected segments has been provided, figure that out
        '''

        # If the selected segments are an array of lists (passed from GUI)
        if isinstance(selected_segs, np.ndarray):
            LOGGER.info('Segment Guiding: Guiding on segments selected from GUI')
            # If there is more than one orientation provided
            if isinstance(selected_segs[0], list):
                self.selected_segment_ids = []
                for i, orientation in enumerate(selected_segs):
                    self.selected_segment_ids.append([s - 1 for s in orientation])

            # If there is only one
            else:
                self.selected_segment_ids = selected_segs - 1

        # Else if they are a regfile, parse it.
        elif os.path.exists(selected_segs):
            LOGGER.info('Segment Guiding: Guiding on segments selected in {}'.format(selected_segs))
            self.parse_regfile(selected_segs)

        # Otherwise, we don't know what it is...
        else:
            raise TypeError('Unrecognized data type passed to selected_segs ({}); must be regfile.txt path or array of indices.'.format(selected_segs))

    def parse_regfile(self, selected_segs):
        # Make sure the file is ascii-readable
        n_segs = len(self.SegIDArray)
        try:
            read_selected_segs = asc.read(selected_segs)
            column_names = read_selected_segs.colnames
            LOGGER.info('Segment Guiding: ' + 'Selected segment coordinates read from {}'.format(selected_segs))
        except:
            raise TypeError('Incompatible regfile type: ', selected_segs)

        # Are the segment positions already in V2/V3?
        if (any(['V2Seg' == c for c in column_names])) and \
           (any(['V3Seg' == c for c in column_names])):
            selected_segment_coords = read_selected_segs

        # If not, convert pixel x/y to V2/V3
        elif (any(['x' == c for c in column_names])) and \
             (any(['y' == c for c in column_names])):
            selected_segment_coords = Table()
            v2, v3 = coordinate_transforms.Raw2Tel(read_selected_segs['x'],
                                                   read_selected_segs['y'],
                                                   self.fgsNum)
            selected_segment_coords['V2Seg'], selected_segment_coords['V3Seg'] = v2, v3

        # If the positions aren't V2/V3 or x/y, what are they??
        else:
            raise TypeError('Incompatible regfile type: ', selected_segs)

        # Match locations of selected segments to IDs of known segments
        selected_segs_ids = []
        for coords in zip(selected_segment_coords['V2Seg', 'V3Seg']):
            for i_seg in range(n_segs):
                if coords[0]['V2Seg'] == self.V2SegArray[i_seg] and \
                   coords[0]['V3Seg'] == self.V3SegArray[i_seg]:
                    selected_segs_ids.append(i_seg)
        if len(selected_segs_ids) == 0:
            raise TypeError('Coordinates of selected segments file do not match those of the provided input file.')

        self.selected_segment_ids = [selected_segs_ids]

    def get_guidestar_params_from_visit_file(self, visit_file):
        # Open provided visit file, verify formatting, and get index of data start
        with open(visit_file) as vf:
            lines = vf.readlines()

            i_start = 0
            for i_line, line in enumerate(lines):
                if 'Corrected RA' in line:
                    i_start = i_line
                    break

            if ('Short Term Schedule Display' not in lines[0]) or \
               (i_start == 0):
                raise ValueError('Provided visit file {} has unknown formatting; cannot parse.'.format(visit_file))

        # Read in as astropy table
        names = ['Order', 'Star IDs', 'FGS', 'Corrected RA', 'Corrected Dec',
                 'Probability', 'ID V2', 'ID V3', 'ID X', 'ID Y', 'ID PA @ star',
                 'FGS Magnitude' ,'FGS Mag Uncert', 'Count Rate', 'Count Rate Uncert']
        selected_guide_stars = asc.read(visit_file, data_start=i_start - 1, names=names)

        # Verify only one set of GS parameters are provided and extract them
        GS_params = []
        for col in ['Corrected RA', 'Corrected Dec', 'ID PA @ star', 'FGS']:
            values = set(selected_guide_stars[col])
            if len(values) > 1:
                raise ValueError('Cannot parse {} from input visit file; too many values provided: {}'.format(col, list(values)))
            param = list(values)[0]
            GS_params.append(param)
        RA, Dec, PA, fgsNum = GS_params

        # Update class attributes
        self.RA = RA
        self.Dec = Dec
        self.PA = PA
        self.fgsNum = fgsNum

        return RA, Dec, PA, fgsNum

    def plot_segments(self):
        # Plot segments in V2/V3 frame
        plt.figure(1)
        plt.clf()
        plt.plot(self.V2SegArray, self.V3SegArray, 'b*')
        plt.plot(self.V2SegArray.mean(), self.V3SegArray.mean(), 'ro')
        plt.grid(True)
        plt.axis([40.0, -40.0, -40.0, 40.0])  # V2 to the left
        plt.title('Segments')
        plt.xlabel('<-- Delta V2')
        plt.ylabel('Delta V3 -->')
        for s in range(len(self.V2SegArray)):
            plt.text(self.V2SegArray[s], self.V3SegArray[s], self.SegIDArray[s])
        plt.savefig(os.path.join(self.out_dir, self.root + '_V2V3segments.png'))

        # Plot calculate segments' RA and Dec
        plt.figure(2)
        plt.clf()
        plt.plot(self.SegRA, self.SegDec, 'b*')
        RAmean = self.SegRA.mean()
        Decmean = self.SegDec.mean()
        plt.plot(RAmean, Decmean, 'ro')
        segN = int(self.segNum)
        if segN > 0:
            plt.plot(self.SegRA[segN - 1], self.SegDec[segN - 1], 'mx', markersize=12)
        for i in range(len(self.V2SegArray)):
            plt.text(self.SegRA[i], self.SegDec[i], str(i + 1))
        plt.grid(True)
        plt.title('Segment RA and Dec')
        plt.xlabel('RA')
        plt.ylabel('Dec')
        plt.gca().invert_xaxis()
        plt.gca().ticklabel_format(useOffset=False)
        plt.savefig(os.path.join(self.out_dir, self.root + '_RADecsegments.png'))

    def check_guidestar_params(self):
        """
        Ensure all guidestar parameters (RA, Dec, PA, and boresight offset) fall
        within appropriate ranges.
        """
        V2B = self.V2Boff
        V3B = self.V3Boff
        gsRA = self.RA
        gsDec = self.Dec
        gsPA = self.PA

        # ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ BORESIGHT OFFSET ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

        # Guide star information
        msg = ["OK", "Boresight parameter conversion error",
               "Boresight parameter out of range"]

        errcode = self.checkout(V2B, -10.0, 10.0)
        if errcode != 0:
            error = msg[errcode]
            raise ValueError(error)
            return error
        else:
            V2B = float(V2B)

        errcode = self.checkout(V3B, -10.0, 10.0)
        if errcode != 0:
            error = msg[errcode]
            raise ValueError(error)
            return error
        else:
            V3B = float(V3B)

        # ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ RA, DEC, AND POSITION ANGLE ~ ~ ~ ~ ~ ~

        msg = ["OK", "Guide Star parameter conversion error",
               "Guide Star parameter out of range"]

        errcode = self.checkout(gsRA, 0.0, 360.0)
        if errcode != 0:
            error = msg[errcode]
            raise ValueError(error)
            return error
        else:
            gsRA = float(gsRA)

        errcode = self.checkout(gsDec, -90.0, 90.0)
        error = msg[errcode]
        if errcode != 0:
            error = msg[errcode]
            raise ValueError(error)
            return error
        else:
            gsDec = float(gsDec)

        errcode = self.checkout(gsPA, -180.0, 180.0)
        error = msg[errcode]
        if errcode != 0:
            raise ValueError(error)
            return
        else:
            gsPA = float(gsPA)

        return V2B, V3B, gsRA, gsDec, gsPA

    def checkout(self, str, low, high):
        """Test conversion from string to float.
        if float conversion works, test range
        return errcode 0 for OK, 1 for conversion error, 2 for range error"""

        try:
            x = float(str)
            if low <= x <= high:
                errcode = 0
            else:
                errcode = 2
        except ValueError:
            errcode = 1
        return errcode

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# MAIN FUNCTION
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


def run_tool(segment_infile, program_id=0, observation_num=0, visit_num=0, root=None,
             GUI=False, GS_params_dict=None, selected_segs=None, vss_infile=None,
             out_dir=None, data=None, masterGUIapp=None, refonly=False,
             ct_uncert_fctr=0.9, countrate_factor=None):
    # if not GS_params_dict and not GUI:
    #     GS_params_dict = {'V2Boff': 0.1,  # V2 boresight offset
    #                       'V3Boff': 0.2,  # V3 boresight offset
    #                       'fgsNum': 1,  # guider number
    #                       'RA': 30.,  # RA of guide star
    #                       'Dec': 50.,  # Dec of guide star
    #                       'PA': 2.,  # position angle of guide star
    #                       'segNum': 0}  # selected segment to guide on

    root = utils.make_root(root, segment_infile)
    utils.create_logger_from_yaml(__name__, root=root, level='DEBUG')

    try:
        if GUI:
            # Parse ALLpsfs.txt for locations of segments
            all_segment_locations = asc.read(segment_infile)
            x = all_segment_locations['x']
            y = all_segment_locations['y']
            coords = [(x_i, y_i) for x_i, y_i in zip(x, y)]

            # Find the minimum distance between PSFs
            if len(coords) < 2:
                # For cases where we only have star, we assume that we are sufficiently
                #isolated from other stars, but also that the guide star's PSF may be
                #distorted enough that it might appear quite large on the detector
                dist = 20
            else:
                dist = np.floor(np.min(utils.find_dist_between_points(coords))) - 1.

            # Run the GUI to select guide and reference stars
            inds, segNum = SegmentGuidingGUI.run_SelectSegmentOverride(data, x, y, dist, masterGUIapp=masterGUIapp)
            LOGGER.info('Segment Guiding: {} segment override commands generated with segNum = {}'.format(len(inds), segNum))
            GS_params_dict['segNum'] = segNum

            # Turn index list into selected segments file
            selected_segs = np.array(inds)

        # Set up guiding calculator object
        sg = SegmentGuidingCalculator(segment_infile, program_id, observation_num,
                                      visit_num, root=root, GUI=GUI,
                                      GS_params_dict=GS_params_dict,
                                      selected_segs=selected_segs,
                                      vss_infile=vss_infile, out_dir=out_dir,
                                      refonly=refonly, ct_uncert_fctr=ct_uncert_fctr,
                                      countrate_factor=countrate_factor)
        sg.ChosenSeg()
        sg.Calculate()

    except Exception as e:
        LOGGER.exception(e)
        raise
