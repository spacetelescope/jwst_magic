"""Generate segment guiding override files.

The JWST MAGIC Segment Guiding Tool is used during early wavefront
commissioning, when the mirror segments are still unstacked and/or unphased.
The tool takes a list of segment locations and guide star parameters,
using them to calculates the effective RA and Dec where each segment
appears in the sky. The user then specifies which segments to select as
the "guide star" and as the "reference stars". Finally, the tool generates
either a segment override file (SOF) or a photometry override file (POF),
in the form of gs-override*.txt file that the Visit Scheduling System
(VSS) will use to generate a visit.

Authors
-------
    - Colin Cox (original creator, May 2017)
    - Lauren Chambers (modifications in 2018)
    - Keira Brooks (modifications in 2018)

Use
---
    To generate a segment override file (SOF), the segment_guiding module
    can be used in the following way:
    ::
        from jwst_magic.segment_guiding import segment_guiding
        segment_guiding.generate_segment_override_file(
            segment_infile, guider, program_id, observation_num, visit_num,
            root=None, out_dir=None, selected_segs=None,
            click_to_select_GUI=True, data=None, guide_star_params_dict=None,
            threshold_factor=0.9, parameter_dialog=True, oss_factor=0.6,
            masterGUIapp=None):


    Or to generate a photometry override file (POF):
    ::
         segment_guiding.generate_photometry_override_file(
            root, program_id, observation_num, visit_num,
            countrate_factor, out_dir=None):
"""

# Standard Library Imports
import getpass
import logging
import os
import time


# Third Party Imports
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.io import ascii as asc
from astropy.table import Table
import matplotlib
jenkins = 'jenkins' in os.getcwd()
if matplotlib.get_backend() != 'Qt5Agg' and not jenkins:
    matplotlib.use("Qt5Agg")
import matplotlib.path as mpltPath
import matplotlib.pyplot as plt
import numpy as np
from PyQt5 import uic
from PyQt5.QtWidgets import QDialog
import pysiaf
from pysiaf.utils import rotations

# Local Imports
from .. import utils
from jwst_magic._utils import coordinate_transforms
from ..segment_guiding import SegmentGuidingGUI

# Establish segment guiding files directory
__location__ = os.path.realpath(os.path.join(os.getcwd(), os.path.dirname(__file__)))
PACKAGE_PATH = os.path.split(__location__)[0]
OUT_PATH = os.path.split(PACKAGE_PATH)[0]  # Location of out/ and logs/ directory

# Open the SIAF with pysiaf
FGS_SIAF = pysiaf.Siaf('FGS')

# Start logger
LOGGER = logging.getLogger(__name__)


class SegmentGuidingCalculator:
    def __init__(self, override_type, program_id, observation_num, visit_num,
                 root, out_dir, segment_infile=None, guide_star_params_dict=None,
                 selected_segs=None, threshold_factor=0.9, countrate_factor=None,
                 oss_factor=0.6):
        """Initialize the segment guiding calculator class.

        Parameters
        ----------
        override_type: str
            What kind of file to generate. Options are "SOF" (segment override
            file) or "POF" (photometry override file)
        program_id : int
            APT program number
        observation_num : int
            Observation number
        visit_num : int
            Visit number
        root : str
            Name used to generate output folder and output filenames.
        out_dir : str
            Location of out/ directory.

        segment_infile : str, optional
            Filepath to all_found_psfs*.txt file with list of all segment locations
            and countrates
            Used for SOF Generation
        guide_star_params_dict : dict, optional
            Dictionary containing guide star parameters, for example:
            {'v2_boff': 0.1,  # boresight offset in V2 (arcsec)
             'v3_boff': 0.2,  # boresight offset in V3 (arcsec)
             'fgs_num': 1,  # guider number
             'ra': '271d 05m 14.85s',  # RA of guide star (Allowed range: 0 - 360 degrees; flexible format)
             'dec': '-29:31:08.9',  # Dec of guide star (Allowed range: 0 - 360 degrees; flexible format)
             'pa': 2.,  # position angle of guide star (Allowed range: 0 - 360 degrees)
             'seg_num': 0}  # selected segment to guide on
             Used for SOF Generation
        selected_segs : str, optional
            Filepath to guiding_selections*.txt file with list of locations and
            countrates for the selected segments (guide and reference stars)
            Used for SOF Generation
        threshold_factor : float, optional
            The factor by which countrates are multiplied to determine
            the countrate uncertainty
            Used for SOF Generation
        countrate_factor : float, optional
            The factor by which countrates are multiplied by to simulate
            diffuse PSFs (e.g. in MIMF)
            Used for POF Generation
        """

        # Initialize parameters into attributes
        self.override_type = override_type
        self.program_id = int(program_id)
        self.observation_num = int(observation_num)
        self.visit_num = int(visit_num)
        self.root = root
        self.out_dir = out_dir
        self.threshold_factor = threshold_factor
        self.countrate_factor = countrate_factor
        self.oss_factor = oss_factor

        # Initialize other attributes
        # Will the override file be written out using the 'ref-only' syntax? Should be yes.
        self._refonly = True

        # Ensure the output directory exists
        utils.ensure_dir_exists(self.out_dir)

        # Set up to do segment override calculations
        if self.override_type == "SOF":
            self.get_gs_params(guide_star_params_dict)
            self.get_guider_aperture()
            self.parse_infile(segment_infile)
            self.get_selected_segs(selected_segs)

    def get_center_pointing(self):
        """Determine the V2/V3 position of the chosen segment.

        Check that the user-provided segment ID number is valid (0 to
        18). Calculate the central V2/V3 point of either the provided
        segment or the center of the segment array.
        """

        # Ensure the provided segment ID is valid
        segment_max = len(self.v2_seg_array)
        if (self.seg_num < 0) or (self.seg_num > segment_max):
            msg = 'Segment number {} out of range (0, {})'.format(self.seg_num, segment_max)
            raise ValueError(msg)

        # Determine the central V2/V3 point from the given segment ID

        # FIXME: is this wrong? adding the V2/V3 reference twice...?

        # If a specific segment was provided, set the V2/V3 ref point to be
        # that segment's location
        if self.seg_num > 0:
            self.v2_seg_n = self.v2_seg_array[self.seg_num - 1]
            self.v3_seg_n = self.v3_seg_array[self.seg_num - 1]
        # Otherwise, if the input segment ID was 0, set the V2/V3 ref point to
        # be the mean of all segments' locations
        else:
            self.v2_seg_n = self.v2_seg_array.mean()
            self.v3_seg_n = self.v3_seg_array.mean()

        # Calculate aim
        self.v2_ref = float(self.v2_ref)   # obtained from FGS setup
        self.v3_ref = float(self.v3_ref)
        dv2_aim = self.v2_seg_n
        dv3_aim = self.v3_seg_n
        self.v2_aim = self.v2_ref + dv2_aim
        self.v3_aim = self.v3_ref + dv3_aim

        # Convert to Ideal coordinates
        self.x_idl_aim, self.y_idl_aim = self.fgs_siaf_aperture.tel_to_idl(self.v2_aim, self.v3_aim)

    def get_guider_aperture(self):
        """Extract needed parameters from the SIAF file for the given FGS.

        Taking the current guider number, extracts V2Ref, V3Ref,
        V3IdlYAngle, and VIdlParity from the respective aperture in the
        FGS SIAF.
        """
        # Ensure the guider number is valid
        if str(self.fgs_num) not in ['1', '2']:
            raise ValueError('Invalid guider number: "{}"'.format(self.fgs_num))

        # Construct the aperture name
        det = 'FGS' + str(self.fgs_num) + '_FULL_OSS'

        # Read SIAF for appropriate guider aperture with pysiaf
        self.fgs_siaf_aperture = FGS_SIAF[det]
        self.v2_ref = self.fgs_siaf_aperture.V2Ref
        self.v3_ref = self.fgs_siaf_aperture.V3Ref
        self.v3_idl_yangle = self.fgs_siaf_aperture.V3IdlYAngle
        self.v_idl_parity = self.fgs_siaf_aperture.VIdlParity

    def calculate_effective_ra_dec(self):
        """Calculate the effective RAs and Decs for each segment.
        """
        self.n_segments = len(self.seg_id_array)

        # Get the attitude matrix
        attitude = rotations.attitude(self.v2_aim + self.v2_boff,
                                      self.v3_aim + self.v3_boff,
                                      self.ra, self.dec, self.pa)

        # Get RA and Dec for each segment.
        self.seg_ra = np.zeros(self.n_segments)
        self.seg_dec = np.zeros(self.n_segments)
        for i in range(self.n_segments):
            v2 = self.v2_ref + self.v2_seg_array[i]
            v3 = self.v3_ref + self.v3_seg_array[i]
            self.seg_ra[i], self.seg_dec[i] = rotations.pointing(attitude, v2, v3,
                                                                 positive_ra=True)

        # Convert V2/V3 coordinates to ideal coordinates and detector frame coordinates
        idl_coords = self.fgs_siaf_aperture.tel_to_idl(self.v2_seg_array + self.v2_ref,
                                                       self.v3_seg_array + self.v3_ref)
        self.x_idl_segs, self.y_idl_segs = idl_coords
        self.x_det_segs, self.y_det_segs = self.fgs_siaf_aperture.idl_to_det(self.x_idl_segs,
                                                                             self.y_idl_segs)

        # Check to make sure all the computed segment locations are within
        # the needed FOV
        self.check_segments_inside_fov(attitude)

    def write_override_file(self, verbose=True):
        """Write the segment guiding override file: {out_dir}/out/{root}/
        gs-override_{program_id}_{observation_num}_{visit_num}.txt

        Parameters
        ----------
        verbose : bool, optional
            Log results of calculations and file content
        """
        # Define path and name of output override file
        out_file = 'gs-override-{}_{}_{}.txt'.format(self.program_id, self.observation_num,
                                                     self.visit_num)
        out_file = os.path.join(self.out_dir, out_file)

        # Print summary of input data (guide star RA, Dec, and PA, etc...)
        if verbose and self.override_type == "SOF":
            # Print guide star and boresight parameters
            summary_output = """Guide Star Parameters
                Aperture FGS: {0}
                V2/V3 Refs: ({1:10.4f}, {2:10.4f}) arc-sec
                Guiding segment number: {3}
                V2/V3 Boresight offset: ({4}, {5}) arc-sec
                Guide star RA & Dec: ({6}, {7}) degrees
                Position angle: {8} degrees""".\
                format(self.fgs_num, self.v2_ref, self.v3_ref, self.seg_num,
                       self.v2_boff, self.v3_boff, self.ra, self.dec, self.pa)
            LOGGER.info('Segment Guiding: ' + summary_output)

            all_segments = 'All Segment Locations'
            all_segments += '\n                Segment     dV2    dV3    xIdl' +\
                            '   yIdl     RA         Dec         xDet     yDet'
            for p in range(self.n_segments):
                all_segments += ('\n                %5s    %6.2f %6.2f  %6.2f %6.2f  %10.6f %10.6f  %8.2f %8.2f'
                                 % (self.seg_id_array[p], self.v2_seg_array[p],
                                    self.v3_seg_array[p], self.x_idl_segs[p],
                                    self.y_idl_segs[p], self.seg_ra[p],
                                    self.seg_dec[p], self.x_det_segs[p], self.x_det_segs[p]))
            LOGGER.info('Segment Guiding: ' + all_segments)

        # Write out override file with RA/Decs of selected segments
        with open(out_file, 'w') as f:
            # Determine whether to include a multiplicative countrate factor
            countrate_qualifier = ' -count_rate_factor={:.3f}'.\
                format(self.countrate_factor) if self.countrate_factor else ''
            out_string = 'sts -gs_select {:4d}:{}:{}{}'.\
                format(self.program_id, self.observation_num,
                       self.visit_num, countrate_qualifier)
            if self.override_type == "SOF":
                # Determine which segments have been selected
                orientations = list(self.selected_segment_ids)
                guide_segments = [s[0] for s in orientations]
                all_selected_segs = list(set(np.concatenate(orientations)))

                if self._refonly:
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

                # If countrates were included in the input file, use them!
                # The OSS factor is multiplied by OSS to account for the 3x3 countrates
                # So we need to divide it out here and for every segment override file
                rate = self.countrate_array / self.oss_factor
                uncertainty = self.countrate_array * self.threshold_factor / self.oss_factor

                # Write the commands for each orientation
                for i_o, orientation in enumerate(orientations):
                    guide_seg_id = orientation[0]

                    label = 'star'
                    seg = i_o + 1
                    # If implementing "ref_only" labels, determine if star is guide
                    # or reference star, and alter label and ID accordingly
                    if self._refonly:
                        if i_o >= n_guide_segments:
                            label = 'ref_only'
                            seg = i_o + 1 - n_guide_segments

                    # Format segment properties (ID, RA, Dec, countrate, uncertainty)
                    star_string = ' -%s%d = %d, %.6f, %.6f, %.1f, %.1f' % (
                        label, seg, guide_seg_id + 1, self.seg_ra[guide_seg_id],
                        self.seg_dec[guide_seg_id], rate[guide_seg_id],
                        uncertainty[guide_seg_id])

                    if not self._refonly or (self._refonly and label == 'star'):
                        # Add list of segment IDs for all reference stars
                        for ref_seg_id in orientation:
                            if ref_seg_id != guide_seg_id:
                                star_string += ', %d' % (ref_seg_id + 1)

                    out_string += star_string

                # Write out the override report
                self.write_override_report(orientations, n_guide_segments)

            f.write(out_string)

            if verbose:
                LOGGER.info('Segment Guiding: Guide Star Override: ' +
                            out_string.replace('-star', '\n                -star').
                            replace('-ref_only', '\n                -ref_only'))
                LOGGER.info('Segment Guiding: Saved override command to {}'.
                            format(out_file))



    def write_override_report(self, orientations, n_guide_segments):
        # Define path and name of output override report
        out_file = 'gs-override-{}_{}_{}_REPORT.txt'.format(self.program_id, self.observation_num,
                                                     self.visit_num)
        out_file = os.path.join(self.out_dir, out_file)

        username = getpass.getuser()

        with open(out_file, 'w') as f:
            f.write('Guide Star Override Report\n')
            f.write('Generated on {} at {} by {}\n'.format(time.strftime("%Y/%m/%d"), time.strftime("%H:%M:%S"), username))
            f.write('\n')
            f.write('{:14s}: {:s}'.format('Program ID', str(self.program_id)) + '\n')
            f.write('{:14s}: {:s}'.format('Observation #', str(self.observation_num)) + '\n')
            f.write('{:14s}: {:s}'.format('Visit #', str(self.visit_num)) + '\n')
            f.write('{:14s}: {:d}'.format('FGS Detector', self.fgs_num) + '\n')
            f.write('{:14s}: {:f}'.format('Guide Star RA', self.ra) + '\n')
            f.write('{:14s}: {:f}'.format('Guide Star Dec', self.dec) + '\n')
            f.write('{:14s}: {:f}'.format('V3 PA @ GS', self.pa) + '\n')
            f.write('\n')
            columns = 'Star Name, MAGIC ID, RA, Dec, Ideal X, Ideal Y, OSS Ideal X, OSS Ideal Y, Detector X, Detector Y'.split(', ')
            header_string_to_format = '{:^13s}| '* len(columns)
            header_string = header_string_to_format[:-2].format(*columns) + '\n'
            f.write(header_string)
            f.write('-'*len(header_string) + '\n')

            for i_o, orientation in enumerate(orientations):
                guide_seg_id = orientation[0]

                label = 'star'
                seg = i_o + 1
                if self._refonly:
                    if i_o >= n_guide_segments:
                        label = 'ref_only'
                        seg = i_o + 1 - n_guide_segments

                values = [label + str(seg), int(guide_seg_id + 1), self.seg_ra[guide_seg_id], self.seg_dec[guide_seg_id],
                           self.x_idl_segs[guide_seg_id], self.y_idl_segs[guide_seg_id],
                           -self.x_idl_segs[guide_seg_id], self.y_idl_segs[guide_seg_id],
                           self.x_det_segs[guide_seg_id], self.y_det_segs[guide_seg_id]]

                row_string_to_format = '{:<13s}| ' + '{:<13d}| ' +  '{:<13f}| '* (len(values) - 2)
                row_string = row_string_to_format[:-2].format(*values) + '\n'
                f.write(row_string)

    def check_segments_inside_fov(self, attitude):
        """Check to make sure that the calculated RA and Dec of each
        segment is within the field of view of the given FGS.

        Parameters
        ----------
        attitude : 3 x 3 numpy array
            Attitude matrix generated by pysiaf.utils.rotations.attitude
        """
        # Check to make sure no segments are off the detector, pixel-wise
        for x, y, i_seg in zip(self.x_det_segs, self.y_det_segs, self.seg_id_array):
            if x < 0.5 or x > 2048.5:
                LOGGER.warning('Segment Guiding: %8s off detector in X direction' % i_seg)
            if y < 0.5 or y > 2048.5:
                LOGGER.warning('Segment Guiding: %8s off detector in Y direction' % i_seg)

        # Get the vertices of the given FGS detector
        vertices = [(self.fgs_siaf_aperture.XIdlVert1, self.fgs_siaf_aperture.YIdlVert1),
                    (self.fgs_siaf_aperture.XIdlVert2, self.fgs_siaf_aperture.YIdlVert2),
                    (self.fgs_siaf_aperture.XIdlVert3, self.fgs_siaf_aperture.YIdlVert3),
                    (self.fgs_siaf_aperture.XIdlVert4, self.fgs_siaf_aperture.YIdlVert4)]

        # Convert the vertices from ideal coordinates to RA & Dec
        vertices_sky = []
        for v_idlx, v_idly in vertices:
            v_v2, v_v3 = self.fgs_siaf_aperture.idl_to_tel(v_idlx, v_idly)
            v_sky = rotations.pointing(attitude, v_v2, v_v3)
            vertices_sky.append(v_sky)

        # Create a matplotlib Path that encompasses the entire detector
        fov_path = mpltPath.Path(vertices_sky)

        # Determine if every segment RA and Dec is within that Path (i.e.
        # within the guider FOV)
        seg_pointings = [(seg_ra, seg_dec) for seg_ra, seg_dec in zip(self.seg_ra, self.seg_dec)]
        segs_in_fov = fov_path.contains_points(seg_pointings)
        if not segs_in_fov.all():
            segments_outside = np.where(segs_in_fov is False)[0]
            raise ValueError(
                'Incorrect segment guiding calculations. Segment(s) {} is outside of the FGS{} FOV. Cannot generate segment override file that will not fail.'
                .format(segments_outside, self.fgs_num)
            )

        # And because I don't trust anything anymore, straight up check that the
        # segments are within a guider ~FOV of the commanded GS RA/Dec
        gs_pointing = SkyCoord(ra=self.ra * u.degree, dec=self.dec * u.degree)
        fgs_fov_length = 2.3 * u.arcmin
        fgs_radius = np.sqrt(2 * (fgs_fov_length / 2) ** 2)
        for i, p in enumerate(seg_pointings):
            p = SkyCoord(ra=p[0] * u.degree, dec=p[1] * u.degree)
            sep = p.separation(gs_pointing)
            if sep > fgs_radius:
                raise ValueError('Segment {} at RA, Dec = ({}, {}) is outside the FGS{} FOV. Cannot generate segment override file that will not fail.'
                                 .format(i + 1, p.ra, p.dec, self.fgs_num))

    def get_gs_params(self, guide_star_params_dict):
        """Map guide_star_params_dict values and keys to attributes
        """

        self.v2_boff = float(guide_star_params_dict['v2_boff'])
        self.v3_boff = float(guide_star_params_dict['v3_boff'])
        self.fgs_num = int(guide_star_params_dict['fgs_num'])

        ra_value = guide_star_params_dict['ra']
        dec_value = guide_star_params_dict['dec']
        gs_coord = SkyCoord(ra_value, dec_value, unit=(u.deg, u.deg))
        self.ra = gs_coord.ra.degree
        self.dec = gs_coord.dec.degree

        self.pa = float(guide_star_params_dict['pa'])
        self.seg_num = int(guide_star_params_dict['seg_num'])

    def parse_infile(self, segment_infile):
        """Get the segment positions and count rates from a file.

        Parameters
        ----------
        segment_infile : str
            File containing segment locations and count rates

        Raises
        ------
        TypeError
            Incompatible file type provided as segment_infile
        """
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
                v2, v3 = coordinate_transforms.Raw2Tel(read_table['x'],
                                                       read_table['y'],
                                                       self.fgs_num)
                segment_coords['V2Seg'], segment_coords['V3Seg'] = v2, v3

            else:
                raise TypeError('Incompatible file type: ', segment_infile)

            # If the countrates are included in the input file, read them!
            if any(['countrate' == c for c in column_names]):
                self.countrate_array = read_table['countrate']

        else:
            raise TypeError('Incompatible file type: ', segment_infile)

        # Define the IDs and coordinates of all segments
        self.seg_id_array = segment_coords['SegID']
        self.v2_seg_array = segment_coords['V2Seg']
        self.v3_seg_array = segment_coords['V3Seg']
        LOGGER.info(
            'Segment Guiding: {} segment coordinates read from {}'.
            format(len(self.seg_id_array), segment_infile)
        )

    def get_selected_segs(self, selected_segs):
        """If a file of selected segments has been provided, get their
        locations and count rates.

        Parameters
        ----------
        selected_segs : str
            File containing locations and count rates of selected segments

        Raises
        ------
        TypeError
            Incompatible file type provided as selected_segs
        """

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

        # Else if they are a guiding_selections*.txt, parse it.
        elif os.path.exists(selected_segs):
            LOGGER.info('Segment Guiding: Guiding on segments selected in {}'.format(selected_segs))
            self.parse_guiding_selections_file(selected_segs)

        # Otherwise, we don't know what it is...
        else:
            raise TypeError(
                'Unrecognized data type passed to selected_segs ({}); must be guiding_selections*.txt path or array of indices.'.
                format(selected_segs)
            )

    def parse_guiding_selections_file(self, selected_segs):
        """Extract the segment positions and count rates from a guiding_selections*.txt

        Parameters
        ----------
        selected_segs : str
            Filepath to guiding_selections*.txt

        Raises
        ------
        TypeError
            Incompatible guiding_selections*.txt file provided as selected_segs
        """
        # Make sure the file is ascii-readable
        n_segs = len(self.seg_id_array)
        try:
            read_selected_segs = asc.read(selected_segs)
            column_names = read_selected_segs.colnames
            LOGGER.info(
                'Segment Guiding: Selected segment coordinates read from {}'.
                format(selected_segs)
            )
        except:
            raise TypeError('Incompatible guiding_selections*.txt type: ', selected_segs)

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
                                                   self.fgs_num)
            selected_segment_coords['V2Seg'], selected_segment_coords['V3Seg'] = v2, v3

        # If the positions aren't V2/V3 or x/y, what are they??
        else:
            raise TypeError('Incompatible guiding_selections*.txt type: ', selected_segs)

        # Match locations of selected segments to IDs of known segments
        selected_segs_ids = []
        for coords in zip(selected_segment_coords['V2Seg', 'V3Seg']):
            for i_seg in range(n_segs):
                if coords[0]['V2Seg'] == self.v2_seg_array[i_seg] and \
                   coords[0]['V3Seg'] == self.v3_seg_array[i_seg]:
                    selected_segs_ids.append(i_seg)
        if len(selected_segs_ids) == 0:
            raise TypeError(
                'Coordinates of selected segments file do not match those of '
                'the provided input file.'
            )

        self.selected_segment_ids = [selected_segs_ids]

    def plot_segments(self):
        """Generate and save plots of segments in V2/V3 and RA/Dec.
        """
        # Plot segments in V2/V3 frame
        plt.figure(1)
        plt.clf()
        plt.plot(self.v2_seg_array, self.v3_seg_array, 'b*')
        plt.plot(self.v2_seg_array.mean(), self.v3_seg_array.mean(), 'ro')
        plt.grid(True)
        plt.axis([40.0, -40.0, -40.0, 40.0])  # V2 to the left
        plt.title('Segments')
        plt.xlabel('<-- Delta V2')
        plt.ylabel('Delta V3 -->')
        for s in range(len(self.v2_seg_array)):
            plt.text(self.v2_seg_array[s], self.v3_seg_array[s], self.seg_id_array[s])
        plt.savefig(os.path.join(self.out_dir, self.root + '_V2V3segments.png'))

        # Plot calculate segments' RA and Dec
        plt.figure(2)
        plt.clf()
        plt.plot(self.seg_ra, self.seg_dec, 'b*')
        ra_mean = self.seg_ra.mean()
        dec_mean = self.seg_dec.mean()
        plt.plot(ra_mean, dec_mean, 'ro')
        seg_n = int(self.seg_num)
        if seg_n > 0:
            plt.plot(self.seg_ra[seg_n - 1], self.seg_dec[seg_n - 1], 'mx', markersize=12)
        for i in range(len(self.v2_seg_array)):
            plt.text(self.seg_ra[i], self.seg_dec[i], str(i + 1))
        plt.grid(True)
        plt.title('Segment RA and Dec')
        plt.xlabel('RA')
        plt.ylabel('Dec')
        plt.gca().invert_xaxis()
        plt.gca().ticklabel_format(useOffset=False)
        plt.savefig(os.path.join(self.out_dir, self.root + '_RADecsegments.png'))

    def check_guidestar_params(self, override_type):
        """Ensure all guidestar parameters (RA, Dec, PA, and boresight
        offset) fall within appropriate ranges.

        Parameters
        ----------
        override_type: str
            What kind of file to generate. Options are "SOF" (segment override
            file) or "POF" (photometry override file)
        """
        if override_type == "SOF":
            # ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ BORESIGHT OFFSET ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
            # Guide star information
            # msg = ["OK", "Boresight parameter conversion error",
            #        "Boresight parameter out of range"]
            #
            # if self.v2_boff:
            #     errcode = self.checkout(self.v2_boff, -10.0, 10.0)
            #     if errcode != 0:
            #         error = msg[errcode]
            #         raise ValueError(error)
            #     else:
            #         self.v2_boff = float(self.v2_boff)
            #
            # if self.v3_boff:
            #     errcode = self.checkout(self.v3_boff, -10.0, 10.0)
            #     if errcode != 0:
            #         error = msg[errcode]
            #         raise ValueError(error)
            #     else:
            #         self.v3_boff = float(self.v3_boff)

            # ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ RA, DEC, AND POSITION ANGLE ~ ~ ~ ~ ~ ~
            # These values can be set to None, otherwise, have to be in range
            msg = ["OK", "Guide Star parameter conversion error",
                   "Guide Star parameter out of range"]

            if self.ra:
                errcode = self.checkout(self.ra, 0.0, 360.0)
                if errcode != 0:
                    error = "{} for RA. Expecting between 0.0 and 360.0 degrees.".format(msg[errcode])
                    raise ValueError(error)
                else:
                    self.ra = float(self.ra)

            if self.dec:
                errcode = self.checkout(self.dec, -90.0, 90.0)
                if errcode != 0:
                    error = "{} for DEC. Expecting between -90.0 and 90.0 degrees.".format(msg[errcode])
                    raise ValueError(error)
                else:
                    self.dec = float(self.dec)

            if self.pa:
                errcode = self.checkout(self.pa, 0.0, 360.0)
                error = "{} for POSITION ANGLE. Expecting between 0.0 and 360.0 degrees.".format(msg[errcode])
                if errcode != 0:
                    raise ValueError(error)
                else:
                    self.pa = float(self.pa)

        elif override_type == "POF":
            # ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ COUNT RATE FACTOR ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
            # Count rate factor has to be between 0 and 1 (required by VSS)
            msg = ["OK", "Countrate factor conversion error",
                   "Countrate factor out of range"]
            if self.countrate_factor:
                errcode = self.checkout(self.countrate_factor, 0.0, 1.0)
                if errcode != 0:
                    error = "{} for count_rate_factor. Expecting between 0.0 and 1.0.".format(msg[errcode])
                    raise ValueError(error)

    @staticmethod
    def checkout(value, low, high):
        """Test conversion from string to float. If float conversion
        works, test range return errcode 0 for OK, 1 for conversion
        error, 2 for range error
        """

        try:
            x = float(value)
            if low <= x <= high:
                errcode = 0
            else:
                errcode = 2
        except ValueError:
            errcode = 1
        return errcode


def generate_segment_override_file(segment_infile, guider,
                                   program_id, observation_num, visit_num,
                                   root=None, out_dir=None, selected_segs=None,
                                   click_to_select_gui=True,
                                   data=None, guide_star_params_dict=None,
                                   threshold_factor=0.9, parameter_dialog=True,
                                   oss_factor=0.6, master_gui_app=None):
    """Run the segment guiding tool to select guide and reference stars and
    generate a segment guiding override file.

    Parameters
    ----------
    segment_infile : str
        File path to all_found_psfs*.txt file with list of all segment locations
        and countrates
    guider : int
        Which guider is being used: 1 or 2
    program_id : int
        APT program number
    observation_num : int
        Observation number
    visit_num : int
        Visit number

    root : str, optional
        Name used to generate output folder and output filenames. If not
        specified, will be derived from the segment_infile name.
    out_dir : str, optional
        Location of out/ directory. If not specified, will be placed
        within the repository: .../tools/fgs_commissioning/out/
    selected_segs : str, optional
        File path to guiding_selections*.txt file with list of locations and
        countrates for the selected segments (guide and reference stars)
        Required if click_to_select_GUI=False
    click_to_select_gui : bool, optional
        Will the tool use the segment guiding GUI?
        Required if selected_segs=None
    data : 2-D numpy array, optional
        Image that will be displayed in the click-to-select GUI
        Required if click_to_select_GUI=True
    guide_star_params_dict : dict, optional
        Dictionary containing guide star parameters, for example:
            {'v2_boff': 0.1,  # boresight offset in V2 (arcsec)
             'v3_boff': 0.2,  # boresight offset in V3 (arcsec)
             'fgs_num': 1,  # guider number
             'ra': 30.,  # RA of guide star (Allowed range: 0 - 360 degrees)
             'dec': 50.,  # Dec of guide star (Allowed range: 0 - 360 degrees)
             'pa': 2.,  # position angle of guide star (Allowed range: 0 - 360 degrees)
             'seg_num': 0}  # selected segment to guide on
        Required if parameter_dialog=False
    threshold_factor : float, optional
        The factor by which countrates are multiplied to determine
        the countrate uncertainty
    parameter_dialog : bool, optional
        Prompt the user to enter parameters (countrate factors, APT
        numbers, RA, Dec, PA, and boresight offset) from a dialog box
        rather than manually providing arguments.
        Required if guide_star_params_dict=None
    oss_factor : float, optional
        The factor that OSS applies to the 3x3 box countrate in order to represent
        the full number of countrate that we care about.
    master_gui_app : qApplication, optional
        qApplication instance of parent GUI
    """

    root = utils.make_root(root, segment_infile)
    utils.create_logger_from_yaml(__name__, root=root, level='DEBUG')
    out_dir = utils.make_out_dir(out_dir, OUT_PATH, root)

    try:
        # Get the guide star parameters
        if parameter_dialog:
            SOF_parameter_dialog = SegmentGuidingDialog(
                "SOF", guider, program_id, observation_num, visit_num
            )
            accepted = SOF_parameter_dialog.exec()
            params = SOF_parameter_dialog.get_dialog_parameters() if accepted else None

            if params is not None:
                guide_star_params_dict, program_id, observation_num, visit_num, \
                    threshold_factor, _ = params
            else:
                LOGGER.warning('Segment Guiding: SOF creation cancelled.')
                return

        elif guide_star_params_dict is None:
            raise ValueError(
                'In order to run the segment guiding tool with '
                '`parameter_dialog=False`, you must provide a dictionary '
                'of guide star parameters to the `guide_star_params_dict` '
                'argument in `segment_guiding.generate_segment_override_file()`.'
            )

        # Determine which segments are the guide and reference segments
        if click_to_select_gui:
            if data is None:
                raise ValueError(
                    'In order to run the segment guiding tool with '
                    '`click_to_select_GUI=True`, you must provide a 2-D '
                    'numpy array of the image to the `data` argument in '
                    '`segment_guiding.generate_segment_override_file()`.'
                )

            guide_star_params_dict, selected_segs = _click_to_select_segments(
                segment_infile, data, guide_star_params_dict,
                master_gui_app, selected_segs=selected_segs
            )
        elif selected_segs is None:
            raise ValueError(
                'In order to run the segment guiding tool with '
                '`click_to_select_GUI=False`, you must provide a file '
                'specifying the locations and count rates of the guide '
                'and reference stars as the `selected_segs` argument in'
                '`segment_guiding.generate_segment_override_file()`.'
            )

        # Set up guiding calculator object
        sg = SegmentGuidingCalculator(
            "SOF", program_id, observation_num, visit_num, root, out_dir,
            segment_infile=segment_infile,
            guide_star_params_dict=guide_star_params_dict,
            selected_segs=selected_segs, threshold_factor=threshold_factor,
            oss_factor=oss_factor
        )
        # Verify all guidestar parameters are valid
        sg.check_guidestar_params("SOF")

        # Determine the V2/V3 of the pointing center
        sg.get_center_pointing()

        # Write a SOF
        sg.calculate_effective_ra_dec()
        sg.write_override_file()  # Print and save final output
        sg.plot_segments()  # Save .pngs of plots

    except Exception as e:
        LOGGER.exception(e)
        raise


def generate_photometry_override_file(root, program_id, observation_num, visit_num,
                                      countrate_factor=None, out_dir=None,
                                      parameter_dialog=True):
    """Generate a photometry override file (used for commissioning activities
    where the PSFs are stacked but unphased, like during MIMF).

    Parameters
    ----------
    root : str
        Name used to generate output folder and output filenames. If not
        specified, will be derived from the segment_infile name
    program_id : int
        APT program number
    observation_num : int
        Observation number
    visit_num : int
        Visit number

    countrate_factor : float, optional
        The factor by which countrates are multiplied by to simulate
        diffuse PSFs (e.g. in MIMF).
        Required if parameter_dialog=False
    out_dir : str, optional
        Location of out/ directory. If not specified, will be placed
        within the repository: .../tools/fgs_commissioning/out/
    parameter_dialog : bool, optional
        Prompt the user to enter parameters (countrate factors, APT
        numbers) from a dialog box rather than manually providing arguments.
        Required if countrate_factor=None
    """

    utils.create_logger_from_yaml(__name__, root=root, level='DEBUG')
    out_dir = utils.make_out_dir(out_dir, OUT_PATH, root)

    try:
        # Get the program parameters and countrate factor
        if parameter_dialog:
            POF_parameter_dialog = SegmentGuidingDialog(
                "POF", None, program_id, observation_num, visit_num
            )
            accepted = POF_parameter_dialog.exec()
            params = POF_parameter_dialog.get_dialog_parameters() if accepted else None

            _, program_id, observation_num, visit_num, _, countrate_factor = params

        elif countrate_factor is None:
            raise ValueError(
                'In order to run the segment guiding tool with '
                '`parameter_dialog=False`, you must provide a countrate factor '
                'between 0 and 1 to the `countrate_factor` argument in '
                '`segment_guiding.generate_photometry_override_file()`.'
            )

        # Set up guiding calculator object
        sg = SegmentGuidingCalculator(
            "POF", program_id, observation_num, visit_num, root, out_dir,
            countrate_factor=countrate_factor
        )
        # Verify all guidestar parameters are valid
        sg.check_guidestar_params("POF")

        # Write the POF
        sg.write_override_file()

    except Exception as e:
        LOGGER.exception(e)
        raise

class SegmentGuidingDialog(QDialog):
    """Define a dialog window to prompt user for guide star parameters
    and other parameters needed to generate the override file.

    Parameters
    ----------
    override_type: str
            What kind of file to generate. Options are "SOF" (segment override
            file) or "POF" (photometry override file)
    guider : int
        Guider number (1 or 2)
    program_id : int
        APT program number
    observation_num : int
        Observation number
    visit_num : int
        Visit number

    Returns
    -------
    tup
        Tuple containing the following arguments: (guide_star_params_dict,
        program_id, observation_num, visit_num, threshold_factor,
        countrate_factor)
    """
    def __init__(self, override_type, guider, program_id, observation_num, visit_num):
        # Initialize attributes
        self.override_type = override_type
        self.guider = guider
        self.program_id = program_id
        self.observation_num = observation_num
        self.visit_num = visit_num

        # Initialize dialog object
        QDialog.__init__(self, modal=True)

        # Import .ui file
        if override_type == "SOF":
            uic.loadUi(os.path.join(__location__, 'segmentOverrideFileDialog.ui'), self)
        elif override_type == "POF":
            uic.loadUi(os.path.join(__location__, 'photometryOverrideFileDialog.ui'), self)

        # Set defaults from parsed header
        self.lineEdit_programNumber.setText(str(program_id))
        self.lineEdit_observationNumber.setText(str(observation_num))
        self.lineEdit_visitNumber.setText(str(visit_num))

    def get_dialog_parameters(self):
        """Parses the user input into the segment guiding dialog box, differentiating
        between input for SOFs and POFs.

        Returns
        -------
        tup
            Tuple containing the following arguments: (guide_star_params_dict,
            program_id, observation_num, visit_num, threshold_factor,
            countrate_factor)
        """

        # Get parameters for dictionary from dialog
        if self.override_type == "SOF":
            # Parse what the boresight offset is
            if self.radioButton_boresightNIRCam.isChecked():
                x_offset = float(self.lineEdit_boresightX.text())
                y_offset = float(self.lineEdit_boresightY.text())
                v2_offset, v3_offset = _convert_nrca3pixel_offset_to_v2v3_offset(x_offset,
                                                                                 y_offset)
                LOGGER.info(
                    'Segment Guiding: Applying boresight offset of {}, {} arcsec (Converted from {}, {} pixels)'.
                        format(v2_offset, v3_offset, x_offset, y_offset)
                )
            else:
                v2_offset = float(self.lineEdit_boresightV2.text())
                v3_offset = float(self.lineEdit_boresightV3.text())
                LOGGER.info(
                    'Segment Guiding: Applying boresight offset of {}, {} arcsec'.
                        format(v2_offset, v3_offset)
                )

            # Parse the RA, Dec, and PA
            ra_value = self.lineEdit_RA.text()
            if self.comboBox_RAUnit.currentText() == 'Degrees':
                ra_unit = u.deg
            elif self.comboBox_RAUnit.currentText() == 'Hours':
                ra_unit = u.hourangle
            dec_value = self.lineEdit_Dec.text()

            gs_coord = SkyCoord(ra_value, dec_value, unit=(ra_unit, u.deg))
            ra = gs_coord.ra.degree
            dec = gs_coord.dec.degree

            pa = float(self.lineEdit_PA.text())

            # Populate the parameter dictionary
            guide_star_params_dict = {
                'v2_boff': v2_offset,
                'v3_boff': v3_offset,
                'fgs_num': self.guider,
                'ra': ra,
                'dec': dec,
                'pa': pa,
                'seg_num': 0
            }

            # Countrate factors
            threshold_factor = float(self.lineEdit_countrateUncertainty.text())
            countrate_factor = None

        elif self.override_type == "POF":
            countrate_factor = float(self.doubleSpinBox_countrateFactor.value())
            threshold_factor = None
            guide_star_params_dict = None

        # Get APT information and other necessary parameters
        program_id = self.lineEdit_programNumber.text()
        observation_num = self.lineEdit_observationNumber.text()
        visit_num = self.lineEdit_visitNumber.text()

        return guide_star_params_dict, program_id, observation_num, visit_num, threshold_factor, countrate_factor


def _click_to_select_segments(segment_infile, data, guide_star_params_dict,
                              master_gui_app, selected_segs=None):
    """Raise the segment guiding GUI and prompt the user to select which segments
    to use as the guide and reference stars.

    Parameters
    ----------
    segment_infile : str
        File path to all_found_psfs*.txt file with list of all segment locations
        and countrates
    data : 2-D numpy array
        Image that will be displayed in the click-to-select GUI
    guide_star_params_dict : dict
        Dictionary containing guide star parameters
    master_gui_app : qApplication or None
        qApplication instance of parent GUI
    selected_segs : str, optional
        File path to guiding_selections*.txt file with list of locations and
        countrates for the pre-selected segments (guide and reference stars)

    Returns
    -------
    guide_star_params_dict : dict
        Dictionary containing guide star parameters with redefined seg_num
        value based on user selection
    selected_segs : str
        Array of indices corresponding to the segments to use as guide and
        reference stars

    """
    if selected_segs is not None:
        # Parse all_found_psfs*.txt for locations of segments
        all_segment_locations = asc.read(segment_infile)
        x = all_segment_locations['x']
        y = all_segment_locations['y']
        coords = [(x_i, y_i) for x_i, y_i in zip(x, y)]

        # Find the minimum distance between PSFs
        if len(coords) < 2:
            # For cases where we only have star, we assume that we are sufficiently
            # isolated from other stars, but also that the guide star's PSF may be
            # distorted enough that it might appear quite large on the detector
            dist = 20
        else:
            dist = np.floor(np.min(utils.find_dist_between_points(coords))) - 1.
    else:
        x = []
        y = []
        dist = 20

    # Run the GUI to select guide and reference stars
    inds, seg_num = SegmentGuidingGUI.run_segment_override_gui(
        data, x, y, dist, selected_segs=selected_segs,
        masterGUIapp=master_gui_app
    )
    LOGGER.info(
        'Segment Guiding: {} segment override commands generated with seg_num = {}'.
        format(len(inds), seg_num)
    )
    guide_star_params_dict['seg_num'] = seg_num

    # Turn index list into selected segments file
    selected_segs = np.array(inds)

    return guide_star_params_dict, selected_segs

def _convert_nrca3pixel_offset_to_v2v3_offset(x_offset, y_offset):
    """Convert a boresight offset from NIRCam A3 pixels to V2/V3 arcsec

    Parameters
    ----------
    x_offset : float
        Boresight offset in NIRCam A3 X pixels
    y_offset : float
        Boresight offset in NIRCam A3 Y pixels

    Returns
    -------
    v2_offset, v3_offset : tup
        Boresight offset in V2/V3 (arcsec)
    """
    # Get pixel scale
    nrc_siaf = pysiaf.Siaf('NIRCam')
    nrca3 = nrc_siaf['NRCA3_FULL_OSS']
    nircam_sw_x_scale = nrca3.XSciScale  # arcsec/pixel
    nircam_sw_y_scale = nrca3.YSciScale  # arcsec/pixel

    # Convert x/y offsets to V2/V3
    v2_offset = x_offset * nircam_sw_x_scale  # arcsec
    v3_offset = y_offset * nircam_sw_y_scale  # arcsec

    return v2_offset, v3_offset
