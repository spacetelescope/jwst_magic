"""Generate segment guiding override files.

The JWST MAGIC Segment Guiding Tool (SGT) is used during early
wavefront commissioning, when the mirror segments are still unstacked.
The tool takes a list of segment locations and guide star parameters,
using them to calculates the effective RA and Dec where each segment
appears in the sky. The user then specifies which segments to select as
the "guide star" and as the "reference stars". Finally, the tool
generates a segment guiding override file (gs-override*.txt) that the
Visit Scheduling System (VSS) will use to generate a visit.

Authors
-------
    - Colin Cox (original creator, May 2017)
    - Lauren Chambers (modifications in 2018)
    - Keira Brooks (modifications in 2018)

Use
---
    The segment_guiding module can be used either with the Segment Guiding GUI:
    ::
        from jwst_magic.segment_guiding import segment_guiding
        segment_guiding.run_tool(segment_infile=segment_infile, guider=guider)

    Or with the segment dialog box:
    ::
        from jwst_magic.segment_guiding import segment_guiding
        segment_guiding.run_tool(program_id=program_id, observation=observation,
                                 visit=visit, parameter_dialog=True)

    Or from dictionary of parameters:
    ::
        from jwst_magic.segment_guiding import segment_guiding
        segment_guiding.run_tool(program_id=program_id, observation=observation,
                                 visit=visit, guide_star_params_dict=guide_star_params_dict,
                                 parameter_dialog=False)

    Optional arguments:
        ``segment_infile`` - filepath to ALLpsfs.txt file with list
            of all segment locations and countrates
        ``guider`` - which guider is being used: 1 or 2
        ``root`` - name used to generate output folder and output
            filenames. If not specified, will be derived from the
            segment_infile name
        ``out_dir`` - location of out/ directory. If not provided,
            the image(s) will be saved within the repository at
            tools/fgs-commissioning/
        ``selected_segs`` - filepath to regfile.txt file with list of
            locations and countrates for the selected segments (guide
            and reference stars). If not provided, will default to
            {out_dir}/out/{root}/{root}_G{guider}_regfile.txt
        ``GUI`` - will the tool use the segment guiding GUI? If not
            specified, set to True.
        ``vss_infile`` - filepath to guide star report provided by VSS
        ``data`` - image that will be displayed in the click-to-select GUI
        ``masterGUIapp`` - qApplication instance of parent GUI
        ``refonly`` - will the override file be written out using the
            'ref-only' syntax?
        ``parameter_dialog`` - prompt the user to enter parameters
            (countrate factors, APT numbers, RA, Dec, PA, and boresight
            offset) from a dialog box rather than manually providing
            arguments. If True, the remaining arguments are not necessary.
        ``program_id`` - APT program number
        ``observation_num`` - observation number
        ``visit_num`` - visit number
        ``guide_star_params_dict`` - dictionary containing guide star
            parameters, for example:
            {'v2_boff': 0.1,  # boresight offset in V2 (arcsec)
             'v3_boff': 0.2,  # boresight offset in V3 (arcsec)
             'fgs_num': 1,  # guider number
             'ra': 30.,  # RA of guide star (Allowed range: 0 - 360 degrees)
             'dec': 50.,  # Dec of guide star (Allowed range: 0 - 360 degrees)
             'pa': 2.,  # position angle of guide star (Allowed range: 0 - 360 degrees)
             'seg_num': 0}  # selected segment to guide on
        ``ct_uncert_fctr`` - the factor by which countrates are
            multiplied to determine the countrate uncertainty
        ``countrate_factor`` - the factor by which countrates are
            multipled by to simulate diffuse PSFs (e.g. in MIMF)
"""

# Standard Library Imports
import os
import logging

# Third Party Imports
import matplotlib
matplotlib.use("Qt5Agg")
import matplotlib.pyplot as plt
import pysiaf
from pysiaf.utils import rotations
import numpy as np
from astropy.io import ascii as asc
from astropy.table import Table
from PyQt5 import uic
from PyQt5.QtWidgets import QDialog

# Local Imports
from .. import coordinate_transforms, utils
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
    def __init__(self, program_id, observation_num, visit_num, segment_infile,
                 root=None, guide_star_params_dict=None, refonly=False,
                 selected_segs=None, vss_infile=None, out_dir=None,
                 ct_uncert_fctr=0.9, countrate_factor=None, oss_factor=0.6):
        """Initialize the class.

        Parameters
        ----------
        program_id : int
            APT program number
        observation_num : int
            Observation number
        visit_num : int
            Visit number
        segment_infile : str, optional
            Filepath to ALLpsfs.txt file with list of all segment locations
            and countrates
        root : str, optional
            Name used to generate output folder and output filenames. If not
            specified, will be derived from the segment_infile name
        guide_star_params_dict : dict, optional
            Dictionary containing guide star parameters, for example:
            {'v2_boff': 0.1,  # boresight offset in V2 (arcsec)
             'v3_boff': 0.2,  # boresight offset in V3 (arcsec)
             'fgs_num': 1,  # guider number
             'ra': 30.,  # RA of guide star (Allowed range: 0 - 360 degrees)
             'dec': 50.,  # Dec of guide star (Allowed range: 0 - 360 degrees)
             'pa': 2.,  # position angle of guide star (Allowed range: 0 - 360 degrees)
             'seg_num': 0}  # selected segment to guide on
        refonly : bool, optional
            Will the override file be written out using the 'ref-only' syntax?
        selected_segs : str, optional
            Filepath to regfile.txt file with list of locations and
            countrates for the selected segments (guide and reference stars)
        vss_infile : str, optional
            Filepath to guide star report provided by VSS
        out_dir : str, optional
            Location of out/ directory
        ct_uncert_fctr : float, optional
            The factor by which countrates are multiplied to determine
            the countrate uncertainty
        countrate_factor : float, optional
            The factor by which countrates are multiplied by to simulate
            diffuse PSFs (e.g. in MIMF)
        """

        # Initialize attributes
        self.root = root
        self.out_dir = utils.make_out_dir(out_dir, OUT_PATH, self.root)
        self.program_id = program_id
        self.observation_num = observation_num
        self.visit_num = visit_num
        self.refonly = refonly  # Implement "refonly" label for reference stars?
        self.ct_uncert_fctr = ct_uncert_fctr
        self.countrate_factor = countrate_factor
        self.oss_factor = oss_factor

        # Ensure the output directory exists
        utils.ensure_dir_exists(self.out_dir)

        # Parse the input file type (ALLpsfs.txt, regfile.txt, and VSS infile)
        self.get_gs_params(vss_infile, guide_star_params_dict)
        if segment_infile:
            self.parse_infile(segment_infile)
        if selected_segs:
            self.get_selected_segs(selected_segs)

        # Get aperture parameters from FGS SIAF
        if self.fgs_num:
            self.get_guider_aperture()
        else:
            LOGGER.warning('Segment Guiding: No guider information provided, creating photometry override file.')

    def get_chosen_segment_position(self):
        """Determine the V2/V3 position of the chosen segment.

        Check that the user-provided segment ID number is valid (0 to
        18). Calculate the central V2/V3 point of either the provided
        segment or the center of the segment array.
        """
        # Try to convert provided segment ID number to an integer
        try:
            seg_n = int(self.seg_num)
        except ValueError:
            raise ValueError('Unrecognized segment number: {}'.format(self.seg_num))

        # Refresh FGS data
        self.get_guider_aperture()

        # Ensure the provided segment ID is valid
        segment_max = len(self.v2_seg_array)
        if (seg_n < 0) or (seg_n > segment_max):
            msg = 'Segment number {} out of range (0, {})'.format(seg_n, segment_max)
            raise ValueError(msg)

        # Determine the central V2/V3 point from the given segment ID

        # If a specific segment was provided, set the V2/V3 ref point to be
        # that segment's location
        if seg_n > 0:
            self.v2_seg_n = self.v2_seg_array[seg_n - 1]
            self.v3_seg_n = self.v3_seg_array[seg_n - 1]
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
        self.v2_aim = self.v3_ref + dv3_aim

        # Convert to Ideal coordinates
        self.x_idl, self.y_idl = self.fgs_siaf_aperture.tel_to_idl(self.v2_aim, self.v2_aim)

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
        """Calculate the effective RAs and Decs for each segment, and
        write the segment guiding override file.
        """
        nseg = len(self.seg_id_array)

        # Convert V2/V3 coordinates to ideal coordinates
        idl_coords = self.fgs_siaf_aperture.tel_to_idl(self.v2_seg_array + self.v2_ref,
                                                       self.v3_seg_array + self.v3_ref)
        self.x_idl_segs, self.y_idl_segs = idl_coords

        # Get the attitude matrix
        attitude = rotations.attitude(self.v2_aim + self.v2_boff,
                                      self.v2_aim + self.v3_boff,
                                      self.ra, self.dec, float(self.pa))

        # Get RA and Dec for each segment.
        self.seg_ra = np.zeros(nseg)
        self.seg_dec = np.zeros(nseg)
        for i in range(nseg):
            V2 = self.v2_ref + self.v2_seg_array[i]
            V3 = self.v3_ref + self.v3_seg_array[i]
            self.seg_ra[i], self.seg_dec[i] = rotations.pointing(attitude, V2, V3,
                                                               positive_ra=True)

        # Convert segment coordinates to detector frame
        self.x_det, self.y_det = self.fgs_siaf_aperture.idl_to_det(self.x_idl_segs, self.y_idl_segs)

        # Check to make sure no segments are off the detector
        for x, y, i_seg in zip(self.x_det, self.y_idl_segs, self.seg_id_array):
            if x < 0.5 or x > 2048.5:
                LOGGER.warning('Segment Guiding: %8s off detector in X direction' % i_seg)
            if y < 0.5 or y > 2048.5:
                LOGGER.warning('Segment Guiding: %8s off detector in Y direction' % i_seg)

        # Check to make sure that RA is between 0 and 360 and Dec is between -90 and 90
        self.check_coords()
        self.nseg = nseg

    def write_override_file(self, nseg=None, verbose=True):
        """Write the segment guiding override file: {out_dir}/out/{root}/
        gs-override_{program_id}_{observation_num}_{visit_num}.txt

        Parameters
        ----------
        nseg : int
            The number of segments in the image
        verbose : bool, optional
            Log results of calculations and file content
        """
        # Define path and name of output override file
        out_file = 'gs-override-{}_{}_{}.txt'.format(self.program_id, self.observation_num,
                                                     self.visit_num)
        out_file = os.path.join(self.out_dir, out_file)

        # Print summary of input data (guide star RA, Dec, and PA, etc...)
        if verbose and self.fgs_num:
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

            if nseg:
                all_segments = 'All Segment Locations'
                all_segments += '\n                Segment     dV2    dV3    xIdl' +\
                                '   yIdl     RA         Dec         xDet     yDet'
                for p in range(nseg):
                    all_segments += ('\n                %5s    %6.2f %6.2f  %6.2f %6.2f  %10.6f %10.6f  %8.2f %8.2f'
                                     % (self.seg_id_array[p], self.v2_seg_array[p],
                                        self.v3_seg_array[p], self.x_idl_segs[p],
                                        self.y_idl_segs[p], self.seg_ra[p],
                                        self.seg_dec[p], self.x_idl_segs[p], self.y_idl_segs[p]))
                LOGGER.info('Segment Guiding: ' + all_segments)

        # Write out override file with RA/Decs of selected segments
        with open(out_file, 'w') as f:
            # Determine whether to include a multiplicative countrate factor
            countrate_qualifier = ' -count_rate_factor={:.3f}'.\
                format(self.countrate_factor) if self.countrate_factor else ''
            out_string = 'sts -gs_select {:4d}:{}:{}{}'.\
                format(int(self.program_id), self.observation_num,
                       self.visit_num, countrate_qualifier)
            if nseg:
                ##FIXME everything below needs to be optional
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
                # The OSS factor is multiplied by OSS to account for the 3x3 countrates
                # So we need to divide it out here and for every segment override file
                try:
                    rate = self.counts_array / self.oss_factor
                    uncertainty = self.counts_array * self.ct_uncert_fctr / self.oss_factor
                except AttributeError:
                    rate = [0.0] * nseg / self.oss_factor
                    uncertainty = [0.0] * nseg / self.oss_factor

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
                        label, seg, guide_seg_id + 1, self.seg_ra[guide_seg_id],
                        self.seg_dec[guide_seg_id], rate[guide_seg_id],
                        uncertainty[guide_seg_id])

                    if not self.refonly or (self.refonly and label == 'star'):
                        # Add list of segment IDs for all reference stars
                        for ref_seg_id in orientation:
                            if ref_seg_id != guide_seg_id:
                                star_string += ', %d' % (ref_seg_id + 1)

                    out_string += star_string

                f.write(out_string)

            if verbose:
                LOGGER.info('Segment Guiding: Guide Star Override: ' +
                            out_string.replace('-star', '\n                -star').
                            replace('-ref_only', '\n                -ref_only'))
                if nseg:
                    LOGGER.info('Segment Guiding: Saved {} segment commands to {}'.
                            format(len(orientations), out_file))

    def check_coords(self):
        """Check to make sure that RA is between 0 and 360 and Dec
        between -90 and 90
        """
        for i, ra in enumerate(self.seg_ra):
            if ra > 360.0:
                LOGGER.warning('Segment Guiding: RA = {}'.format(ra))
                self.seg_ra -= self.seg_ra
            elif ra < 0.0:
                LOGGER.warning('Segment Guiding: RA = {}'.format(ra))
                self.seg_ra += 360.0
            else:
                continue

        for i, dec in enumerate(self.seg_dec):
            if dec > 90.0:
                LOGGER.warning('Segment Guiding: Dec = {}'.format(dec))
                self.seg_dec -= 180.0
            elif dec < -90.0:
                LOGGER.warning('Segment Guiding: Dec = {}'.format(dec))
                self.seg_dec += 180.0
            else:
                continue

    def get_gs_params(self, vss_infile, guide_star_params_dict):
        """Get guide star parameters from dictionary or VSS file
        """
        if vss_infile and guide_star_params_dict:
            LOGGER.info(
                'Segment Guiding: Reading RA, Dec, and PA from VSS file {}'.
                format(vss_infile)
            )
            LOGGER.info(
                'Segment Guiding: Reading boresight offset and segment '
                'number from user-provided dictionary.'
            )
            self.get_guidestar_params_from_visit_file(vss_infile)
            self.v2_boff = guide_star_params_dict['v2_boff']
            self.v3_boff = guide_star_params_dict['v3_boff']
            self.seg_num = guide_star_params_dict['seg_num']

        elif vss_infile:
            LOGGER.info(
                'Segment Guiding: Reading RA, Dec, and PA from VSS file {}'.
                format(vss_infile)
            )
            LOGGER.info(
                'Segment Guiding: Setting boresight offset = 0 and segment '
                'number = 0.'
            )
            self.get_guidestar_params_from_visit_file(vss_infile)
            self.v2_boff = 0
            self.v3_boff = 0
            self.seg_num = 0

        elif guide_star_params_dict:
            LOGGER.info(
                'Segment Guiding: Reading all GS parameters from user-provided dictionary.'
            )
            # Map guide_star_params_dict keys to attributes
            for attr_name in guide_star_params_dict.keys():
                if attr_name == "seg_num" or attr_name == "fgs_num":
                    try:
                        setattr(self, attr_name, int(guide_star_params_dict[attr_name]))
                    except TypeError:
                        setattr(self, attr_name, guide_star_params_dict[attr_name])
                else:
                    try:
                        setattr(self, attr_name, float(guide_star_params_dict[attr_name]))
                    except TypeError:
                        setattr(self, attr_name, guide_star_params_dict[attr_name])

        else:
            LOGGER.info('If running the tool outside of the GUI, and '
                         'no dictionary of parameters in provided, will create'
                         'photometry override file.')

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
                self.counts_array = read_table['countrate']

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

        # Else if they are a regfile, parse it.
        elif os.path.exists(selected_segs):
            LOGGER.info('Segment Guiding: Guiding on segments selected in {}'.format(selected_segs))
            self.parse_regfile(selected_segs)

        # Otherwise, we don't know what it is...
        else:
            raise TypeError(
                'Unrecognized data type passed to selected_segs ({}); must be regfile.txt path or array of indices.'.
                format(selected_segs)
            )

    def parse_regfile(self, selected_segs):
        """Extract the segment positions and count rates from a regfile.txt

        Parameters
        ----------
        selected_segs : str
            Filepath to regfile.txt

        Raises
        ------
        TypeError
            Incompatible regfile.txt file provided as selected_segs
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
                                                   self.fgs_num)
            selected_segment_coords['V2Seg'], selected_segment_coords['V3Seg'] = v2, v3

        # If the positions aren't V2/V3 or x/y, what are they??
        else:
            raise TypeError('Incompatible regfile type: ', selected_segs)

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

    def get_guidestar_params_from_visit_file(self, visit_file):
        """Parse the Short Term Schedule file provided by VSS to
        determine the RA, Dec, and PA of the guide star.

        Parameters
        ----------
        visit_file : str
            Short term schedule file

        Raises
        ------
        ValueError
            Incompatible VSS file provided
        """
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
                raise ValueError(
                    'Provided visit file {} has unknown formatting; cannot parse.'.
                    format(visit_file)
                )

        # Read in as Astropy table
        names = ['Order', 'Star IDs', 'FGS', 'Corrected RA', 'Corrected Dec',
                 'Probability', 'ID V2', 'ID V3', 'ID X', 'ID Y', 'ID PA @ star',
                 'FGS Magnitude', 'FGS Mag Uncert', 'Count Rate', 'Count Rate Uncert']
        selected_guide_stars = asc.read(visit_file, data_start=i_start - 1, names=names)

        # Verify only one set of GS parameters are provided and extract them
        guide_star_params = []
        for col in ['Corrected RA', 'Corrected Dec', 'ID PA @ star', 'FGS']:
            values = set(selected_guide_stars[col])
            if len(values) > 1:
                raise ValueError(
                    'Cannot parse {} from input visit file; too many values provided: {}'.
                    format(col, list(values))
                )
            param = list(values)[0]
            guide_star_params.append(param)
        ra, dec, pa, fgs_num = guide_star_params

        # Update class attributes
        self.ra = ra
        self.dec = dec
        self.pa = pa
        self.fgs_num = fgs_num


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
        RAmean = self.seg_ra.mean()
        Decmean = self.seg_dec.mean()
        plt.plot(RAmean, Decmean, 'ro')
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

    def check_guidestar_params(self):
        """Ensure all guidestar parameters (RA, Dec, PA, and boresight
        offset) fall within appropriate ranges.
        """
        # ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ BORESIGHT OFFSET ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        # Guide star information
        msg = ["OK", "Boresight parameter conversion error",
               "Boresight parameter out of range"]

        if self.v2_boff:
            errcode = self.checkout(self.v2_boff, -10.0, 10.0)
            if errcode != 0:
                error = msg[errcode]
                raise ValueError(error)
            else:
                self.v2_boff = float(self.v2_boff)

        if self.v3_boff:
            errcode = self.checkout(self.v3_boff, -10.0, 10.0)
            if errcode != 0:
                error = msg[errcode]
                raise ValueError(error)
            else:
                self.v3_boff = float(self.v3_boff)

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

        # ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ COUNT RATE FACTOR ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        # Count rate factor has to be between 0 and 1 (required by VSS)
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


def convert_nrca3pixel_offset_to_v2v3_offset(x_offset, y_offset):
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


def open_segment_guiding_dialog(guider, program_id, observation_num, visit_num):
    """Raise a dialog window to prompt user for guide star parameters
    and other parameters needed to generate the override file.

    Parameters
    ----------
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
        program_id, observation_num, visit_num, ct_uncert_fctr,
        countrate_factor)
    """
    # Initialize dialog widget
    segment_guiding_dialog = QDialog()

    # Import .ui file
    uic.loadUi(os.path.join(__location__, 'segmentGuidingDialog.ui'), segment_guiding_dialog)

    # Set defaults from parsed header
    segment_guiding_dialog.lineEdit_programNumber.setText(str(program_id))
    segment_guiding_dialog.lineEdit_observationNumber.setText(str(observation_num))
    segment_guiding_dialog.lineEdit_visitNumber.setText(str(visit_num))

    # Run window and wait for response
    segment_guiding_dialog.exec()

    # Get parameters for dictionary from dialog
    try:
        # Parse what the boresight offset is
        if segment_guiding_dialog.radioButton_boresightNIRCam.isChecked():
            x_offset = float(segment_guiding_dialog.lineEdit_boresightX.text())
            y_offset = float(segment_guiding_dialog.lineEdit_boresightY.text())
            v2_offset, v3_offset = convert_nrca3pixel_offset_to_v2v3_offset(x_offset,
                                                                            y_offset)
            LOGGER.info(
                'Segment Guiding: Applying boresight offset of {}, {} arcsec (Converted from {}, {} pixels)'.
                format(v2_offset, v3_offset, x_offset, y_offset)
            )
        else:
            v2_offset = float(segment_guiding_dialog.lineEdit_boresightV2.text())
            v3_offset = float(segment_guiding_dialog.lineEdit_boresightV3.text())
            LOGGER.info(
                'Segment Guiding: Applying boresight offset of {}, {} arcsec'.
                format(v2_offset, v3_offset)
            )

        if segment_guiding_dialog.lineEdit_RA.text():
            ra = float(segment_guiding_dialog.lineEdit_RA.text())
        else:
            ra = None
        if segment_guiding_dialog.lineEdit_Dec.text():
            dec = float(segment_guiding_dialog.lineEdit_Dec.text())
        else:
            dec = None
        if segment_guiding_dialog.lineEdit_PA.text():
            pa = float(segment_guiding_dialog.lineEdit_PA.text())
        else:
            pa = None

        # Populate the parameter dictionary
        guide_star_params_dict = {
            'v2_boff': v2_offset,
            'v3_boff': v3_offset,
            'fgs_num': guider,
            'ra': ra,
            'dec': dec,
            'pa': pa,
            'seg_num': 0
        }

        # Get APT information and other necessary parameters
        program_id = segment_guiding_dialog.lineEdit_programNumber.text()
        observation_num = segment_guiding_dialog.lineEdit_observationNumber.text()
        visit_num = segment_guiding_dialog.lineEdit_visitNumber.text()
        ct_uncert_fctr = float(segment_guiding_dialog.lineEdit_countrateUncertainty.text())
        if segment_guiding_dialog.checkBox_countrateFactor.isChecked():
            countrate_factor = float(segment_guiding_dialog.doubleSpinBox_countrateFactor.value())
        else:
            countrate_factor = None
    except ValueError as e:
        if "could not convert string to float:" not in str(e):
            raise
        else:
            return

    return guide_star_params_dict, program_id, observation_num, visit_num, ct_uncert_fctr, countrate_factor

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# MAIN FUNCTION
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


def run_tool(segment_infile=None, guider=None, root=None, program_id=0, observation_num=0,
             visit_num=0, click_to_select_GUI=True, guide_star_params_dict=None,
             selected_segs=None, vss_infile=None, out_dir=None, data=None,
             masterGUIapp=None, refonly=False, ct_uncert_fctr=0.9,
             countrate_factor=None, parameter_dialog=True, oss_factor=0.6):
    """Run the segment guiding tool to generate a segment guiding
    override file.

    Parameters
    ----------
    segment_infile : str, optional
        Filepath to ALLpsfs.txt file with list of all segment locations
        and countrates
    guider : int, optional
        Which guider is being used: 1 or 2
    root : str, optional
        Name used to generate output folder and output filenames. If not
        specified, will be derived from the segment_infile name
    program_id : int, optional
        APT program number
    observation_num : int, optional
        Observation number
    visit_num : int, optional
        Visit number
    click_to_select_GUI : bool, optional
        Will the tool use the segment guiding GUI?
    guide_star_params_dict : dict, optional
        Reuired if parameter_dialog=False
        Dictionary containing guide star parameters, for example:
        {'v2_boff': 0.1,  # boresight offset in V2 (arcsec)
         'v3_boff': 0.2,  # boresight offset in V3 (arcsec)
         'fgs_num': 1,  # guider number
         'ra': 30.,  # RA of guide star (Allowed range: 0 - 360 degrees)
         'dec': 50.,  # Dec of guide star (Allowed range: 0 - 360 degrees)
         'pa': 2.,  # position angle of guide star (Allowed range: 0 - 360 degrees)
         'seg_num': 0}  # selected segment to guide on
    selected_segs : str, optional
        Filepath to regfile.txt file with list of locations and
        countrates for the selected segments (guide and reference stars)
    vss_infile : str, optional
        Filepath to guide star report provided by VSS
    out_dir : str, optional
        Location of out/ directory
    data : 2-D numpy array, optional
        Image that will be displayed in the click-to-select GUI
    masterGUIapp : qApplication, optional
        qApplication instance of parent GUI
    refonly : bool, optional
        Will the override file be written out using the 'ref-only' syntax?
    ct_uncert_fctr : float, optional
        The factor by which countrates are multiplied to determine
        the countrate uncertainty
    countrate_factor : float, optional
        The factor by which countrates are multiplied by to simulate
        diffuse PSFs (e.g. in MIMF)
    parameter_dialog : bool, optional
        Prompt the user to enter parameters (countrate factors, APT
        numbers, RA, Dec, PA, and boresight offset) from a dialog box
        rather than manually providing arguments
    oss_factor : float
        The factor that OSS applies to the 3x3 box counts in order to represent
        the full number of counts that we care about. This needs to be
        compensated for in the Sement Override file.
    """

    root = utils.make_root(root, segment_infile)
    utils.create_logger_from_yaml(__name__, root=root, level='DEBUG')

    if not guide_star_params_dict:
        guide_star_params_dict = {'v2_boff': None,
                                  'v3_boff': None,
                                  'fgs_num': None,
                                  'ra': None,
                                  'dec': None,
                                  'pa': None,
                                  'seg_num': None}
    try:
        # Run the SGT dialog to get the rest of the parameters from the user
        if parameter_dialog:
            guide_star_params_dict, program_id, observation_num, visit_num, \
                ct_uncert_fctr, countrate_factor = open_segment_guiding_dialog(
                    guider, program_id, observation_num, visit_num)

        if click_to_select_GUI:
            if data is None:
                raise ValueError(
                    'In order to run the segment guiding tool with '
                    '`click_to_select_GUI=True`, you must provide a 2-D '
                    'numpy array of the image to the `data` argument in '
                    '`segment_guiding.run_tool`.'
                )
            # Parse ALLpsfs.txt for locations of segments
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

            # Run the GUI to select guide and reference stars
            inds, seg_num = SegmentGuidingGUI.run_segment_override_gui(
                data, x, y, dist, selected_segs=selected_segs,
                masterGUIapp=masterGUIapp
            )
            LOGGER.info(
                'Segment Guiding: {} segment override commands generated with seg_num = {}'.
                format(len(inds), seg_num)
            )
            guide_star_params_dict['seg_num'] = seg_num

            # Turn index list into selected segments file
            selected_segs = np.array(inds)

        # Set up guiding calculator object
        sg = SegmentGuidingCalculator(program_id, observation_num, visit_num,
                                      segment_infile=segment_infile, root=root,
                                      guide_star_params_dict=guide_star_params_dict,
                                      selected_segs=selected_segs,
                                      vss_infile=vss_infile, out_dir=out_dir,
                                      refonly=refonly, ct_uncert_fctr=ct_uncert_fctr,
                                      countrate_factor=countrate_factor, oss_factor=oss_factor)
        # Verify all guidestar parameters are valid
        sg.check_guidestar_params()
        if guide_star_params_dict['ra'] and guide_star_params_dict['dec']: #FIXME
            sg.get_chosen_segment_position()
            sg.calculate_effective_ra_dec()
            sg.write_override_file(nseg=sg.nseg) # Print and save final output
            sg.plot_segments() # Save .pngs of plots
        else:
            sg.write_override_file()


    except Exception as e:
        LOGGER.exception(e)
        raise
