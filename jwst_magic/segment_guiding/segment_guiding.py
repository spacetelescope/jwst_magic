"""Generate segment guiding override files.

The JWST MAGIC Segment Guiding Tool is used during early wavefront
commissioning, when the mirror segments are still unstacked and/or unphased.
The tool takes a list of segment locations and guide star parameters,
using them to calculates the effective RA and Dec where each segment
appears in the sky. The user then specifies which segments to select as
the "guide star" and as the "reference stars". Finally, the tool generates
either a segment override file (SOF) or a photometry override file (POF),
in the form of gs_override*.txt file that the Visit Scheduling System
(VSS) will use to generate a visit.

Authors
-------
    - Colin Cox (original creator, May 2017)
    - Lauren Chambers (modifications in 2018)
    - Keira Brooks (modifications in 2018)
    - Shannon Osborne (modifications in 2020)

Use
---
    To generate a segment override file (SOF), the segment_guiding module
    can be used in the following way:
    ::
        from jwst_magic.segment_guiding import segment_guiding
        segment_guiding.generate_segment_override_file(
            segment_infile, guider, program_id, observation_num, visit_num,
            root=None, out_dir=None, selected_segs=None,
            guide_star_params_dict=None,
            threshold_factor=0.9, parameter_dialog=True,
            masterGUIapp=None):


    Or to generate a photometry override file (POF):
    ::
         segment_guiding.generate_photometry_override_file(
            root, program_id, observation_num, visit_num,
            countrate_factor, out_dir=None):
"""

# Standard Library Imports
from datetime import datetime
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
JENKINS = 'jenkins' in os.getcwd()
if matplotlib.get_backend() != 'Qt5Agg' and not JENKINS:
    matplotlib.use("Qt5Agg")
import matplotlib.path as mpltPath
import matplotlib.pyplot as plt
import numpy as np
import pysiaf
from pysiaf.utils import rotations

# Local Imports
if not JENKINS:
    from jwst_magic.segment_guiding import SegmentGuidingGUI
from jwst_magic.utils import coordinate_transforms, utils

# Establish segment guiding files directory
__location__ = os.path.realpath(os.path.join(os.getcwd(), os.path.dirname(__file__)))
PACKAGE_PATH = os.path.split(__location__)[0]
OUT_PATH = os.path.split(PACKAGE_PATH)[0]  # Location of out/ and logs/ directory

# Open the SIAF with pysiaf
FGS_SIAF = pysiaf.Siaf('FGS')


class SegmentGuidingCalculator:
    def __init__(self, override_type, program_id, observation_num, visit_num,
                 root, out_dir, segment_infile_list=None, guide_star_params_dict=None,
                 selected_segs_list=None, threshold_factor=0.9, countrate_factor=None,
                 countrate_uncertainty_factor=None, log=None):
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
        segment_infile_list : list of str, optional
            Filepath(s) to all_found_psfs*.txt file with list of all
            segment locations and countrates
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
        selected_segs_list : list of str, optional
            List of filepath(s) to guiding_selections*.txt files with list of
            locations and countrates for the selected segments (guide and
            reference stars). Used for SOF Generation
        threshold_factor : float, optional
            The factor by which countrates are multiplied to determine
            the countrate uncertainty
            Used for SOF Generation
        countrate_factor : float, optional
            The factor by which countrates are multiplied by to simulate
            diffuse PSFs (e.g. in MIMF)
            Used for POF Generation
        countrate_uncertainty_factor : float, optional
            The factor by which countrate uncertainties are multiplied by to simulate
            diffuse PSFs (e.g. in MIMF)
            Used for POF Generation
        log : logger object
            Pass a logger object (output of tils.create_logger_from_yaml) or a new log
            will be created
        """

        # Initialize parameters into attributes
        self.override_type = override_type
        self.program_id = int(program_id)
        self.observation_num = observation_num if observation_num != '' else None
        self.visit_num = visit_num  if visit_num != '' else None
        self.root = root
        self.out_dir = out_dir
        self.threshold_factor = threshold_factor
        self.countrate_factor = countrate_factor
        self.countrate_uncertainty_factor = countrate_uncertainty_factor

        # Start logger
        if log is None:
            self.log = logging.getLogger(__name__)
        else:
            self.log = log

        # Initialize other attributes
        # Will the override file be written out using the 'ref-only' syntax? Should be yes.
        self._refonly = True

        # Ensure the output directory exists
        utils.ensure_dir_exists(self.out_dir)

        # Set up to do segment override calculations
        if self.override_type == "SOF":
            self.get_gs_params(guide_star_params_dict)
            self.get_guider_aperture()
            self.parse_infile(segment_infile_list)
            self.get_selected_segs(selected_segs_list)

    def get_center_pointing(self):
        """Determine the V2/V3 position of the chosen segment.

        Check that the user-provided segment ID number is valid (0 to
        18). Calculate the central V2/V3 point of either the provided
        segment or the center of the segment array.
        """

        # Ensure the provided segment ID is valid
        flat_v2_array = [val for l in self.v2_seg_array for val in l]
        segment_max = int(len(flat_v2_array) / self._num_infiles) # same number of founds psfs for each config
        if isinstance(self.seg_num, int):
            if (self.seg_num < 0) or (self.seg_num > segment_max):
                msg = 'Segment number {} out of range (0, {})'.format(self.seg_num, segment_max)
                raise ValueError(msg)
        elif not isinstance(self.seg_num, list):
            raise ValueError('Center of pointing {} but either be of type int or list, not {}'.format(self.seg_num,
                                                                                                    type(self.seg_num)))

        # Determine the central V2/V3 point from the given segment ID

        # FIXME: is this wrong? adding the V2/V3 reference twice...?

        # If a specific segment was provided, set the V2/V3 ref point to be
        # that segment's location
        self.v2_seg_n, self.v3_seg_n = [], []
        self.v2_aim, self.v3_aim = [], []
        self.x_idl_aim, self.y_idl_aim = [], []

        for i, (v2_seg_array, v3_seg_array) in enumerate(zip(self.v2_seg_array, self.v3_seg_array)):
            if isinstance(self.seg_num, list):
                v2_seg_n, v3_seg_n = coordinate_transforms.Raw2Tel(self.seg_num[0], self.seg_num[1], self.fgs_num)
            elif self.seg_num > 0:
                v2_seg_n = v2_seg_array[self.seg_num - 1]
                v3_seg_n = v3_seg_array[self.seg_num - 1]
            # Otherwise, if the input segment ID was 0, set the V2/V3 ref point to
            # be the mean of all segments' locations
            else:
                v2_seg_n = v2_seg_array.mean()
                v3_seg_n = v3_seg_array.mean()

            # Calculate aim
            v2_aim = self.v2_ref + v2_seg_n # v2_seg_n = dv2_aim
            v3_aim = self.v3_ref + v3_seg_n # v3_seg_n = dv3_aim

            # Convert to Ideal coordinates
            x_idl_aim, y_idl_aim = self.fgs_siaf_aperture.tel_to_idl(v2_aim, v3_aim)

            # Save variables to list attributes
            self.v2_seg_n.append(v2_seg_n)
            self.v3_seg_n.append(v3_seg_n)
            self.v2_aim.append(v2_aim)
            self.v3_aim.append(v3_aim)
            self.x_idl_aim.append(x_idl_aim)
            self.y_idl_aim.append(y_idl_aim)

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
        self.v2_ref = float(self.fgs_siaf_aperture.V2Ref)
        self.v3_ref = float(self.fgs_siaf_aperture.V3Ref)
        self.v3_idl_yangle = self.fgs_siaf_aperture.V3IdlYAngle
        self.v_idl_parity = self.fgs_siaf_aperture.VIdlParity

    def calculate_effective_ra_dec(self):
        """Calculate the effective RAs and Decs for each segment.
        Loop through different configurations and save attributes in lists
        """
        self.n_segments = []
        self.seg_ra, self.seg_dec = [], []
        self.x_idl_segs, self.y_idl_segs = [], []
        self.x_det_segs, self.y_det_segs = [], []

        for i, seg_id_array in enumerate(self.seg_id_array):
            n_segments = len(seg_id_array)

            # Get the attitude matrix
            attitude = rotations.attitude(self.v2_aim[i] + self.v2_boff,
                                          self.v3_aim[i] + self.v3_boff,
                                          self.ra, self.dec, self.pa)

            # Get RA and Dec for each segment.
            seg_ra = np.zeros(n_segments)
            seg_dec = np.zeros(n_segments)
            for j in range(n_segments):
                v2 = self.v2_ref + self.v2_seg_array[i][j]
                v3 = self.v3_ref + self.v3_seg_array[i][j]
                seg_ra[j], seg_dec[j] = rotations.pointing(attitude, v2, v3,
                                                                     positive_ra=True)

            # Convert V2/V3 coordinates to ideal coordinates and detector frame coordinates
            idl_coords = self.fgs_siaf_aperture.tel_to_idl(self.v2_seg_array[i] + self.v2_ref,
                                                           self.v3_seg_array[i] + self.v3_ref)
            x_idl_segs, y_idl_segs = idl_coords
            x_det_segs, y_det_segs = self.fgs_siaf_aperture.idl_to_det(x_idl_segs,y_idl_segs)

            # Check to make sure all the computed segment locations are within
            # the needed FOV
            self.check_segments_inside_fov(attitude, x_det_segs, y_det_segs, seg_id_array, seg_ra, seg_dec)

            self.n_segments.append(n_segments)
            self.seg_ra.append(seg_ra)
            self.seg_dec.append(seg_dec)
            self.x_idl_segs.append(x_idl_segs)
            self.y_idl_segs.append(y_idl_segs)
            self.x_det_segs.append(x_det_segs)
            self.y_det_segs.append(y_det_segs)

    def write_override_file(self, verbose=True):
        """Write the segment guiding override file: {out_dir}/out/{root}/
        gs_override_{program_id}_{observation_num}_{visit_num}.txt

        Parameters
        ----------
        verbose : bool, optional
            Log results of calculations and file content
        """
        # Split multiple specified observations up
        obs_num_list, obs_list_name = self._split_obs_num(self.observation_num)

        # Define path and name of output override file
        # Only add observations and visits to the list if an observation is specified
        out_file = '{}_gs_override_{}'.format(datetime.now().strftime('%Y%m%d'),
                                              self.program_id)
        if obs_list_name is not None:
            out_file += "_{}".format(obs_list_name)

            if len(obs_num_list) > 1:
                out_file += "_1"
                if self.visit_num is not None:
                    if int(self.visit_num) != 1:
                        self.log.warning('Visit number set to 001. You cannot specify '
                                       'visit numbers if specifying multiple observations.')
                self.visit_num = 1

            elif self.visit_num is not None:
                out_file += "_{}".format(self.visit_num)

        out_file += ('.txt')

        out_file = os.path.join(self.out_dir, out_file)

        # Update variables
        if self.override_type == "SOF":

            # Flatten data into 1 list of all PSFs
            self.countrate_array_flat = [x for n in self.countrate_array for x in n.tolist()]
            self.n_segments_flat = sum(self.n_segments)
            self.seg_id_array_flat = np.arange(len(
                [x for n in self.seg_id_array for x in n.tolist()])) + 1  # Re-do numbering so all have unique numbers
            self.v2_seg_array_flat = [x for n in self.v2_seg_array for x in n.tolist()]
            self.v3_seg_array_flat = [x for n in self.v3_seg_array for x in n.tolist()]
            self.x_idl_segs_flat = [x for n in self.x_idl_segs for x in n.tolist()]
            self.y_idl_segs_flat = [x for n in self.y_idl_segs for x in n.tolist()]
            self.seg_ra_flat = np.array(self.seg_ra).flatten()
            self.seg_dec_flat = np.array(self.seg_dec).flatten()
            self.x_det_segs_flat = [x for n in self.x_det_segs for x in n.tolist()]
            self.y_det_segs_flat = [x for n in self.y_det_segs for x in n.tolist()]

            for radec in zip(self.seg_ra_flat, self.seg_dec_flat):
                if list(zip(self.seg_ra_flat, self.seg_dec_flat)).count(radec) > 1:
                    inds = [i for i, x in enumerate(list(zip(self.seg_ra_flat, self.seg_dec_flat))) if x == radec]
                    for i in inds[1:]:
                        self.seg_id_array_flat[i] = int(self.seg_id_array_flat[inds[0]])

            # Re-match numbering in self.selected_segment_ids
            self.selected_segment_ids_flat = []
            for i, config in enumerate(self.selected_segment_ids):
                new_config = []
                for psf_ind in config:
                    shift = sum(self.n_segments[:i])
                    new_config.append(self.seg_id_array_flat[shift + psf_ind] - 1)
                self.selected_segment_ids_flat.append(new_config)

            #  Print summary of input data (guide star RA, Dec, and PA, etc...)
            if verbose:
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
                self.log.info('Segment Guiding: ' + summary_output)

                all_segments = 'All Segment Locations'
                all_segments += '\n                Segment     dV2    dV3    xIdl' +\
                                '   yIdl     RA         Dec         xDet     yDet'
                for p in range(self.n_segments_flat):
                    all_segments += ('\n                %5s    %6.2f %6.2f  %6.2f %6.2f  %10.6f %10.6f  %8.2f %8.2f'
                                     % (self.seg_id_array_flat[p], self.v2_seg_array_flat[p],
                                        self.v3_seg_array_flat[p], self.x_idl_segs_flat[p],
                                        self.y_idl_segs_flat[p], self.seg_ra_flat[p],
                                        self.seg_dec_flat[p], self.x_det_segs_flat[p], self.y_det_segs_flat[p]))
                self.log.info('Segment Guiding: ' + all_segments)

        # Determine what the commands should be:
        if obs_num_list is None:
            obs_list = [None]
            visit_list = [None]
        elif len(obs_num_list) == 1:
            obs_list = [self.observation_num]
            visit_list = [self.visit_num]
        else:
            obs_list = obs_num_list
            visit_list = [1] * len(obs_num_list)

        # Write out override file with RA/Decs of selected segments
        with open(out_file, 'w') as f:
            # Determine whether to include a multiplicative countrate factor
            countrate_qualifier = ' -count_rate_factor={:.4f}'.\
                format(self.countrate_factor) if self.countrate_factor else ''
            countrate_uncertainty_qualifier = ' -count_rate_uncertainty_factor={:.4f}'.\
                format(self.countrate_uncertainty_factor) if self.countrate_uncertainty_factor else ''
            out_string = 'sts -gs_select '

            for o, v in zip(obs_list, visit_list):
                if out_string[-2:] == 't ':
                    out_string += "{:05d}".format(self.program_id)
                else:
                    out_string += ", {:05d}".format(self.program_id)

                if o is not None:
                    out_string +="{:03d}".format(int(o))
                    if v is not None:
                        out_string +="{:03d}".format(int(v))

            out_string += (countrate_qualifier)
            out_string += (countrate_uncertainty_qualifier)

            if self.override_type == "SOF":
                # Determine which segments have been selected
                orientations = list(self.selected_segment_ids_flat)
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
                # Note: This used to be where the oss factor was divided out, but that is now removed
                rate = self.countrate_array_flat
                uncertainty = np.array(self.countrate_array_flat) * self.threshold_factor

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
                        label, seg, guide_seg_id + 1, self.seg_ra_flat[guide_seg_id],
                        self.seg_dec_flat[guide_seg_id], rate[guide_seg_id],
                        uncertainty[guide_seg_id])

                    if not self._refonly or (self._refonly and label == 'star'):
                        # Add list of segment IDs for all reference stars
                        for ref_seg_id in orientation:
                            if ref_seg_id != guide_seg_id:
                                star_string += ', %d' % (ref_seg_id + 1)

                    out_string += star_string

                # Write out the override report
                self.write_override_report(out_file, orientations, n_guide_segments, obs_list_name)

            f.write(out_string)

            if verbose:
                self.log.info('Segment Guiding: Guide Star Override: ' +
                            out_string.replace('-star', '\n                -star').
                            replace('-ref_only', '\n                -ref_only'))
                self.log.info('Segment Guiding: Saved override command to {}'.
                            format(out_file))


    def write_override_report(self, filename, orientations, n_guide_segments, obs_list_name):
        # Define path and name of output override report
        file_root = filename.split('.txt')[0]
        out_file = '{}_REPORT.txt'.format(file_root)
        out_file = os.path.join(self.out_dir, out_file)

        username = getpass.getuser()

        with open(out_file, 'w') as f:
            f.write('Guide Star Override Report\n')
            f.write('Generated on {} at {} by {}\n'.format(time.strftime("%Y/%m/%d"), time.strftime("%H:%M:%S"), username))
            f.write('\n')
            f.write('{:14s}: {:s}'.format('Program ID', str(self.program_id)) + '\n')
            f.write('{:14s}: {:s}'.format('Observation #', str(obs_list_name)) + '\n')
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

                values = [label + str(seg), int(guide_seg_id + 1), self.seg_ra_flat[guide_seg_id], self.seg_dec_flat[guide_seg_id],
                           self.x_idl_segs_flat[guide_seg_id], self.y_idl_segs_flat[guide_seg_id],
                           -self.x_idl_segs_flat[guide_seg_id], self.y_idl_segs_flat[guide_seg_id],
                           self.x_det_segs_flat[guide_seg_id], self.y_det_segs_flat[guide_seg_id]]

                row_string_to_format = '{:<13s}| ' + '{:<13d}| ' +  '{:<13f}| '* (len(values) - 2)
                row_string = row_string_to_format[:-2].format(*values) + '\n'
                f.write(row_string)

    def check_segments_inside_fov(self, attitude, x_det_segs, y_det_segs, seg_id_array, seg_ra, seg_dec):
        """Check to make sure that the calculated RA and Dec of each
        segment is within the field of view of the given FGS.

        Parameters
        ----------
        attitude : 3 x 3 numpy array
            Attitude matrix generated by pysiaf.utils.rotations.attitude
        """
        # Check to make sure no segments are off the detector, pixel-wise
        for x, y, i_seg in zip(x_det_segs, y_det_segs, seg_id_array):
            if x < 0.5 or x > 2048.5:
                self.log.warning('Segment Guiding: %8s off detector in X direction' % i_seg)
            if y < 0.5 or y > 2048.5:
                self.log.warning('Segment Guiding: %8s off detector in Y direction' % i_seg)

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

        # Determine if every segment RA and Dec is within that Path (i.e. within the guider FOV)
        seg_pointings = [(seg_ra, seg_dec) for seg_ra, seg_dec in zip(seg_ra, seg_dec)]
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
        self.seg_num = guide_star_params_dict['seg_num']

    def parse_infile(self, segment_infile_list):
        """Get the segment positions and count rates from a file.

        Parameters
        ----------
        segment_infile_list : list of str
            List of files containing all segment locations and count rates

        Raises
        ------
        TypeError
            Incompatible file type provided as segment_infile
        """
        self.countrate_array = []
        self.seg_id_array = []
        self.v2_seg_array = []
        self.v3_seg_array = []
        self._num_infiles = 0

        for segment_infile in segment_infile_list:
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
                    self.countrate_array.append(read_table['countrate'])

            else:
                raise TypeError('Incompatible file type: ', segment_infile)

            # Define the IDs and coordinates of all segments
            self.seg_id_array.append(segment_coords['SegID'])
            self.v2_seg_array.append(segment_coords['V2Seg'])
            self.v3_seg_array.append(segment_coords['V3Seg'])
            self._num_infiles += 1

            self.log.info(
                'Segment Guiding: {} segment coordinates read from {}'.
                format(len(segment_coords['SegID']), segment_infile)
            )

    def get_selected_segs(self, selected_segs):
        """If a file of selected segments has been provided, get their
        locations and count rates. Output variable self.selected_segment_ids
        will count ids from 0->17

        Parameters
        ----------
        selected_segs_list : list of str and/or an array
            list of file(s) containing locations and count rates of
            selected segments. If an array, reminder that the ids
            should go from 1 -> 18

        Raises
        ------
        TypeError
            Incompatible file type provided as selected_segs
        """

        # If the selected segments are an array of lists
        if isinstance(selected_segs, np.ndarray):
            self.log.info('Segment Guiding: Guiding on an array of segments')
            # If there is more than one orientation provided
            if isinstance(selected_segs[0], list):
                self.selected_segment_ids = []
                for i, orientation in enumerate(selected_segs):
                    self.selected_segment_ids.append([s - 1 for s in orientation])

            # If there is only one
            else:
                self.selected_segment_ids = selected_segs - 1

        # Else if they are a list of guiding_selections*.txt files, parse them.
        elif bool(selected_segs) and all(isinstance(elem, str) for elem in selected_segs):
            self.selected_segment_ids =[]
            for i, path in enumerate(selected_segs):
                if os.path.exists(path):
                    selected_segs_ids = self.parse_guiding_selections_file(path, i)
                    self.selected_segment_ids.append(selected_segs_ids)
                else:
                    # remove that path from the selected segs
                    selected_segs.pop(i)

        # Otherwise, we don't know what it is...
        else:
            raise TypeError(
                'Unrecognized data type passed to selected_segs ({}); must be guiding_selections*.txt path or array of indices.'.
                format(selected_segs)
            )

    def parse_guiding_selections_file(self, selected_segs, num_file):
        """Extract the segment positions and count rates from a guiding_selections*.txt

        Parameters
        ----------
        selected_segs : str
            Filepath to guiding_selections*.txt
        num_file : int
            The index of the guiding_selections*.txt file in the input list.
            Used to match PSFs in selections file to v2_seg_array and v3_seg_array
            which were read in from the all_found_psfs file

        Raises
        ------
        TypeError
            Incompatible guiding_selections*.txt file provided as selected_segs
        """
        # Make sure the file is ascii-readable
        n_segs = len(self.seg_id_array[num_file])
        try:
            read_selected_segs = asc.read(selected_segs)
            column_names = read_selected_segs.colnames
            self.log.info(
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
                if coords[0]['V2Seg'] == self.v2_seg_array[num_file][i_seg] and \
                   coords[0]['V3Seg'] == self.v3_seg_array[num_file][i_seg]:
                    selected_segs_ids.append(i_seg)
        if len(selected_segs_ids) == 0:
            raise TypeError(
                'Coordinates of selected segments file do not match those of '
                'the provided input file.'
            )

        return selected_segs_ids

    def plot_segments(self):
        """Generate and save plots of segments in V2/V3 and RA/Dec.
        """
        # Plot segments in V2/V3 frame
        for j, (v2_seg_array, v3_seg_array, seg_id_array, seg_ra, seg_dec) in \
                enumerate(zip(self.v2_seg_array, self.v3_seg_array, self.seg_id_array, self.seg_ra, self.seg_dec)):
            plt.figure(1)
            plt.clf()
            plt.plot(v2_seg_array, v3_seg_array, 'b*')
            plt.plot(v2_seg_array.mean(), v3_seg_array.mean(), 'ro')
            plt.grid(True)
            plt.axis([40.0, -40.0, -40.0, 40.0])  # V2 to the left
            plt.title('Segments')
            plt.xlabel('<-- Delta V2')
            plt.ylabel('Delta V3 -->')
            for s in range(len(v2_seg_array)):
                plt.text(v2_seg_array[s], v3_seg_array[s], seg_id_array[s])
            plt.savefig(os.path.join(self.out_dir, self.root + '_V2V3segments_config{}.png'.format(j+1)))

            # Plot calculate segments' RA and Dec
            plt.figure(2)
            plt.clf()
            plt.plot(seg_ra, seg_dec, 'b*')
            ra_mean = seg_ra.mean()
            dec_mean = seg_dec.mean()
            plt.plot(ra_mean, dec_mean, 'ro')
            for i in range(len(v2_seg_array)):
                plt.text(seg_ra[i], seg_dec[i], str(i + 1))
            plt.grid(True)
            plt.title('Segment RA and Dec')
            plt.xlabel('RA')
            plt.ylabel('Dec')
            plt.gca().invert_xaxis()
            plt.gca().ticklabel_format(useOffset=False)
            plt.savefig(os.path.join(self.out_dir, self.root + '_RADecsegments_config{}.png'.format(j+1)))

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

            for val in ['observation_num', 'visit_num']:
                attr = getattr(self, val)
                if attr is not None:
                    try:
                        int(attr)
                    except ValueError:
                        raise ValueError('Invalid input for SOF: multiple {} numbers not allowed'.format(val.split('_')[0]))
                else:
                    raise ValueError('Invalid input for SOF: {} number must be set. It cannot be None.'.format(val.split('_')[0]))

        elif override_type == "POF":
            # ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ COUNT RATE FACTOR ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
            # Count rate factor has to be between 0 and 1 (required by VSS)
            msg = ["OK", "Countrate factor conversion error",
                   "Countrate factor out of range"]
            if self.countrate_factor is not None:
                errcode = self.checkout(self.countrate_factor, 0.0, 1.0)
                if errcode != 0:
                    error = "{} for count_rate_factor. Expecting between 0.0 and 1.0.".format(msg[errcode])
                    raise ValueError(error)
            if self.countrate_uncertainty_factor is not None:
                errcode = self.checkout(self.countrate_uncertainty_factor, 0.01, 1.0, low_inclusive=True)
                if errcode != 0:
                    error = "{} for count_rate_uncertainty_factor. Expecting between 0.01 (inclusive) and 1.0.".format(msg[errcode])
                    raise ValueError(error)

    @staticmethod
    def checkout(value, low, high, low_inclusive=False):
        """Test conversion from string to float. If float conversion
        works and is within test range, return errcode 0 for OK, 1 for conversion
        error, 2 for range error. The 'high' bounds are exclusive.

        Parameters
        ----------
        value: float
            The value to be checked
        low: float
            Lower bounds
        high: float
            Upper bounds
        low_inclusive: boolean, False
            Whether or not to include the lower boundary condition
        """

        try:
            x = float(value)
            if low_inclusive:
                if low <= x < high:
                    errcode = 0
                else:
                    errcode = 2
            else:
                if low < x < high:
                    errcode = 0
                else:
                    errcode = 2
        except ValueError:
            errcode = 1
        return errcode

    @staticmethod
    def _split_obs_num(obs_num):
        """Parse strings listing one or more observation numbers, with single numbers
        separated by commas and ranges separated by hyphens.

        Parameters
        ----------
        obs_num : str
            String of one or more observation numbers and/or ranges

        Returns
        -------
        final_num_list : list
            List of all observation numbers denoted in obs_num
        obs_list_string : str
            String containing ordered list of all observations and ranges
        """
        # If there is not observation number specific, don't bother
        if obs_num is None:
            return None, None

        # First, divide things up by commas
        comma_split_obs_num_list = str(obs_num).split(',')

        # Then, deal with ranges
        final_num_list = []
        for obs_str in comma_split_obs_num_list:
            if '-' in obs_str:
                bounds = obs_str.split('-')
                obs_range = list(np.arange(int(bounds[0]), int(bounds[1]) + 1))
                final_num_list.extend(obs_range)
            else:
                # And be sure to make everything an int!
                final_num_list.append(int(obs_str))

        final_num_list = sorted(np.unique(final_num_list))

        # Finally, create a properly-formatted string including all observations
        obs_list_string = ' '
        for i, obs_int in enumerate(final_num_list):
            incremental = obs_int == final_num_list[i - 1] + 1
            last_number = obs_int == final_num_list[-1]
            in_range = obs_list_string[-1] == '-'

            if i == 0:
                obs_list_string = str(obs_int)
            elif incremental and not in_range and last_number:
                obs_list_string += '-{}'.format(obs_int)
            elif incremental and not in_range:
                obs_list_string += '-'
            elif incremental and not last_number:
                continue
            elif incremental and last_number:
                obs_list_string += str(obs_int)
            elif in_range:
                obs_list_string += str(final_num_list[i - 1])
                obs_list_string += ',{}'.format(obs_int)
            else:
                obs_list_string += ',{}'.format(obs_int)

        return final_num_list, obs_list_string


def generate_segment_override_file(segment_infile_list, guider,
                                   program_id, observation_num, visit_num,
                                   ra=None, dec=None,
                                   root=None, out_dir=None, selected_segs_list=None,
                                   guide_star_params_dict=None, threshold_factor=0.9,
                                   parameter_dialog=True, dialog_obj=None,
                                   master_gui_app=None, log=None):
    """Run the segment guiding tool to select guide and reference stars and
    generate a segment guiding override file.

    Parameters
    ----------
    segment_infile_list : list of str
        File path(s) to all_found_psfs*.txt file with list of all segment locations
        and countrates
    guider : int
        Which guider is being used: 1 or 2
    program_id : int or str
        APT program number
    observation_num : int or str
        Observation number
    visit_num : int or str
        Visit number
    ra :  float, optional
        RA of guide star
    dec :  float, optional
        DEC of guide star
    root : str, optional
        Name used to generate output folder and output filenames. If not
        specified, will be derived from the segment_infile name.
    out_dir : str, optional
        Location of out/ directory. If not specified, will be placed
        within the repository: .../jwst_magic/out/
    selected_segs_list : list of str or array
        List of file path(s) to guiding_selections*.txt files with list of
        locations and countrates for the selected segments (guide and
        reference stars). Can also be a array of lists, 1 list per command.
        Reminder that the length of selected_segs_list should match the length
        of segment_infile_list, both are the length of the number of commands.
        If
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
    dialog_obj : SegmentGuidingDialog object, optional
        If parameter_dialog is True, can pass a pre-set dialog object or
        re-create it if dialog_obj=None
    master_gui_app : qApplication, optional
        qApplication instance of parent GUI
    log : logger object
        Pass a logger object (output of utils.create_logger_from_yaml) or a new log
        will be created
    """

    root = utils.make_root(root, segment_infile_list[0])  # pull the root from the 0th all psfs file
    out_dir = utils.make_out_dir(out_dir, OUT_PATH, root)

    if log is None:
        log = logging.getLogger(__name__)

    try:
        # Get the guide star parameters
        if parameter_dialog:
            if dialog_obj is None:
                dialog_obj = SegmentGuidingGUI.SegmentGuidingDialog(
                    "SOF", guider, program_id, observation_num, visit_num, ra, dec, log=log
                )
            accepted = dialog_obj.exec()
            params = dialog_obj.get_dialog_parameters() if accepted else None

            if params is not None:
                guide_star_params_dict, program_id, observation_num, visit_num, \
                    threshold_factor, _, _ = params

                # Overwrite default seg_num information
                seg_num_path = os.path.join(out_dir, 'center_pointing_{}_G{}.txt'.format(root, guider))
                if os.path.exists(seg_num_path):
                    log.info(
                        'Segment Guiding: Pulling center of pointing information from {}'.format(seg_num_path))
                    in_table = asc.read(seg_num_path, format='commented_header', delimiter=',')
                    try:
                        guide_star_params_dict['seg_num'] = int(in_table['segnum'][0])
                    except ValueError:
                        guide_star_params_dict['seg_num'] = [float(i) for i in in_table['segnum'][0].split(' ')]
                else:
                    utils.write_cols_to_file(seg_num_path, labels=['segnum'], cols=[0], log=log)
                    log.warning(
                        "Segment Guiding: Couldn't find center of pointing file center_pointing_{}_G{}.txt. "
                        "Assuming seg_num = 0 and writing out the file.".format(root, guider))

            else:
                log.warning('Segment Guiding: SOF creation cancelled.')
                return

        elif guide_star_params_dict is None:
            raise ValueError(
                'In order to run the segment guiding tool with '
                '`parameter_dialog=False`, you must provide a dictionary '
                'of guide star parameters to the `guide_star_params_dict` '
                'argument in `segment_guiding.generate_segment_override_file()`.'
            )

        # Check if there is an existing file with the same prog/obs/visit
        if not JENKINS:
            overwrite_existing_file = SegmentGuidingGUI.check_override_overwrite(
                out_dir, program_id, observation_num, visit_num, logger=log
            )
            if overwrite_existing_file:
                return

        # Determine which segments are the guide and reference segments
        if selected_segs_list is None:
            raise ValueError(
                'In order to run the segment guiding tool, you must '
                'provide a list of files specifying the locations and count rates '
                'of the guide and reference stars as the `selected_segs_list` '
                'argument in `segment_guiding.generate_segment_override_file()`.'
            )

        # Set up guiding calculator object
        sg = SegmentGuidingCalculator(
            "SOF", program_id, observation_num, visit_num, root, out_dir,
            segment_infile_list=segment_infile_list,
            guide_star_params_dict=guide_star_params_dict,
            selected_segs_list=selected_segs_list, threshold_factor=threshold_factor,
            log=log
        )
        # Verify all guidestar parameters are valid
        sg.check_guidestar_params("SOF")

        # Determine the V2/V3 of the pointing center
        sg.get_center_pointing()

        # Write a SOF
        sg.calculate_effective_ra_dec()
        sg.write_override_file()  # Print and save final output

        if not JENKINS:
            sg.plot_segments()  # Save .pngs of plots

    except Exception as e:
        log.exception(e)
        raise


def generate_photometry_override_file(root, program_id, observation_num, visit_num,
                                      countrate_factor=None,
                                      countrate_uncertainty_factor=None,
                                      out_dir=None, parameter_dialog=True,
                                      dialog_obj=None, log=None):
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
        within the repository: .../jwst_magic/out/
    parameter_dialog : bool, optional
        Prompt the user to enter parameters (countrate factors, APT
        numbers) from a dialog box rather than manually providing arguments.
        Required if countrate_factor=None
    dialog_obj : SegmentGuidingDialog object, optional
        If parameter_dialog is True, can pass a pre-set dialog object or
        re-create it if dialog_obj=None
    log : logger object
        Pass a logger object (output of tils.create_logger_from_yaml) or a new log
        will be created
    """

    if log is None:
        log = logging.getLogger(__name__)
    out_dir = utils.make_out_dir(out_dir, OUT_PATH, root)

    try:
        # Get the program parameters and countrate factor
        if parameter_dialog:
            if dialog_obj is None:
                dialog_obj = SegmentGuidingGUI.SegmentGuidingDialog(
                    "POF", None, program_id, observation_num, visit_num, log=log
                )
            accepted = dialog_obj.exec()
            params = dialog_obj.get_dialog_parameters() if accepted else None

            _, program_id, observation_num, visit_num, _, countrate_factor, countrate_uncertainty_factor = params

        elif countrate_factor is None:
            raise ValueError(
                'In order to run the segment guiding tool with '
                '`parameter_dialog=False`, you must provide a countrate factor '
                'between 0.0 and 1.0 to the `countrate_factor` argument in '
                '`segment_guiding.generate_photometry_override_file()`.'
            )
        elif countrate_uncertainty_factor is None:
            raise ValueError(
                'In order to run the segment guiding tool with '
                '`parameter_dialog=False`, you must provide a countrate factor '
                'between 0.01 (inclusive) and 1.0 to the `countrate_uncertainty_factor` argument in '
                '`segment_guiding.generate_photometry_override_file()`.'
            )

        # Set up guiding calculator object
        sg = SegmentGuidingCalculator(
            "POF", program_id, observation_num, visit_num, root, out_dir,
            countrate_factor=countrate_factor, countrate_uncertainty_factor=countrate_uncertainty_factor, log=log
        )
        # Verify all guidestar parameters are valid
        sg.check_guidestar_params("POF")

        # Write the POF
        sg.write_override_file()

    except Exception as e:
        log.exception(e)
        raise
