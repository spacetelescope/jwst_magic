"""Interactive segment guiding dialog for creating SOFs and POFs

Builds a GUI with PyQt5 that prompts the user input information and
create segment guiding override command(s) and written as segment
override files in ``segment_guiding.py``.
Authors
-------
    - Lauren Chambers
    - Shannon Osborne


Notes
-----
For the GUI to run successfully, the QtAgg matplotlib backend should
be used. This can be set by declaring:
    ::
    import matplotlib
    matplotlib.use('Qt5Agg')

Note that this declaration must occur before pyplot or any other
matplotlib-dependent packages are imported.

"""

# Standard Library Imports
from __future__ import unicode_literals
import glob
import logging
import os

# Third Party Imports
from astropy import units as u
from astropy.coordinates import SkyCoord
from PyQt5 import uic
from PyQt5.QtCore import Qt
from PyQt5.QtWidgets import (QDialog, QMessageBox, QWidget)

# Local Imports
from jwst_magic.utils.coordinate_transforms import nrcpixel_offset_to_v2v3_offset

# Paths
__location__ = os.path.realpath(os.path.join(os.getcwd(), os.path.dirname(__file__)))


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
    program_id : int or str
        APT program number
    observation_num : optional, int or str
        Observation number
    visit_num : optional, int or str
        Visit number
    ra :  optional, float
        RA of guide star
    dec :  optional, float
        DEC of guide star
    threshold : optional, float
        Count rate uncertainty factor
    detector : optional, str
        NIRCam detector used for default combobox

    Returns
    -------
    tup
        Tuple containing the following arguments: (guide_star_params_dict,
        program_id, observation_num, visit_num, threshold_factor,
        countrate_factor)
    """
    def __init__(self, override_type, guider, program_id, observation_num, visit_num, ra=None, dec=None,
                 log=None, threshold=None, detector=None):
        # Initialize attributes
        self.override_type = override_type
        self.guider = guider
        self.program_id = program_id
        self.observation_num = observation_num
        self.visit_num = visit_num
        self.ra = ra
        self.dec = dec
        self.detector = detector if detector is not None else 'A3'
        threshold = None if threshold == '' else threshold

        # Start logger
        if log is None:
            self.log = logging.getLogger(__name__)
        else:
            self.log = log

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

        # Setting only for SOF, not POF
        try:
            self.lineEdit_RA.setText(str(ra if ra is not None else ''))
            self.lineEdit_Dec.setText(str(dec if dec is not None else ''))
            self.lineEdit_countrateUncertainty.setText(str(threshold if threshold is not None else 0.6))
            index = self.comboBox_detector.findText(f'NRC{self.detector}', Qt.MatchFixedString)
            self.comboBox_detector.setCurrentIndex(index)
        except AttributeError:
            pass

    def get_dialog_parameters(self):
        """Parses the user input into the segment guiding dialog box, differentiating
        between input for SOFs and POFs.

        Returns
        -------
        tup
            Tuple containing the following arguments: (guide_star_params_dict,
            program_id, observation_num, visit_num, threshold_factor,
            countrate_factor, countrate_uncertainty_factor)
        """

        # Get parameters for dictionary from dialog
        if self.override_type == "SOF":
            # Parse what the boresight offset is
            if self.radioButton_boresightNIRCam.isChecked():
                detector = str(self.comboBox_detector.currentText())
                x_offset = float(self.lineEdit_boresightX.text())
                y_offset = float(self.lineEdit_boresightY.text())
                v2_offset, v3_offset = nrcpixel_offset_to_v2v3_offset(x_offset, y_offset,
                                                                      detector=detector)
                self.log.info(
                    'Segment Guiding: Applying boresight offset of {}, {} arcsec (Converted from {}, {} {} pixels)'.
                        format(v2_offset, v3_offset, x_offset, y_offset, detector)
                )
            else:
                v2_offset = float(self.lineEdit_boresightV2.text())
                v3_offset = float(self.lineEdit_boresightV3.text())
                self.log.info(
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
                'center_of_pointing': 0
            }

            # Countrate factors
            threshold_factor = float(self.lineEdit_countrateUncertainty.text())
            countrate_factor = None
            countrate_uncertainty_factor = None

        elif self.override_type == "POF":
            countrate_factor = float(self.doubleSpinBox_countrateFactor.value())
            countrate_uncertainty_factor = float(self.doubleSpinBox_countrateUncertaintyFactor.value())
            threshold_factor = None
            guide_star_params_dict = None

        # Get APT information and other necessary parameters
        program_id = self.lineEdit_programNumber.text()
        observation_num = self.lineEdit_observationNumber.text()
        visit_num = self.lineEdit_visitNumber.text()

        return guide_star_params_dict, program_id, observation_num, visit_num, threshold_factor, countrate_factor, countrate_uncertainty_factor


def check_override_overwrite(out_dir, program_id, observation_num, visit_num,
                             logger=None):
    """Check if there is an existing override file with the same program ID,
    observation number, and visit number. If yes, raise a dialog box to prompt
    the user whether or not to overwrite the existing file.

    Parameters
    ----------
    out_dir : str
        Location of out/ directory. If not specified, will be placed
        within the repository: .../jwst_magic/out/
    program_id : int
        APT program number
    observation_num : int
        Observation number
    visit_num : int
        Visit number
    logger : logging object, optional
        The desired logger to write out messages to

    Returns
    -------
    overwrite_existing_file : boolean
        User's response whether or not to overwrite the existing file.
    """
    plural_obs_num = ',' in str(observation_num) or '-' in str(observation_num)

    # Handle null values
    if observation_num == '':
        observation_num = 0
    if visit_num == '':
        visit_num = 0

    # If there is more than one observation, set the visit to 1
    if plural_obs_num:
        visit_num = 1

    # Go through every override file in the output directory and check if it matches
    existing_files = glob.glob(os.path.join(out_dir, '*gs_override*.txt'))
    for file_path in existing_files:

        file_name = os.path.basename(file_path)
        # Don't check reports
        if 'REPORT' in file_name:
            continue

        # Handle the fact that there may or may not be obs/visit numbers
        file_root = file_name.split('.')[0]
        date, _, _, *ids = file_root.split('_')
        if len(ids) >= 3:
            prog, obs, visit, *_ = ids
        elif len(ids) >= 2:
            prog, obs, *_ = ids
            visit = 0
        elif len(ids) >= 1:
            prog, *_ = ids
            obs = 0
            visit = 0

        if isinstance(obs, str):
            plural_obs_num_match = ',' in obs or '-' in obs
        else:
            plural_obs_num_match = None

        # Determine if the current file matches the new one
        # Compare integers if only 1 (or no) obs
        if not plural_obs_num and not plural_obs_num_match:
            match = (int(prog) == int(program_id) and
                    int(obs) == int(observation_num) and
                    int(visit) == int(visit_num))
        # Compare the observation string if more than one obs specified
        else:
            match = (int(prog) == int(program_id) and
                     obs == observation_num and
                     int(visit) == int(visit_num))

        # If there is a match, ask the user whether to overwrite
        if match:
            buttonReply = QMessageBox.question(
                QWidget(), 'Existing Override File',
                "A file already exists from {} with the same program, ".format(date) + \
                "observation, and visit numbers: \n{} \n".format(file_path) + \
                "Do you want to write another?",
                QMessageBox.Yes | QMessageBox.No, QMessageBox.No
            )
            if buttonReply == QMessageBox.Yes:
                logger.info('Segment Guiding: Overwriting file at {}.'.format(file_path))
                return False
            else:
                logger.info(
                    'Segment Guiding: User chose not to overwrite existing file at {}.'.format(
                        file_path))
                return True
