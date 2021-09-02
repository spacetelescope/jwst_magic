"""GUI interface to JWST MAGIC

The primary interface to JWST MAGIC, from which almost
all of the tool functions can be operated.

Authors
-------
    - Keira Brooks
    - Lauren Chambers
    - Shannon Osborne

Use
---
This GUI can be run in the python shell as such:
    ::
    import jwst_magic
    jwst_magic.run_tool_GUI()

References
----------
Standard output/error stream handling:
    - https://stackoverflow.com/questions/8356336/how-to-capture-output-
      of-pythons-interpreter-and-show-in-a-text-widget
    - https://stackoverflow.com/questions/616645/how-do-i-duplicate-sys-
      stdout-to-a-log-file-in-python

Notes
-----
1. For the GUI to run successfully, the QtAgg matplotlib backend should
be used. This can be set by declaring:
    ::
    import matplotlib
    matplotlib.use('Qt5Agg')

Note that this declaration must occur before pyplot or any other
matplotlib-dependent packages are imported.

2. Because this code is run in a suite that also uses pyplot, there
will already by instances of the QApplication object floating around
when this GUI is called. However, only one instance of QApplication can
be run at once without things crashing terribly. In all GUIs within the
JWST MAGIC package, be sure to use the existing instance
of QApplication (access it at QtCore.QCoreApplication.instance()) when
calling the QApplication instance to run a window/dialog/GUI.
"""

# Standard Library Imports
import fnmatch
import glob
import io
import logging
import os
import re
import shutil
import sys
import urllib.request
import zipfile

# Third Party Imports
from astropy.io import ascii as asc
import fgscountrate
from lxml import etree
import matplotlib
import numpy as np
if matplotlib.get_backend() != 'Qt5Agg':
    matplotlib.use('Qt5Agg')  # Make sure that we are using Qt5
from PyQt5 import QtCore, uic, QtGui
from PyQt5.QtCore import Qt, pyqtSlot
from PyQt5.QtWidgets import (QApplication, QMainWindow, QMessageBox, QFileDialog,
                             QDialog, QComboBox, QInputDialog)
from PyQt5.QtGui import QStandardItemModel

# Local Imports
from jwst_magic import run_magic
from jwst_magic.convert_image import renormalize, background_stars_GUI
from jwst_magic.fsw_file_writer import rewrite_prc
from jwst_magic.segment_guiding import segment_guiding
from jwst_magic.star_selector.SelectStarsGUI import StarClickerMatplotlibCanvas, run_SelectStars
from jwst_magic.utils import utils

# Define all needed paths
__location__ = os.path.realpath(os.path.join(os.getcwd(), os.path.dirname(__file__)))
PACKAGE_PATH = os.path.dirname(os.path.realpath(__file__))
OUT_PATH = os.path.split(PACKAGE_PATH)[0]  # Location of out/ and logs/ directory
SOGS_PATH = '***REMOVED***/guiding/'

LOGGER = logging.getLogger(__name__)

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# OUTPUT HANDLER
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


class EmittingStream(QtCore.QObject):
    """Define a class to tee sys.stdout to textEdit_log
    """

    def __init__(self, textEdit_log):
        # Initialize the super class (QObject)
        super(EmittingStream, self).__init__()

        # Redefine sys.stdout to call this class
        self.stdout = sys.stdout
        sys.stdout = self

        # Connect to the main window's textEdit_log widget
        self.textEdit_log = textEdit_log

    def __del__(self):
        """Upon deletion, set sys.stdout back to normal
        """
        sys.stdout = self.stdout

    def write(self, text):
        """Write output to both stdout and to textEdit_log
        """
        self.stdout.write(text)
        try:
            if 'warning' in text.lower():
                text = "<span style=\" font-size:13pt; font-weight:600; color:#B8AE17;\" >" + text + "</span>"
            if 'error' in text.lower():
                text = "<span style=\" font-size:13pt; font-weight:600; color:#B84317;\" >" + text + "</span>"
            self.write_output(text)
        except RuntimeError:
            pass

    def flush(self):
        """In reality, empty a buffer. Here, do nothing.
        """
        pass

    def write_output(self, text):
        """Format and append stdout text to the log QTextEdit.
        """
        current_text = self.textEdit_log.toHtml()
        text_noANSI = re.sub(r'(\[.+?m)', '', text).replace('\n', '<br>')
        self.textEdit_log.setHtml(current_text + text_noANSI)

        # Move the cursor to the end so that the most recent text is visible
        cursor = self.textEdit_log.textCursor()
        cursor.movePosition(QtGui.QTextCursor.End)
        self.textEdit_log.setTextCursor(cursor)

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# GUI CLASS DEFINITION
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


class MasterGui(QMainWindow):
    def __init__(self, root=None, norm_value=12.0, norm_units='FGS Magnitude',
                 nircam_det=None, nircam=True, smoothing='default',
                 steps=None, in_file=None, bkgd_stars=False, out_dir='', convert_im=True,
                 star_selection=True, star_selection_gui=True, file_writer=True,
                 segment_guiding=False, app=None, itm=False):

        # Initialize attributes
        self.app = app
        self.commissioning_dict = {}
        self.root_default = root
        self.converted_im_circles = []
        self.shifted_im_circles = []
        self.bkgd_stars = None
        self.bkgrdstars_hdr = None
        self._bkgdstars_dialog = None
        self.itm = itm
        self.program_id = ''
        self.observation_num = ''
        self.visit_num = ''
        self.gs_id = ''
        self.apt_guider = ''
        self.gs_ra = ''
        self.gs_dec = ''
        self._test_sg_dialog = None
        self.log = None
        self.log_filename = None

        # Initialize main window object
        QMainWindow.__init__(self)

        # Import .ui file
        uic.loadUi(os.path.join(__location__, 'masterGUI.ui'), self)

        # Set up the custom output stream
        EmittingStream(self.textEdit_log)

        # Create and load GUI session
        self.setWindowTitle('JWST MAGIC')
        self.adjust_screen_size_mainGUI()
        self.init_matplotlib()
        self.define_MainGUI_connections()
        self.setup_commissioning_naming()
        self.show()

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # GUI CONSTRUCTION
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    def adjust_screen_size_mainGUI(self):
        """ Adjust the GUI sizes for laptop screens
        """
        # Determine screen size
        screen_size = self.app.desktop().screenGeometry()
        width, height = screen_size.width(), screen_size.height()

        # Adjust the scroll window size
        if width - 200 < self.scrollArea_mainGUI.minimumWidth():
            # Window is too wide
            self.scrollArea_mainGUI.setMinimumWidth(width - 200)
        if height - 200 < self.scrollArea_mainGUI.minimumHeight():
            # Window is too tall
            self.scrollArea_mainGUI.setMinimumHeight(height - 200)

    def init_matplotlib(self):
        """Set up the three matplotlib canvases that will preview the
        input, converted, and shifted images in the "Image Preview" section.
        """

        # Connect matplotlib canvases to the tab widgets
        self.canvas_input = StarClickerMatplotlibCanvas(
            parent=self.canvas_input_slot, data=None, x=None, y=None,
            left=0.12, bottom=0, right=0.87, top=1
        )
        self.canvas_converted = StarClickerMatplotlibCanvas(
            parent=self.canvas_converted_slot, data=None, x=None, y=None,
            left=0.12, bottom=0, right=0.87, top=1
        )
        self.canvas_shifted = StarClickerMatplotlibCanvas(
            parent=self.canvas_shifted_slot, data=None, x=None, y=None,
            left=0.12, bottom=0, right=0.87, top=1
        )

        # Set the dimensions to be square and as big as possible
        max_dim = max(self.canvas_converted_slot.width(), self.canvas_converted_slot.height())
        self.canvas_input_slot.setMinimumSize(max_dim, max_dim)
        self.canvas_input.setMinimumSize(max_dim, max_dim)
        self.canvas_converted_slot.setMinimumSize(max_dim, max_dim)
        self.canvas_converted.setMinimumSize(max_dim, max_dim)
        self.canvas_shifted_slot.setMinimumSize(max_dim, max_dim)
        self.canvas_shifted.setMinimumSize(max_dim, max_dim)

    def define_MainGUI_connections(self):
        # Standard output and error
        self.textEdit_log.setFont(QtGui.QFont("Courier New"))
        self.buttonGroup_name.buttonClicked.connect(self.update_naming_method)

        # Main window widgets
        self.pushButton_run.clicked.connect(self.run_tool)
        self.pushButton_quit.clicked.connect(self.close_application)

        # General input widgets
        self.pushButton_inputImage.clicked.connect(self.update_input)
        self.lineEdit_inputImage.editingFinished.connect(self.update_input)
        self.buttonGroup_guider.buttonClicked.connect(self.update_filepreview)
        self.buttonGroup_guider.buttonClicked.connect(self.check_guider_against_apt)
        self.lineEdit_root.editingFinished.connect(self.update_filepreview)
        self.pushButton_root.clicked.connect(self.on_click_root)
        self.pushButton_defaultout.clicked.connect(self.on_click_default_out)
        self.pushButton_out.clicked.connect(self.on_click_out)
        self.textEdit_out.installEventFilter(self)
        self.pushButton_manualid.clicked.connect(self.update_apt_gs_values)

        # Image converter widgets
        self.pushButton_backgroundStars.clicked.connect(self.on_click_bkgdstars)
        self.pushButton_delbackgroundStars.clicked.connect(self.on_click_del_bkgrdstars)
        self.horizontalSlider_coarsePointing.sliderReleased.connect(self.on_change_jitter)
        self.lineEdit_coarsePointing.editingFinished.connect(self.on_change_jitter)

        # Star selector widgets
        self.pushButton_regfileStarSelector.clicked.connect(self.on_click_infile)
        self.comboBox_regfileStarSelector.activated.connect(self.on_combobox_change)

        # FSW File writer widgets
        self.groupBox_fileWriter.toggled.connect(self.update_groupBox_fileWriter)
        self.checkBox_OSS.toggled.connect(self.turn_off_threshold)

        # Segment guiding widgets
        self.buttonGroup_segmentGuiding_idAttitude.buttonClicked.connect(self.update_segment_guiding_shift)
        self.groupBox_segmentGuiding.toggled.connect(self.update_segment_guiding_shift)
        self.radioButton_regfileSegmentGuiding.toggled.connect(self.enable_segment_guiding)
        self.radioButton_photometryOverride.toggled.connect(self.enable_segment_guiding)
        self.lineEdit_regfileSegmentGuiding.editingFinished.connect(self.update_checkable_combobox)

        # Image preview widgets
        self.checkBox_showStars.toggled.connect(self.on_click_showstars)
        self.checkBox_showStars_shifted.toggled.connect(self.on_click_showstars)
        self.comboBox_showcommandsconverted.currentIndexChanged.connect(self.update_converted_image_preview)
        self.comboBox_showcommandsshifted.currentIndexChanged.connect(self.update_shifted_image_preview)

    def setup_commissioning_naming(self):
        # If not on SOGS:
        if not utils.on_sogs_network():
            # Shut down selection boxes
            self.comboBox_practice.removeItem(0)
            self.comboBox_practice.setDisabled(True)
            self.comboBox_practice.addItem('* NOT ON SOGS NETWORK *')

            # Switch to the manual naming pane
            self.stackedWidget.setCurrentIndex(1)
            self.radioButton_name_manual.setChecked(True)
            return

        # If on SOGS, pull out practice names from existing practice directories
        else:
            sogs_search = os.path.join(SOGS_PATH, '*')
            sogs_dirs = [os.path.basename(dir) for dir in glob.glob(sogs_search)]
            for d in ['data', 'processing', 'MAGIC_logs']:
                sogs_dirs.remove(d)

        # Load all OTE cars from commissioning activities YAML file
        self.commissioning_dict = utils.get_car_data()
        cars_list = list(self.commissioning_dict.keys())

        # Use to populate practice and CAR dropdown boxes
        for practice in sorted(sogs_dirs):
            self.comboBox_practice.addItem(practice)
        for car in sorted(cars_list):
            self.comboBox_car.addItem(car.upper())

        # Connect combo boxes to one another
        self.comboBox_practice.currentIndexChanged.connect(self.update_commissioning_name)
        self.comboBox_car.currentIndexChanged.connect(self.update_commissioning_name)
        self.lineEdit_obs.editingFinished.connect(self.update_commissioning_name)
        self.pushButton_commid.clicked.connect(self.update_commissioning_name)

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # WIDGET CONNECTIONS
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    @pyqtSlot()
    def eventFilter(self, source, event):
        """"EdittingFinished" event filter for out directory textbox
        """
        if hasattr(self, 'textEdit_out'):
            # Parse out the focus leaving the out_dir box, and update
            # all other textboxes accordingly
            if event.type() == QtCore.QEvent.FocusOut:
                # (event.type() == QtCore.QEvent.KeyPress and event.key() == QtCore.Qt.Key_Return):

                # Read the current out directory name
                dirname = self.textEdit_out.toPlainText()

                # Remove any new lines from the dirname
                dirname = dirname.rstrip()
                self.textEdit_out.setPlainText(dirname)

                # Only continue if that directory actually exists
                if dirname is not None:
                    if not os.path.exists(dirname):
                        raise FileNotFoundError('Output directory {} does not exist.'.format(dirname))

                self.update_filepreview(new_guiding_selections=True)

        return False

    def close_application(self):
        """ Close the window
        """
        self.close()
        self.app.quit()

    def run_tool(self):
        """
        Takes inputs provided by user and runs run_magic
        """
        # Required
        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        if self.lineEdit_inputImage.text() == "" and not self.radioButton_photometryOverride.isChecked():
            self.no_inputImage_dialog()
            return
        if not self.buttonGroup_guider.checkedButton():
            self.no_guider_dialog()
            return
        if self.lineEdit_root.text() == "" and not self.radioButton_name_commissioning.isChecked():
            self.no_root_dialog()
            return

        # General input
        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        input_image = self.lineEdit_inputImage.text()
        guider = int(self.buttonGroup_guider.checkedButton().text())
        if self.radioButton_name_manual.isChecked():
            root = self.lineEdit_root.text()
            out_dir = self.textEdit_out.toPlainText().rstrip()
        elif self.radioButton_name_commissioning.isChecked():
            root = 'for_obs{:02d}'.format(int(self.lineEdit_obs.text()))
            out_dir = os.path.join(SOGS_PATH,
                                   self.comboBox_practice.currentText(),
                                   self.comboBox_car.currentText().lower().replace('-', ''),
                                   )
        copy_original = True

        # Check for mis-matched guider
        list_of_files = [[os.path.join(dirpath, file) for file in filenames] for (dirpath, _, filenames) in
                         os.walk(os.path.join(out_dir, 'out', root))]
        list_of_files = [item for sublist in list_of_files for item in sublist]
        opposite_guider = [2 if guider == 1 else 1][0]
        if True in [True if ('_G{}_'.format(opposite_guider) in file  or '_G{}.'.format(opposite_guider) in
                file) else False for file in list_of_files]:

            raise ValueError('Data from GUIDER {} found in the root path: {}, which does not match the chosen '
                             'GUIDER {}. Delete all contents from this directory before writing data with a new '
                             'guider.'.format(opposite_guider, os.path.join(out_dir, 'out', root), guider))

        # Convert image
        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        convert_im = True
        nircam = self.radioButton_NIRCam.isChecked()
        nircam_det = str(self.comboBox_detector.currentText())
        normalize = self.checkBox_normalize.isChecked()
        norm_unit = self.comboBox_normalize.currentText()
        try:
            norm_value = float(self.lineEdit_normalize.text())
        except ValueError:
            norm_value = self.lineEdit_normalize.text()
        coarse_point = self.checkBox_coarsePointing.isChecked()
        jitter_rate_arcsec = float(self.lineEdit_coarsePointing.text())
        bkgd_stars = self.bkgd_stars
        bkgrdstars_hdr = self.bkgrdstars_hdr
        itm = self.itm

        # Handle the case where we want to use a pre-existing converted image
        pre_existing_im = self.checkBox_useConvertedImage.isChecked() and convert_im
        if pre_existing_im:
            convert_im = False
            input_image = self.converted_im_file
            copy_original = False
        if self.groupBox_segmentGuiding.isChecked() and self.radioButton_photometryOverride.isChecked():
            convert_im = False

        if convert_im:
            # Check normalization information is defined
            if normalize and norm_value == '':
                raise ValueError('Image Normalization box checked, but no value given.')

       # Set smoothing
        if self.checkBox_globalAlignment.isChecked():
            smoothing = 'high'
        elif self.checkBox_noSmoothing.isChecked():
            smoothing = 'low'
        else:
            smoothing = 'default'

        # Star selection
        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        star_selection = self.groupBox_starSelector.isChecked()
        star_selectiongui = self.radioButton_starSelectorGUI.isChecked()
        if not self.radioButton_regfileStarSelector.isChecked():
            in_file = None
        else:
            in_file = [self.comboBox_regfileStarSelector.itemText(i) for i in
                       range(self.comboBox_regfileStarSelector.count())]

            if len(in_file) == 0 and not self.groupBox_segmentGuiding.isChecked():
                raise ValueError('"Load Unshifted Data from File" option chosen for Star Selector section, '
                                 'but no files were chosen. Either choose files to be loaded or switch to the '
                                 'Click-to-Select GUI option.')

        # File writer
        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        file_writer = self.groupBox_fileWriter.isChecked()

        if file_writer and not star_selection and not self.groupBox_segmentGuiding.isChecked():
            raise ValueError('Cannot run FSW File Writer section on its own. Need to define the input guiding '
                             'selections files in the Star Selector section either through the GUI or loading '
                             'files.')

        steps = []
        if self.checkBox_CAL.isChecked():
            steps.append('CAL')
        if self.checkBox_ID.isChecked():
            steps.append('ID')
        if self.checkBox_ACQ1.isChecked():
            steps.append('ACQ1')
        if self.checkBox_ACQ2.isChecked():
            steps.append('ACQ2')
        if self.checkBox_TRK.isChecked():
            steps.append('TRK')
        if self.checkBox_LOSTRK.isChecked():
            steps.append('LOSTRK')

        # Check threshold type can be a float
        try:
            threshold = float(self.lineEdit_threshold.text())
        except ValueError:
            raise ValueError('Threshold value must be a float')

        # Shift image to ID attitude:
        shift_id_attitude = self.checkBox_id_attitude.isChecked()

        # For the POF case, create DHAS files using the default OSS numbers
        use_oss_defaults = self.checkBox_OSS.isChecked()

        # Rewrite .prc and guiding_selections*.txt ONLY
        if self.checkBox_rewritePRC.isChecked():
            # If use_oss_numbers if also checked, raise an error (needs catalog countrate information)
            if self.checkBox_OSS.isChecked():
                raise ValueError('Cannot use Default OSS Numbers and Rewrite PRC functionality at the same time.')

            # Open converted FGS file
            data, _ = utils.get_data_and_header(self.converted_im_file)

            # Update array for showing in LogNorm
            data[data <= 0] = 1

            # Open all_found_psfs*.txt (list of all identified segments)
            all_psfs = self.all_found_psfs_file
            out_path = os.path.join(out_dir, 'out', self.lineEdit_root.text())

            if self.all_found_psfs_file != '':
                all_rows = asc.read(all_psfs)
                x = all_rows['x'].data
                y = all_rows['y'].data
            else:
                raise FileNotFoundError('Cannot find an all found PSFs file in this directory {}. '
                             'This file can be created by running the star selection section '
                             'of the GUI.'.format(out_path))

            # Run the select stars GUI to determine the new orientation
            inds_list, center_of_pointing = run_SelectStars(data, x, y, 20, guider,
                                                out_dir=out_path,
                                                print_output=False,
                                                masterGUIapp=self.app)

            # Print indices of each guiding configuration
            for i in range(len(inds_list)):
                ind = inds_list[i]
                LOGGER.info('Master GUI: Guiding Configuration {} - GS = {}, RS = {}'.format(i, ind[0],
                                                                                 ', '.join([str(c) for c in ind[1:]])))

            # Rewrite the id.prc and acq.prc files
            rewrite_prc.rewrite_prc(inds_list, center_of_pointing, guider, root, out_dir, threshold=threshold,
                                    shifted=shift_id_attitude)

            # Update converted image preview
            self.update_filepreview(new_guiding_selections=True)
            return

        # Segment guiding
        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        if self.groupBox_segmentGuiding.isChecked():
            # If Id/Obs/Visit and/or GS info hasn't already been set by parse_header() or commissioning section
            if not all(hasattr(self, attr) for attr in ["program_id", "observation_num", "visit_num"]):
                self.program_id, self.observation_num, self.visit_num = '', '', ''
            if not all(hasattr(self, attr) for attr in ["gs_id", "gs_ra", "gs_dec"]):
                self.gs_id, self.gs_ra, self.gs_dec = '', '', ''

            # Check if this is a photometry only override file or segment override file
            if self.radioButton_photometryOverride.isChecked():
                # Raise error if normalization information isn't set
                if norm_value == '':
                    raise ValueError('Missing a normalization value (guide star ID, count rate, or magnitude) which '
                                     'is required to make a photometry override file')

                # Initialize the dialog
                self._test_sg_dialog = segment_guiding.SegmentGuidingGUI.SegmentGuidingDialog(
                                       "POF", None, self.program_id, self.observation_num, self.visit_num, log=None
                )
                # Generate the file
                segment_guiding.generate_photometry_override_file(
                    root, self.program_id, self.observation_num, self.visit_num, guider=guider,
                    norm_value=norm_value, norm_unit=norm_unit, out_dir=out_dir,
                    parameter_dialog=True, dialog_obj=self._test_sg_dialog, log=LOGGER
                )

            else:
                # Define location of all_found_psfs catalog file(s)
                if self.radioButton_shifted.isChecked():
                    guiding_files = self.shifted_guiding_selections_file_list
                    all_psf_files = self.shifted_all_found_psfs_file_list
                    center_pointing_files = [file.replace('guiding_selections', 'center_pointing')
                                             for file in self.shifted_guiding_selections_file_list]
                else:
                    guiding_files = self.guiding_selections_file_list
                    all_psf_files = [self.all_found_psfs_file] * len(guiding_files)
                    center_pointing_files = [self.all_found_psfs_file.replace('unshifted_all_found_psfs_',
                                                                              'center_pointing_')] * len(guiding_files)

                LOGGER.info('Master GUI: Pulling center of pointing information from same location as guiding files')

                # Load selected guiding_selections*.txt
                if len(self.comboBox_guidingcommands.checkedItems()) == 0:
                    raise ValueError('No guiding commands chosen from dropdown box.')

                combobox_choices = [item.text() for item in self.comboBox_guidingcommands.checkedItems()]
                combobox_filenames = [file.split(': ')[-1] for file in combobox_choices]

                # Re-order the files via dialog box if checked
                if self.checkBox_configorder.isChecked():
                    combobox_filenames = self.change_config_order_dialog(combobox_filenames)

                selected_segs_list = []
                segment_infile_list = []
                center_pointing_list = []
                for file in combobox_filenames:
                    ind = np.where([file in filepath for filepath in guiding_files])[0][0]
                    selected_segs_list.append(guiding_files[ind])
                    segment_infile_list.append(all_psf_files[ind])
                    center_pointing_list.append(center_pointing_files[ind])

                # Verify that the all_found_psfs*.txt file(s) exist
                if False in [os.path.exists(path) for path in segment_infile_list]:
                    raise OSError('Missing all founds psfs file(s) {}.'.format(segment_infile_list))

                # Run the tool and generate the file
                # Initialize the dialog
                self._test_sg_dialog = segment_guiding.SegmentGuidingGUI.SegmentGuidingDialog(
                    "SOF", guider, self.program_id, self.observation_num, self.visit_num,
                    ra=self.gs_ra, dec=self.gs_dec, threshold=float(self.lineEdit_threshold.text()),
                    detector=str(self.comboBox_detector.currentText())[:2], log=None
                )

                # Generate the file
                segment_guiding.generate_segment_override_file(
                    segment_infile_list, guider, self.program_id, self.observation_num,
                    self.visit_num, ra=self.gs_ra, dec=self.gs_dec,
                    root=root, out_dir=out_dir, selected_segs_list=selected_segs_list,
                    center_pointing_list=center_pointing_list,
                    master_gui_app=self.app, parameter_dialog=True,
                    dialog_obj=self._test_sg_dialog, log=LOGGER
                )

            # Update converted image preview
            self.update_filepreview()
            return
        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

        # Check for commissioning naming usage not matching intended use case
        if self.radioButton_name_commissioning.isChecked():
            check_cases = [not (input_image.endswith('_cal.fits') or input_image.endswith('_cal_full_frame.fits')),
                           norm_unit == 'FGS countrate']
            if True in check_cases and pre_existing_im is False:
                if check_cases[0]:
                    errmsg1 = 'read in a file that is not a cal.fits file or a padded trk image ' \
                               '(ending in cal_full_frame.fits)'
                    errmsg2 = 'change the file'
                elif check_cases[1]:
                    errmsg1 = 'set the normalization to FGS Countrate'
                    errmsg2 = 'reset the normalization unit'
                raise ValueError("Given the current inputs, MAGIC cannot be run with commissioning naming. "
                                 "MAGIC's commissioning naming is meant to be used for the nominal commissioning case, "
                                 f"but the user has {errmsg1}. Either {errmsg2}, or switch to manual naming.")

        # Run MAGIC
        if convert_im or star_selection or file_writer:
            run_magic.run_all(input_image, guider,
                              root=root,
                              norm_value=norm_value,
                              norm_unit=norm_unit,
                              nircam_det=nircam_det,
                              nircam=nircam,
                              smoothing=smoothing,
                              steps=steps,
                              guiding_selections_file=in_file,
                              bkgd_stars=bkgd_stars,
                              bkgrdstars_hdr=bkgrdstars_hdr,
                              out_dir=out_dir,
                              convert_im=convert_im,
                              star_selection=star_selection,
                              file_writer=file_writer,
                              masterGUIapp=self.app,
                              copy_original=copy_original,
                              normalize=normalize,
                              coarse_pointing=coarse_point,
                              jitter_rate_arcsec=jitter_rate_arcsec,
                              itm=itm,
                              shift_id_attitude=shift_id_attitude,
                              thresh_factor=threshold,
                              use_oss_defaults=use_oss_defaults,
                              logger_passed=LOGGER,
                              )

            # Update converted image preview
            self.update_filepreview(new_guiding_selections=True)

    def update_groupBox_fileWriter(self):
        """Enable/disable items in FSW group box"""
        if self.sender() == self.groupBox_fileWriter:
            if self.groupBox_fileWriter.isChecked():
                self.groupBox_starSelector.setChecked(True)

    def on_change_jitter(self):
        """If the coarse pointing slider or text box controlling the
        jitter rate are changed, update the other one accordingly.
        """
        slider_range = np.linspace(0, 0.3, 101)
        if self.sender() == self.horizontalSlider_coarsePointing:
            # Get slider value
            slider_value = int(self.horizontalSlider_coarsePointing.value())

            # Calculate matching textbox value
            jitter = slider_range[slider_value]

        elif self.sender() == self.lineEdit_coarsePointing:
            # Get jitter textbox value
            jitter = float(self.lineEdit_coarsePointing.text())

            # Make sure the input is not out of bounds
            jitter = min(jitter, 0.3)
            jitter = max(jitter, 0)

            # Calculate matching slider value
            slider_value = (np.abs(slider_range - jitter)).argmin()

        # Update both to show the same value
        self.horizontalSlider_coarsePointing.setValue(slider_value)
        self.lineEdit_coarsePointing.setText('{:.3f}'.format(jitter))

    def update_input(self):
        """ Using the Input Image Open button (open file) and textbox
        """
        # If coming from the button, open the dialog
        if self.sender() == self.pushButton_inputImage:
            # Read selected filename
            filename = self.open_filename_dialog("NIRCam or FGS image", file_type="FITS files (*.fits)")
            self.lineEdit_inputImage.setText(filename)

        # If coming from the textbox, just read the new value
        elif self.sender() == self.lineEdit_inputImage:
            # Read selected filename
            filename = self.lineEdit_inputImage.text()

        # Only continue if the entered image path actually exists:
        if filename is not None and filename != '':
            if not os.path.exists(filename):
                raise FileNotFoundError('Input image {} does not exist.'.format(filename))

        # Derive the root from the filename and confirm a blank output directory
        if self.radioButton_name_manual.isChecked():
            if self.textEdit_out.toPlainText() == "":
                self.textEdit_out.setEnabled(True)
            if self.textEdit_out.toPlainText() == '':
                self.textEdit_out.setEnabled(True)

        # Update the example filepath (and converted image preview, if possible)
        self.update_filepreview()

        # Show input image preview
        self.load_input_image_data(filename)

        # Show converted and shifted image preview, if possible
        self.update_converted_image_preview()
        self.update_shifted_image_preview()

        return filename

    def on_click_out(self):
        """ Using the Out Dir Open button (open directory)
        """
        # Open the Finder directory dialog
        dirname = self.open_dirname_dialog()

        # Only continue if that directory actually exists
        if dirname is not None:
            # Remove any new lines from the dirname
            dirname = dirname.rstrip()
            self.textEdit_out.setText(dirname)

            if not os.path.exists(dirname):
                raise FileNotFoundError('Output directory {} does not exist.'.format(dirname))

            self.update_filepreview(new_guiding_selections=True)

    def on_click_infile(self):
        """ Using the Infile Open button (open file) """
        # Determine which infile is being edited
        if self.sender() == self.pushButton_regfileStarSelector:
            filename_list = self.open_filename_dialog('In/Reg file(s)', multiple_files=True,
                                                 file_type="Input file (*.txt *.incat);;All files (*.*)")

            self.update_regfile_starselector_combobox(filename_list)
            self.update_guiding_selections(new_selections=filename_list)

        return filename_list

    def on_click_root(self):
        """ Using the Set Default Root button"""
        root = utils.make_root(None, self.lineEdit_inputImage.text())
        self.lineEdit_root.setText(root)
        self.update_filepreview()
        return root

    def on_click_default_out(self):
        """ Using the Set Default Out Dir button"""
        out_dir = OUT_PATH
        self.textEdit_out.setText(out_dir)
        self.update_filepreview()
        return out_dir

    def on_click_bkgdstars(self):
        if self.lineEdit_inputImage.text() == "":
            self.no_inputImage_dialog()
            return

        if not self.buttonGroup_guider.checkedButton():
            self.no_guider_dialog()
            return

        if self.radioButton_name_manual.isChecked():
            root = self.lineEdit_root.text()
            out_dir = self.textEdit_out.toPlainText().rstrip()
        elif self.radioButton_name_commissioning.isChecked():
            root = 'for_obs{:02d}'.format(int(self.lineEdit_obs.text()))
            out_dir = os.path.join(SOGS_PATH,
                                   self.comboBox_practice.currentText(),
                                   self.comboBox_car.currentText().lower().replace('-', ''),
                                   )

        # Enable the textbox
        self.textEdit_backgroundStars.setEnabled(True)

        guider = int(self.buttonGroup_guider.checkedButton().text())

        # If normalization is turned on, read the normalization value & unit
        # and calculate FGS Mag of the guide star
        if self.checkBox_normalize.isChecked():
            try:
                norm_value = float(self.lineEdit_normalize.text())
            except ValueError:
                norm_value = self.lineEdit_normalize.text()
            norm_unit = self.comboBox_normalize.currentText()
            fgs_countrate, fgs_mag = renormalize.convert_to_countrate_fgsmag(norm_value, norm_unit, guider)

        # If not, determine the FGS counts of the input image
        else:
            input_image = self.lineEdit_inputImage.text()
            data, _ = utils.get_data_and_header(input_image)
            fgs_countrate = np.sum(data[data > np.median(data)])
            fgs_mag = fgscountrate.convert_cr_to_fgs_mag(fgs_countrate, guider)

        # Run background stars window
        self._bkgdstars_dialog = background_stars_GUI.BackgroundStarsDialog(guider, fgs_mag,
                                                                            out_dir=out_dir, root=root,
                                                                            ra=self.gs_ra, dec=self.gs_dec,
                                                                            in_master_GUI=True)
        accepted = self._bkgdstars_dialog.exec()

        # Pull dict of (x,y) and mag values for each star
        self.bkgd_stars = self._bkgdstars_dialog.return_dict() if accepted else None
        if self.bkgd_stars is None:
            if accepted:  # click ok without finishing the process of adding background stars
                raise ValueError('Background Stars GUI missing information. No background stars selected')
            else:  # click cancel
                # delete any file saved out
                bkgrd_image = os.path.join(out_dir, 'out', root, 'background_stars_{}_G{}.png'.format(root, guider))
                if os.path.exists(bkgrd_image):
                    os.remove(bkgrd_image)
                raise ValueError('Background Stars GUI closed. No background stars selected')

        # Record the method used to generate the background stars and populate the main GUI with that
        method = self._bkgdstars_dialog.method
        method_adverb = {'random': 'randomly',
                         'user-defined': 'as defined by the user',
                         'catalog': 'from a GSC query'}
        self.bkgrdstars_hdr = {}
        self.bkgrdstars_hdr['BACKMETH'] = method
        if method == 'catalog':
            self.bkgrdstars_hdr['BACK_RA'] = self._bkgdstars_dialog.ra_gs
            self.bkgrdstars_hdr['BACK_DEC'] = self._bkgdstars_dialog.dec_gs
            self.bkgrdstars_hdr['BACK_PA'] = self._bkgdstars_dialog.position_angle

        if isinstance(self.bkgd_stars, dict) and method is not None:
            self.textEdit_backgroundStars.setText('{} background stars added {}'.
                                                  format(len(self.bkgd_stars['x']), method_adverb[method]))

    def on_click_del_bkgrdstars(self):
        """Reset background stars information"""
        self.bkgd_stars = None
        self.bkgrdstars_hdr = None
        self.textEdit_backgroundStars.setText('No background stars added')

    def on_click_showstars(self, show):
        """Show or hide plots of star positions and selected stars.
        """
        if self.sender() == self.checkBox_showStars:
            for line in self.converted_im_circles:
                line[0].set_visible(show)
            self.canvas_converted.peaks.set_visible(show)

            self.canvas_converted.draw()

        elif self.sender() == self.checkBox_showStars_shifted:
            for line in self.shifted_im_circles:
                line[0].set_visible(show)
            self.canvas_shifted.peaks.set_visible(show)

            self.canvas_shifted.draw()

    def on_combobox_change(self):
        """
        Clicking comboBox_regfileStarSelector doesn't do anything
        """
        if self.sender() == self.comboBox_regfileStarSelector:
            self.comboBox_regfileStarSelector.setCurrentIndex(0)

    def toggle_convert_im(self):
        # TODO: it's unclear why I (KJB) set this to false. but we want to be able
        # to normalize ITM images so I have commented this out until I better understand it
        # if self.itm:
        #     self.checkBox_normalize.setEnabled(False)
        pass

    def update_segment_guiding_shift(self):
        """This gets triggered when you select the segment guiding
        box or when you click one of the id attitude radio buttons"""
        if self.groupBox_segmentGuiding.isChecked():
            if self.radioButton_name_manual.isChecked():
                root_dir = os.path.join(self.textEdit_out.toPlainText(), 'out', self.lineEdit_root.text())
            else:
                root_dir = self.textEdit_name_preview.toPlainText()

            self.lineEdit_regfileSegmentGuiding.setText(root_dir)

            self.update_checkable_combobox()

    def update_naming_method(self):
        if self.radioButton_name_manual.isChecked():
            self.stackedWidget.setCurrentIndex(1)
            if self.textEdit_out.toPlainText() == "":
                self.textEdit_out.setEnabled(True)
                self.textEdit_out.setText('')
        elif self.radioButton_name_commissioning.isChecked():
            self.stackedWidget.setCurrentIndex(0)

        self.update_filepreview()

    def check_guider_against_apt(self):
        if self.apt_guider != '' and self.buttonGroup_guider.checkedButton() is not None:
            # Check the guider in the APT file matches what's chosen in the GUI
            if int(self.apt_guider) != int(self.buttonGroup_guider.checkedButton().text()):
                self.mismatched_apt_guider_dialog()

    def update_apt_gs_values(self):
        # Query APT + call FGSCountrate tool for Guide Star Information
        if self.radioButton_name_manual.isChecked():
            if self.lineEdit_manualid.text() == '' and self.lineEdit_manualobs.text() == '':
                self.program_id, self.observation_num, self.visit_num = '', '', ''
                self.gs_id, self.gs_ra, self.gs_dec = '', '', ''
            elif self.lineEdit_manualid.text() != '' and self.lineEdit_manualobs.text() != '':
                self.program_id = int(self.lineEdit_manualid.text())
                self.observation_num = int(self.lineEdit_manualobs.text())
                self.visit_num = 1  # Will we ever have a visit that's not 1?
                self.gs_id, self.apt_guider, self.gs_ra, self.gs_dec = self.query_apt_for_gs(self.program_id,
                                                                                             self.observation_num)
            else:
                raise ValueError('Must set both program ID and observation number to use APT')
        elif self.radioButton_name_commissioning.isChecked():
            self.program_id = int(self.lineEdit_commid.text())
            self.observation_num = int(self.lineEdit_obs.text())
            self.visit_num = 1  # Will we ever have a visit that's not 1?
            self.gs_id, self.apt_guider, self.gs_ra, self.gs_dec = self.query_apt_for_gs(self.program_id,
                                                                                         self.observation_num)

        # Check the guider in the APT file matches what's chosen in the GUI
        self.check_guider_against_apt()

        # Update GSID in image normalization
        self.lineEdit_normalize.setText(str(self.gs_id))
        index = self.comboBox_normalize.findText('Guide Star ID', Qt.MatchFixedString)
        self.comboBox_normalize.setCurrentIndex(index)

    def update_commissioning_name(self):
        # Check which values have been selected already
        valid_practice = self.comboBox_practice.currentText() != '- Select Practice -'
        valid_car = self.comboBox_car.currentText() != '- Select CAR -'

        # When the CAR step is changed...
        if valid_car and self.sender() == self.comboBox_car:
            # Add/Update the current APT program number
            self.lineEdit_commid.setText(str(self.commissioning_dict[self.comboBox_car.currentText().lower()]))

        valid_obs = self.lineEdit_obs.text() != ""
        valid_all = valid_practice and valid_car and valid_obs

        # Update which boxes are enabled and disabled accordingly
        self.comboBox_car.setEnabled(valid_practice)
        self.lineEdit_commid.setEnabled(valid_practice)
        self.lineEdit_obs.setEnabled(valid_practice & valid_car)
        self.pushButton_commid.setEnabled(valid_practice & valid_car)

        # Update the preview output path
        if valid_all:
            path = os.path.join(SOGS_PATH,
                                self.comboBox_practice.currentText(),
                                self.comboBox_car.currentText().lower().replace('-', ''),
                                'out',
                                'for_obs{:02d}'.format(int(self.lineEdit_obs.text()))
                                )
            self.textEdit_name_preview.setText(path)
        else:
            self.textEdit_name_preview.setText('')

        # Update file previews
        self.update_filepreview()

        # Update population of guide star information
        if valid_all and any([self.sender() == self.comboBox_car, self.sender() == self.lineEdit_obs,
                              self.sender() == self.pushButton_commid]):
            self.update_apt_gs_values()

    def enable_segment_guiding(self):
        """Set up connection to enable/disable segment guiding options between photometry and regfile radio buttons"""
        if self.sender() == self.radioButton_photometryOverride:
            self.lineEdit_regfileSegmentGuiding.setEnabled(False)
            self.label_guidingcommands.setEnabled(False)
            self.comboBox_guidingcommands.setEnabled(False)
            self.checkBox_configorder.setEnabled(False)
            self.radioButton_shifted.setEnabled(False)
            self.radioButton_unshifted.setEnabled(False)
        elif self.sender() == self.radioButton_regfileSegmentGuiding:
            self.lineEdit_regfileSegmentGuiding.setEnabled(True)
            self.label_guidingcommands.setEnabled(True)
            self.comboBox_guidingcommands.setEnabled(True)
            self.checkBox_configorder.setEnabled(True)
            self.radioButton_shifted.setEnabled(True)
            self.radioButton_unshifted.setEnabled(True)

    def turn_off_threshold(self):
        """Disable threshold lineEdit when OSS default box is checked"""
        if self.checkBox_OSS.isChecked():
            self.label_threshold.setEnabled(False)
            self.lineEdit_threshold.setEnabled(False)
        else:
            self.label_threshold.setEnabled(True)
            self.lineEdit_threshold.setEnabled(True)

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # DIALOG BOXES
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    @staticmethod
    def no_inputImage_dialog():
        no_input_im_dialog = QMessageBox()
        no_input_im_dialog.setText('No input image' + ' ' * 50)
        no_input_im_dialog.setInformativeText('The tool will not be able to continue. Please pass in an input image.')
        no_input_im_dialog.setStandardButtons(QMessageBox.Ok)
        no_input_im_dialog.exec()

    @staticmethod
    def no_guider_dialog():
        no_guider_dialog = QMessageBox()
        no_guider_dialog.setText('No guider chosen' + ' ' * 50)
        no_guider_dialog.setInformativeText('The tool will not be able to continue. Please select Guider 1 or 2.')
        no_guider_dialog.setStandardButtons(QMessageBox.Ok)
        no_guider_dialog.exec()

    def mismatched_apt_guider_dialog(self):
        mismatched_apt_guider_dialog = QMessageBox()
        mismatched_apt_guider_dialog.setText('Current guider does not match APT file' + ' ' * 30)
        mismatched_apt_guider_dialog.setInformativeText('The GUI is currently set to use GUIDER{}, but the APT file '
                                                        'you have selected has the guider set as GUIDER{}. '
                                                        'Please change your guider selection to match the '
                                                        'APT file.'.format(
                                                         int(self.buttonGroup_guider.checkedButton().text()),
                                                         self.apt_guider))
        mismatched_apt_guider_dialog.setStandardButtons(QMessageBox.Ok)
        mismatched_apt_guider_dialog.exec()

    def no_root_dialog(self):
        self.no_root_dialog_box = QMessageBox()
        self.no_root_dialog_box.setText('No root defined' + ' ' * 50)
        self.no_root_dialog_box.setInformativeText('The tool will not be able to continue. Please define a root.')
        self.no_root_dialog_box.setStandardButtons(QMessageBox.Ok)
        self.no_root_dialog_box.exec()

    def open_filename_dialog(self, title, multiple_files=False, file_type="All Files (*)"):
        """ Dialog box for opening a NIRCam of FGS image"""
        # options = QFileDialog.Options()
        if multiple_files==False:
            fileName, _ = QFileDialog.getOpenFileName(self, "Open {}".format(title),
                                                      "", file_type)
            # options=options)
            if fileName:
                return fileName
        else:
            fileName_list, _ = QFileDialog.getOpenFileNames(self, "Open {}".format(title),
                                                      "", file_type)
            # options=options)
            if fileName_list:
                return fileName_list

    def open_dirname_dialog(self):
        """ Dialog box for opening a directory"""
        options = QFileDialog.Options()
        dirName = QFileDialog.getExistingDirectoryUrl(self, "Select Directory").toString()[7:]
        if dirName:
            return dirName

    def segmentGuiding_dialog(self):
        # Initialize dialog widget
        SGT_dialog = QDialog()

        # Import .ui file
        uic.loadUi(os.path.join(__location__, 'segment_guiding', 'segmentGuidingDialog.ui'), SGT_dialog)

        # Set defaults from parsed header or commissioning section
        SGT_dialog.lineEdit_programNumber.setText(self.program_id)
        SGT_dialog.lineEdit_observationNumber.setText(self.observation_num)
        SGT_dialog.lineEdit_visitNumber.setText(self.visit_num)

        # Setting only for SOF, not POF
        try:
            SGT_dialog.lineEdit_RA.setText(self.gs_ra)
            SGT_dialog.lineEdit_DEC.setText(self.gs_dec)
        except AttributeError:
            pass

        # Run window and wait for response
        SGT_dialog.exec()

        return SGT_dialog

    def change_config_order_dialog(self, config_order):
        """Dialog box to change the order of the selected guiding configurations.
        Brought up only if self.checkBox_configorder is checked
        """
        self.change_config_order_dialog_box = QInputDialog()
        text, ok = self.change_config_order_dialog_box.getMultiLineText(self, 'Update Config Order',
                                                               'Re-order guiding selection files:' + ' ' * 100,
                                                               '\n'.join(config_order))
        if ok:
            config_list = [x.strip() for x in text.replace('\n', ",").split(',') if len(x.strip()) != 0]
            if len(config_list) != len(set(config_list)):
                raise ValueError('One or more guiding selections files were repeated in the "Change Config '
                                 'Order Dialog Box." Re-run the Segment Guiding section of this GUI and carefully '
                                 're-order the files.')
            else:
                return config_list
        else:
            raise ValueError('Change Config Order Dialog Box Canceled.')

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # HELPER FUNCTIONS
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    def load_input_image_data(self, filename):
        # Does the input image exist? If so, show it!
        if os.path.exists(filename) and filename.endswith('.fits'):
            # Switch to "Input Image" tab
            self.tabWidget.setCurrentIndex(0)

            # Load data and prep for showing on log scale
            data, header = utils.get_data_and_header(filename)
            data[data <= 0] = 1

            # Try to get information out of the header
            self.parse_header(filename)

            # Plot data
            self.canvas_input.compute_initial_figure(self.canvas_input.fig, data,
                                                     None, None, xlabel='X [pixels]',
                                                     ylabel='Y [pixels]')

            # Update filepath in textbox
            self.textEdit_showingInput.setEnabled(True)
            self.textEdit_showingInput.setText(filename)

        # Otherwise, show nothing
        else:
            # Update textbox showing filepath
            self.textEdit_showingInput.setText('No input image found at {}.'.format(filename))
            self.textEdit_showingInput.setEnabled(False)

            # Update plot to not show anything
            dummy_img = self.canvas_input.axes.imshow(
                np.array([[1e4, 1e4], [1e4, 1e4]]), cmap='bone', clim=(1e-1, 1e2)
            )

        self.canvas_input.draw()

    def parse_header(self, filename):
        # Open the header
        _, header = utils.get_data_and_header(filename)

        # Parse instrument
        try:
            instrument = header['INSTRUME'].upper()
            if "NIRCAM" in instrument:
                self.radioButton_NIRCam.setChecked(True)
            elif "FGS" in instrument:
                self.radioButton_FGS.setChecked(True)
        except KeyError:
            pass

        # Parse detector
        try:
            detector = header['DETECTOR'].upper()
            # Handle NIRCam case
            if "NRC" in detector:
                # Handle longwave case
                if "LONG" in detector:
                    module = detector[3]
                    if module == 'A':
                        index = 5
                    elif module == 'B':
                        index = 10

                # Handle shortwave case
                else:
                    detector = detector.strip()[-2:]
                    index = self.comboBox_detector.findText(detector, Qt.MatchFixedString)

                # Update dropdown menu
                if index >= 0:
                    self.comboBox_detector.setEnabled(True)
                    self.comboBox_detector.setCurrentIndex(index)

            # Handle FGS case
            elif "GUIDER" in detector:
                self.radioButton_FGS.setChecked(True)
        # If there is not DETECTOR keyword, set NIRCam detector back to parse
        except KeyError:
            self.comboBox_detector.setCurrentIndex(0)

        # Parse for ITM
        try:
            origin = header['ORIGIN'].strip()
            if origin == 'ITM':
                self.itm = True
                # self.checkBox_normalize.setEnabled(False) # Why is this set to FALSE? What did past self know that present self doesn't?
            else:
                self.itm = False
        except KeyError:
            pass

        # Parse APT program, observation, and visit information
        if self.radioButton_name_manual.isChecked() and not all(hasattr(self, attr) for attr in
                                                                ["program_id", "observation_num"]):
            keywords = ['PROGRAM', 'OBSERVTN', 'VISIT']
            attributes = ['program_id', 'observation_num', 'visit_num']
            for keyword, attr in zip(keywords, attributes):
                try:
                    setattr(self, attr, str(header[keyword]))
                except KeyError:
                    pass

    def is_valid_path_defined(self):
        if self.radioButton_name_manual.isChecked():
            if not (self.textEdit_out.toPlainText() != "" and self.lineEdit_root.text() != ""
                    and self.buttonGroup_guider.checkedButton()):
                return False
            else:
                return True
        elif self.radioButton_name_commissioning.isChecked():
            if self.textEdit_name_preview.toPlainText() == "" or not self.buttonGroup_guider.checkedButton():
                return False
            else:
                return True

    def update_guiding_selections(self, new_selections=None):
        """
        Update guiding selections in the converted preview combobox, the shifted preview combobox,
        and the guiding selections combobox in segment guiding section of the GUI

        Parameters
        ----------
        new_selections : list of str
            List of new unshifted guiding selection file paths that will be loaded into the GUI
            and need to show up as options in the above locations
        """
        # Update combobox in converted preview section if file(s) are loaded
        if isinstance(new_selections, list):
            self.guiding_selections_file_list = new_selections  # re-define to update chosen selections

            # Clear and reset 0th index in converted combo box
            self.comboBox_showcommandsconverted.blockSignals(True)
            self.comboBox_showcommandsconverted.clear()
            self.comboBox_showcommandsconverted.addItem('- Guiding Command -')
            self.comboBox_showcommandsconverted.blockSignals(False)

            # Remove circles on converted canvas
            try:
                for line in self.converted_im_circles:
                    line[0].set_visible(False)
                self.canvas_converted.peaks.set_visible(False)
                self.canvas_converted.draw()
            except AttributeError:
                pass

            # Clear and reset 0th index in shifted combo box
            self.comboBox_showcommandsshifted.blockSignals(True)
            self.comboBox_showcommandsshifted.clear()
            self.comboBox_showcommandsshifted.addItem('- Guiding Command -')
            self.comboBox_showcommandsshifted.blockSignals(False)

            # Remove circles on shifted canvas
            try:
                for line in self.shifted_im_circles:
                    line[0].set_visible(False)
                self.canvas_shifted.peaks.set_visible(False)
                self.canvas_shifted.draw()
            except AttributeError:
                pass

        # Add chosen files to converted image combobox to choose from
        if len(self.comboBox_showcommandsconverted) == 1 and "Guiding Command" in \
                self.comboBox_showcommandsconverted.currentText():
            for i, command_file in enumerate(self.guiding_selections_file_list):
                item = "Command {}: {}".format(i + 1, command_file.split('/')[-1])
                self.comboBox_showcommandsconverted.addItem(item)

        # Add chosen files to shifted image combobox to choose from
        if len(self.comboBox_showcommandsshifted) == 1 and "Guiding Command" in \
                self.comboBox_showcommandsshifted.currentText():
            self.comboBox_showcommandsshifted.clear()
            for i, command_file in enumerate(self.shifted_guiding_selections_file_list):
                item = "Command {}: {}".format(i + 1, command_file.split('/')[-1])
                self.comboBox_showcommandsshifted.addItem(item)

        self.update_checkable_combobox()

    def update_checkable_combobox(self):
        if self.groupBox_segmentGuiding.isChecked():
            # Define path to look for files based on contents of lineEdit_regfileSegmentGuiding
            if self.sender() == self.lineEdit_regfileSegmentGuiding:
                path = self.lineEdit_regfileSegmentGuiding.text()
                root = path.split('/out/')[-1].split('/')[0]
                if root == '':
                    root = '*'
            else:
                if self.radioButton_name_manual.isChecked():
                    root_dir = os.path.join(self.textEdit_out.toPlainText(), 'out', self.lineEdit_root.text())
                    path = root_dir
                    root = self.lineEdit_root.text()
                else:
                    root_dir = self.textEdit_name_preview.toPlainText()
                    path = root_dir
                    root = root_dir.split('out/')[-1].split('/')[0]

                self.lineEdit_regfileSegmentGuiding.setText(root_dir)

            # Re-define files based on path found above
            txt_files = glob.glob(os.path.join(path, "**/*.txt"), recursive=True)
            acceptable_guiding_files_list, acceptable_all_psf_files_list = \
                self.search_acceptable_files(path, root, '*', shifted=self.radioButton_shifted.isChecked())

            # This would occur in the case of a POF
            if len(txt_files) == 0:
                LOGGER.warning('Master GUI: No guiding_selections and/or all_found_psf files found. This may be okay '
                               'depending on the situation (e.g. when making a POF).')

                # Clear the guiding selections combo box
                self.comboBox_guidingcommands.clear()

                # Clear and reset 0th index in shifted combo box
                self.comboBox_showcommandsshifted.blockSignals(True)
                self.comboBox_showcommandsshifted.clear()
                self.comboBox_showcommandsshifted.addItem('- Guiding Command -')
                self.comboBox_showcommandsshifted.blockSignals(False)

                pass

            else:
                try:
                    if not self.radioButton_shifted.isChecked():
                        self.guiding_selections_file_list = sorted([file for f in acceptable_guiding_files_list for
                                                                    file in fnmatch.filter(txt_files, f)],
                                                                   key=utils.natural_keys)
                        self.all_found_psfs_file = sorted([file for f in acceptable_all_psf_files_list for file in
                                                           fnmatch.filter(txt_files, f)], key=utils.natural_keys)[0]
                        guiding_selections = self.guiding_selections_file_list
                    else:
                        self.shifted_guiding_selections_file_list = sorted(fnmatch.filter(
                            txt_files, acceptable_guiding_files_list[0]), key=utils.natural_keys)
                        self.shifted_all_found_psfs_file_list = sorted(fnmatch.filter(
                            txt_files, acceptable_all_psf_files_list[0]), key=utils.natural_keys)
                        guiding_selections = self.shifted_guiding_selections_file_list

                    # Clear and re-populate checkable combobox in segment guiding section
                    self.comboBox_guidingcommands.clear()
                    i = 0
                    for command_file in guiding_selections:
                        # If the wrong guider number is present in the name, don't list the file
                        if '_G{}_'.format([2 if self.buttonGroup_guider.checkedButton().text() == '1' else 1][0]) not in command_file:
                            item = "Command {}: {}".format(i + 1, command_file.split('/')[-1])
                            self.comboBox_guidingcommands.addItem(item)
                            item = self.comboBox_guidingcommands.model().item(i, 0)
                            item.setFlags(QtCore.Qt.ItemIsUserCheckable | QtCore.Qt.ItemIsEnabled)
                            item.setCheckState(Qt.Unchecked)
                            i += 1

                except IndexError:
                    LOGGER.warning(
                        'Master GUI: Missing guiding_selections and/or all_found_psf files. This may be okay depending '
                        'on the situation (e.g. when making a POF).')

                    # Clear the guiding selections combo box
                    self.comboBox_guidingcommands.clear()

                    # Clear and reset 0th index in shifted combo box
                    self.comboBox_showcommandsshifted.blockSignals(True)
                    self.comboBox_showcommandsshifted.clear()
                    self.comboBox_showcommandsshifted.addItem('- Guiding Command -')
                    self.comboBox_showcommandsshifted.blockSignals(False)

    def update_regfile_starselector_combobox(self, files):
        if len(files) != 0:
            # Clear and populate
            self.comboBox_regfileStarSelector.clear()
            for item in files:
                self.comboBox_regfileStarSelector.addItem(item)

            self.comboBox_regfileStarSelector.lineEdit().setReadOnly(True)

            # Set width to match the longest entry name
            w = self.comboBox_regfileStarSelector.fontMetrics().boundingRect(max(files, key=len)).width()
            self.comboBox_regfileStarSelector.view().setFixedWidth(w + 10)

    def update_converted_image_preview(self):
        # Are all the necessary fields filled in? If not, don't even try.
        if not self.is_valid_path_defined():
            return

        # Does a converted image exist? If so, show it!
        if os.path.exists(self.converted_im_file):
            # Prepare to show converted image
            self.canvas_converted.axes.set_visible(True)
            self.tabWidget.setCurrentIndex(1)
            self.textEdit_showingConverted.setEnabled(True)

            # Update filepath
            self.textEdit_showingConverted.setText(self.converted_im_file)

            # Toggle the "use converted image" buttons
            self.checkBox_useConvertedImage.setEnabled(True)
            self.checkBox_useConvertedImage.setChecked(True)

            # Enable the "show stars" button
            self.checkBox_showStars.setEnabled(True)

            # Enable guiding commands button
            self.comboBox_showcommandsconverted.setEnabled(True)

            # Load data
            data, hdr = utils.get_data_and_header(self.converted_im_file)
            data[data <= 0] = 1

            # Load all_found_psfs*.text
            x, y = [None, None]
            if os.path.exists(self.all_found_psfs_file):
                psf_list = asc.read(self.all_found_psfs_file)
                x = psf_list['x']
                y = psf_list['y']

            # Plot data image and peak locations from all_found_psfs*.txt
            self.canvas_converted.compute_initial_figure(self.canvas_converted.fig, data, x, y)

            # Check for distortion keyword and add banner
            if 'DISTORT' not in hdr.keys():
                self.canvas_converted.axes.text(400., 100.,
                                                'File loaded is incorrectly distorted',
                                                fontsize=11,
                                                bbox={'facecolor': 'lightcoral', 'pad': 7})
            else:
                for txt in self.canvas_converted.axes.texts:
                    txt.set_visible(False)

            # If possible, plot the selected stars in the guiding_selections*.txt
            if self.comboBox_showcommandsconverted.currentIndex() != 0:
                guiding_selections_file = self.guiding_selections_file_list[
                    self.comboBox_showcommandsconverted.currentIndex()-1]
                if os.path.exists(guiding_selections_file):
                    selected_psf_list = asc.read(guiding_selections_file)
                    x_selected = selected_psf_list['x']
                    y_selected = selected_psf_list['y']

                    # Remove old circles
                    for line in self.converted_im_circles:
                        self.canvas_converted.axes.lines.remove(line[0])
                    self.converted_im_circles = [
                        self.canvas_converted.axes.plot(x_selected[0], y_selected[0],
                                                        'o', ms=25, mfc='none',
                                                        mec='yellow', mew=2, lw=0)
                    ]
                    self.converted_im_circles.append(
                        self.canvas_converted.axes.plot(x_selected[1:], y_selected[1:],
                                                        'o', ms=25, mfc='none',
                                                        mec='darkorange', mew=2, lw=0)
                    )
            else:
                for line in self.converted_im_circles:
                    line[0].set_visible(False)
                self.canvas_converted.peaks.set_visible(False)

                self.canvas_converted.draw()

        # If not, show nothing.
        else:
            # Update textbox showing filepath
            self.textEdit_showingConverted.setText(
                'No converted guider {} image found at {}.'.format(self.buttonGroup_guider.checkedButton().text(),
                                                                   self.converted_im_file))
            self.textEdit_showingConverted.setEnabled(False)

            # Clear and reset 0th index in shifted combo box
            self.comboBox_showcommandsconverted.blockSignals(True)
            self.comboBox_showcommandsconverted.clear()
            self.comboBox_showcommandsconverted.addItem('- Guiding Command -')
            self.comboBox_showcommandsconverted.blockSignals(False)

            # Disable the "use converted image" buttons
            self.checkBox_useConvertedImage.setChecked(False)
            self.checkBox_useConvertedImage.setEnabled(False)

            # Disable the "show stars" button
            self.checkBox_showStars.setEnabled(False)
            self.comboBox_showcommandsconverted.setEnabled(False)

            # Update plot to not show anything
            dummy_img = self.canvas_converted.axes.imshow(
                np.array([[1e4, 1e4], [1e4, 1e4]]), cmap='bone', clim=(1e-1, 1e2)
            )

        return self.canvas_converted.draw()

    def update_shifted_image_preview(self):
        # Are all the necessary fields filled in? If not, don't even try.
        if not self.is_valid_path_defined():
            return

        # Is the self.shifted_im_file_list variable defined yet?
        noshow = False
        try:
            if len(self.shifted_im_file_list) == 0:
                noshow = True
        except AttributeError:
            noshow = True

        # Do all the shifted images exist? If so, show them!
        if noshow is False and False not in [os.path.exists(shifted_im_file) for
                                             shifted_im_file in self.shifted_im_file_list]:
            # Prepare to show shifted image
            self.canvas_shifted.axes.set_visible(True)
            self.tabWidget.setCurrentIndex(2)
            self.textEdit_showingShifted.setEnabled(True)

            # Toggle the "use shifted image" buttons
            self.radioButton_shifted.setChecked(True)

            # Enable the "show stars" button
            self.checkBox_showStars_shifted.setEnabled(True)

            # Enable and populate guiding commands button
            self.comboBox_showcommandsshifted.setEnabled(True)
            i = self.comboBox_showcommandsshifted.currentIndex()

            # Update filepath
            self.textEdit_showingShifted.setText(self.shifted_im_file_list[i])

            # Load data
            data, hdr = utils.get_data_and_header(self.shifted_im_file_list[i])
            data[data <= 0] = 1

            # Load all_found_psfs*.text
            x, y = [None, None]
            if os.path.exists(self.shifted_all_found_psfs_file_list[i]):
                psf_list = asc.read(self.shifted_all_found_psfs_file_list[i])
                x, y = np.array([(xi, yi) for (xi, yi) in zip(psf_list['x'], psf_list['y'])
                                 if (xi < 2048) & (yi < 2028)]).T

            # Plot data image and peak locations from all_found_psfs*.txt
            self.canvas_shifted.compute_initial_figure(self.canvas_shifted.fig, data, x, y)

            # Check for distortion keyword and add banner
            if 'DISTORT' not in hdr.keys():
                self.canvas_shifted.axes.text(400., 100.,
                                              'File loaded is incorrectly distorted',
                                              fontsize=11,
                                              bbox={'facecolor': 'lightcoral', 'pad': 7})
            else:
                for txt in self.canvas_shifted.axes.texts:
                    txt.set_visible(False)

            # If possible, plot the selected stars in the guiding_selections*.txt
            shifted_guiding_selections_file = self.shifted_guiding_selections_file_list[i]

            if os.path.exists(shifted_guiding_selections_file):
                selected_psf_list = asc.read(shifted_guiding_selections_file)
                x_selected, y_selected = np.array([(xi, yi) for (xi, yi) in
                                                   zip(selected_psf_list['x'], selected_psf_list['y'])
                                                   if (xi < 2048) & (yi < 2028)]).T

                # Remove old circles
                for line in self.shifted_im_circles:
                    self.canvas_shifted.axes.lines.remove(line[0])

                self.shifted_im_circles = [
                    self.canvas_shifted.axes.plot(x_selected[0], y_selected[0],
                                                  'o', ms=25, mfc='none',
                                                  mec='yellow', mew=2, lw=0)
                ]
                self.shifted_im_circles.append(
                    self.canvas_shifted.axes.plot(x_selected[1:], y_selected[1:],
                                                  'o', ms=25, mfc='none',
                                                  mec='darkorange', mew=2, lw=0)
               )
            else:
                for line in self.shifted_im_circles:
                    line[0].set_visible(False)
                try:
                    self.canvas_shifted.peaks.set_visible(False) # won't be defined until after first pass
                except AttributeError:
                    pass
                self.canvas_shifted.draw()

        # If not, show nothing.
        else:
            # Update textbox showing filepath
            if len(self.shifted_im_file_list) == 0:
                self.textEdit_showingShifted.setText(
                    'No shifted guider {} image found at {}.'.format(self.buttonGroup_guider.checkedButton().text(),
                        '/'.join(self.converted_im_file.split('/')[:-2] + \
                                 ['guiding_config_*/FGS_imgs/shifted_*.fits'])))

            self.textEdit_showingShifted.setEnabled(False)

            # Clear and reset 0th index in shifted combo box
            self.comboBox_showcommandsshifted.blockSignals(True)
            self.comboBox_showcommandsshifted.clear()
            self.comboBox_showcommandsshifted.addItem('- Guiding Command -')
            self.comboBox_showcommandsshifted.blockSignals(False)

            # Uncheck the "use shifted image" button
            self.radioButton_unshifted.setChecked(True)

            # Disable the "show stars" button
            self.checkBox_showStars_shifted.setEnabled(False)
            self.comboBox_showcommandsshifted.setEnabled(False)

            # Update plot to not show anything
            dummy_img = self.canvas_shifted.axes.imshow(
                np.array([[1e4, 1e4], [1e4, 1e4]]), cmap='bone', clim=(1e-1, 1e2)
            )
        return self.canvas_shifted.draw()

    @staticmethod
    def search_acceptable_files(root_dir, root, guider, shifted):
        if shifted is False:
            acceptable_guiding_files_list = [
                os.path.join(root_dir, 'guiding_config_*',
                             'unshifted_guiding_selections_{}_G{}_config*.txt'.format(root, guider)),  # newest
                os.path.join(root_dir, 'guiding_config_*', 'guiding_selections_{}_G{}.txt'.format(root, guider)),
                os.path.join(root_dir, 'guiding_config_*', '{}_G{}_regfile.txt'.format(root, guider)),
                os.path.join(root_dir, 'unshifted_guiding_selections_{}_G{}.txt'.format(root, guider)),
                os.path.join(root_dir, 'guiding_selections_{}_G{}.txt'.format(root, guider)),
                os.path.join(root_dir, '{}_G{}_regfile.txt'.format(root, guider))]  # oldest

            acceptable_all_psf_files_list = [
                os.path.join(root_dir, 'unshifted_all_found_psfs_{}_G{}.txt'.format(root, guider)),
                os.path.join(root_dir, 'all_found_psfs_{}_G{}.txt'.format(root, guider)),
                os.path.join(root_dir, '{}_G{}_ALLpsfs.txt'.format(root, guider))]

        else:
            acceptable_guiding_files_list = [
                os.path.join(root_dir, 'guiding_config_*', 'shifted_guiding_selections_{}_G{}_config*.txt'.format(
                    root, guider)),
                os.path.join(root_dir,'shifted_guiding_selections_{}_G{}_config*.txt'.format(root, guider)),
                os.path.join(root_dir,'shifted_guiding_selections_{}_G{}.txt'.format(root, guider)),
            ]

            acceptable_all_psf_files_list = [
                os.path.join(root_dir, 'guiding_config_*', 'shifted_all_found_psfs_{}_G{}_config*.txt'.format(
                    root, guider)),
                os.path.join(root_dir, 'shifted_all_found_psfs_{}_G{}_config*.txt'.format(root, guider)),
                os.path.join(root_dir, 'shifted_all_found_psfs_{}_G{}.txt'.format(root, guider))
            ]

        return acceptable_guiding_files_list, acceptable_all_psf_files_list

    def update_filepreview(self, new_guiding_selections=False):
        """
        If either: 
          1) manual naming is selected and the root, out_dir, and guider have been defined, or
          2) commissioning naming is selected and the practice, CAR, and observation have
              been selected
        show an example filepath to a simulated image, and auto-populate the guiding_selections*.txt.text
        filepath. Also, auto-generate important filenames as class attributes.
        """
        # If you change the main path, update the guiding selections based on that information
        if self.sender() in [self.pushButton_inputImage, self.lineEdit_inputImage, self.buttonGroup_guider,
                             self.lineEdit_root, self.textEdit_out]:
            new_guiding_selections = True

        if self.is_valid_path_defined():
            manual_naming = self.radioButton_name_manual.isChecked()
            guider = self.buttonGroup_guider.checkedButton().text()

            # Determine root directory
            if manual_naming:
                root = self.lineEdit_root.text()
                root_dir = os.path.join(self.textEdit_out.toPlainText(), 'out',
                                        self.lineEdit_root.text())
            else:
                root = 'for_obs{:02d}'.format(int(self.lineEdit_obs.text()))
                root_dir = self.textEdit_name_preview.toPlainText()

            # Create root directory if it doesn't exist
            utils.ensure_dir_exists(root_dir)

            # Set log if not already set (for first file created with MAGIC)
            if self.log is None:
                self.log, self.log_filename = utils.create_logger_from_yaml('magic', out_dir_root=root_dir,
                                                                            root=root, level='DEBUG')
            # If path has changed, need to create a new log file
            if root_dir != os.path.dirname(self.log_filename):
                self.log, self.log_filename = utils.create_logger_from_yaml('magic', out_dir_root=root_dir,
                                                                            root=root, level='DEBUG')

            # Note: maintaining if statements and "old" file names for backwards compatibility
            txt_files = glob.glob(os.path.join(root_dir, "**/*.txt"), recursive=True)
            fits_files = glob.glob(os.path.join(root_dir, "**/*.fits"), recursive=True)
            acceptable_guiding_files_list, acceptable_all_psf_files_list = \
                self.search_acceptable_files(root_dir, root, guider, shifted=False)

            # Pull every possible guiding selections file and 1 all found psfs file
            try:
                self.guiding_selections_file_list = sorted([file for f in acceptable_guiding_files_list
                                                     for file in fnmatch.filter(txt_files, f)], key=utils.natural_keys)
                self.all_found_psfs_file = sorted([f for f in acceptable_all_psf_files_list if f in txt_files],
                                                  key=utils.natural_keys)[0]
            except IndexError:
                self.guiding_selections_file_list = []
                self.all_found_psfs_file = ''

            self.guiding_selections_file_list_default = self.guiding_selections_file_list  # to go back to default found

            # Update converted FGS image filepath
            self.converted_im_file = os.path.join(root_dir, 'FGS_imgs', 'unshifted_{}_G{}.fits'.format(root, guider))

            # Update shifted FGS image & catalog filepaths
            shifted_acceptable_guiding_files_list, shifted_acceptable_all_psf_files_list = \
                self.search_acceptable_files(root_dir, root, guider, shifted=True)

            self.shifted_all_found_psfs_file_list = sorted(fnmatch.filter(
                txt_files, shifted_acceptable_all_psf_files_list[0]), key=utils.natural_keys)

            self.shifted_guiding_selections_file_list = sorted(fnmatch.filter(
                txt_files, shifted_acceptable_guiding_files_list[0]), key=utils.natural_keys)

            # Update shifted FGS image(s) filepath
            self.shifted_im_file_list = sorted(fnmatch.filter(
                fits_files, os.path.join(root_dir, 'guiding_config_*', 'FGS_imgs',
                                         'shifted_{}_G{}_config*.fits'.format(root, guider))), key=utils.natural_keys)

            # Update guiding_selections*.txt paths in GUI
            if False not in [os.path.exists(file) for file in self.guiding_selections_file_list] and len(self.guiding_selections_file_list) != 0:
                if new_guiding_selections:
                    new_selections = self.guiding_selections_file_list
                    self.update_guiding_selections(new_selections=new_selections)
            else:
                self.comboBox_regfileStarSelector.clear()
                self.lineEdit_regfileSegmentGuiding.setText(root_dir)
                self.comboBox_guidingcommands.clear()

            # If possible, show converted and shifted image previews, too
            self.update_converted_image_preview()
            self.update_shifted_image_preview()

    def query_apt_for_gs(self, program_id, obs_number):
        """
        Function to take in the desired APT Program ID and observation number and
        output the ID, RA, and DEC of the guide star set under the Special Requirements
        tab of that APT file for that ID/Obs. The RA and Dec is found by querying the GSC
        Catalog using the jwst-fgs-countrate module.

        Code adapted from:
            https://github.com/spacetelescope/jwst_magic/blob/master/fgs-commissioning/
                notebooks/generate_commissioning_activities_yaml.ipynb
            https://grit.stsci.edu/jsahlmann/aptxml/blob/master/aptxml/manipulate.py

        Parameters
        ----------
        program_id : int or str
            The APT program ID of interest
        obs_number : int or str
            The observation number of interest

        Returns
        -------
        gs_id : str
            ID of the guide star from the program ID/Obs
        guider : int
            The guider number from the program ID/Obs
        ra : float
            RA of the guide star from the program ID/Obs
        dec : float
            DEC of the guide star from the program ID/Obs

        """
        # Check program_id and obs_number are ints
        if not isinstance(program_id, int):
            program_id = int(program_id)
        if not isinstance(obs_number, int):
            obs_number = int(obs_number)

        LOGGER.info("Master GUI: Checking Program ID {} and Obs #{}".format(program_id, obs_number))

        # Build temporary directory
        apt_file_path = os.path.join(__location__, 'data', 'temp_apt')
        if os.path.exists(apt_file_path):
            shutil.rmtree(apt_file_path)
        os.mkdir(apt_file_path)

        # Download the APT file (.aptx file) from online using the program ID
        urllib.request.urlretrieve(
            'http://www.stsci.edu/jwst/phase2-public/{}.aptx'.format(program_id),
            '{}/{}.aptx'.format(apt_file_path,program_id)
        )

        # Pull XML data from APT file
        with zipfile.ZipFile('{}/{}.aptx'.format(apt_file_path,program_id), 'r') as zf:
            data = zf.read('{}.xml'.format(program_id))
        tree = etree.parse(io.BytesIO(data))

        # Pull: List of Observations for the CAR > Specific Obs > Special Requirements Info (sr)
        namespace_tag = '{http://www.stsci.edu/JWST/APT}'
        observation_list = tree.find(namespace_tag + 'DataRequests').findall('.//' + namespace_tag + 'Observation')

        # Find the right observation number and pull that observation
        for i, obs in enumerate(observation_list):
            num = [x for x in obs.iterchildren() if x.tag.split(namespace_tag)[1] == "Number"][0].text
            if int(num) == obs_number:
                break
            if i == len(observation_list) - 1:  # if you get to the end of the list and don't find a matching obs number
                shutil.rmtree(apt_file_path)
                raise ValueError("This program doesn't have an observation {}".format(obs_number))
        observation = observation_list[i]
        if int([x for x in observation.iterchildren() if x.tag.split(namespace_tag)[1] == "Number"][0].text) \
                != obs_number:
            raise ValueError('Failed to find the right observation. The user will need to check the APT file by hand.')

        sr = [x for x in observation.iterchildren() if x.tag.split(namespace_tag)[1] == "SpecialRequirements"][0]

        # Try to pull the Guide Star information
        try:
            gs = [x for x in sr.iterchildren() if x.tag.split(namespace_tag)[1] == "GuideStarID"][0]
        except IndexError:
            self.lineEdit_normalize.setText('')
            shutil.rmtree(apt_file_path)
            raise ValueError("This observation doesn't have a Guide Star Special Requirement")

        # Pull out the guide star ID and the guider number
        gs_id = [x for x in gs.iterchildren() if x.tag.split(namespace_tag)[1] == "GuideStar"][0].text
        guider = [x for x in gs.iterchildren() if x.tag.split(namespace_tag)[1] == "Guider"][0].text

        # Account for if the guider is written as "guider 1" or "guider1"
        if not isinstance(guider, int):
            guider = guider.lower().replace(' ', '').split('guider')[1]

        LOGGER.info('Master GUI: APT has been queried and found guide star {} and guider {}'.format(gs_id, guider))

        # Tear down temporary directory
        shutil.rmtree(apt_file_path)

        # Use Guide Star ID to get RA/DEC using default GSC in fgscountrate module
        try:
            ra, dec = renormalize.query_guide_star_catalog(gs_id=gs_id)
        except ValueError as err:
            self.lineEdit_normalize.setText('')
            raise ValueError(str(err))

        LOGGER.info('Master GUI: The Guide Star Catalog have been queried and found RA of {} and DEC of {} '.format(ra,
                                                                                                                dec))

        return gs_id, guider, ra, dec


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# MAIN FUNCTION
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

def run_MasterGui(root=None, norm_value=12.0, norm_unit="FGS Magnitude", nircam_det=None,
                  nircam=True, smoothing='default', steps=None, in_file=None,
                  bkgd_stars=False, out_dir='', convert_im=True,
                  star_selection=True, star_selection_gui=True, file_writer=True,
                  segment_guiding=False, itm=False):
    # RUN GUI
    app = QtCore.QCoreApplication.instance()  # Use existing instance, if there is one
    if app is None:
        app = QApplication(sys.argv)

    ex = MasterGui(root, norm_value, norm_unit, nircam_det, nircam, smoothing,
                   steps, in_file, bkgd_stars, out_dir, convert_im, star_selection_gui,
                   file_writer, segment_guiding, app=app, itm=itm)
    # #return ex.settings
    out = app.exec()

    return out
