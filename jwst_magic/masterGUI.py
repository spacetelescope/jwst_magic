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
                             QDialog)
import yaml

# Local Imports
from jwst_magic import run_magic
from jwst_magic.convert_image import renormalize, background_stars
from jwst_magic.fsw_file_writer import rewrite_prc
from jwst_magic.segment_guiding import segment_guiding
from jwst_magic.star_selector.SelectStarsGUI import StarClickerMatplotlibCanvas, run_SelectStars
from jwst_magic.utils import utils

# Define all needed paths
__location__ = os.path.realpath(os.path.join(os.getcwd(), os.path.dirname(__file__)))
PACKAGE_PATH = os.path.dirname(os.path.realpath(__file__))
OUT_PATH = os.path.split(PACKAGE_PATH)[0]  # Location of out/ and logs/ directory
SOGS_PATH = '/data/jwst/wss/guiding/'

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

        # Move the cursor to the end so that the most recent  text is visible
        cursor = self.textEdit_log.textCursor()
        cursor.movePosition(QtGui.QTextCursor.End)
        self.textEdit_log.setTextCursor(cursor)


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# GUI CLASS DEFINITION
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


class MasterGui(QMainWindow):
    def __init__(self, root=None, norm_value=12.0, norm_units='FGS Magnitude',
                 nircam_det=None, nircam=True, smoothing='default',
                 steps=None, in_file=None, bkgd_stars=False, out_dir=OUT_PATH, convert_im=True,
                 star_selection=True, star_selection_gui=True, file_writer=True,
                 segment_guiding=False, app=None, itm=False):

        # Initialize attributes
        self.app = app
        self.commissioning_dict = {}
        self.root_default = root
        self.converted_im_circles = []
        self.shifted_im_circles = []
        self.bkgd_stars = None
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
        """Set up the two matplotlib canvases that will preview the
        input image and converted image in the "Image Preview" section.
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
        self.pushButton_out.clicked.connect(self.on_click_out)
        self.textEdit_out.installEventFilter(self)
        self.pushButton_manualid.clicked.connect(self.update_apt_gs_values)

        # Image convertor widgets
        self.pushButton_backgroundStars.clicked.connect(self.on_click_bkgdstars)
        self.horizontalSlider_coarsePointing.sliderReleased.connect(self.on_change_jitter)
        self.lineEdit_coarsePointing.editingFinished.connect(self.on_change_jitter)

        # Star selector widgets
        self.pushButton_regfileStarSelector.clicked.connect(self.on_click_infile)

        # Segment guiding widgets
        self.pushButton_regfileSegmentGuiding.clicked.connect(self.on_click_infile)
        self.buttonGroup_segmentGuiding_idAttitude.buttonClicked.connect(self.update_segment_guiding_shift)
        self.groupBox_segmentGuiding.toggled.connect(self.update_segment_guiding_shift)

        # Image preview widgets
        self.checkBox_showStars.toggled.connect(self.on_click_showstars)
        self.checkBox_showStars_shifted.toggled.connect(self.on_click_showstars)

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
        self.comboBox_obs.currentIndexChanged.connect(self.update_commissioning_name)
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

                self.update_filepreview()

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
            root = 'for_obs{:02d}'.format(int(self.comboBox_obs.currentText()))
            out_dir = os.path.join(SOGS_PATH,
                                   self.comboBox_practice.currentText(),
                                   self.comboBox_car.currentText().lower().replace('-', ''),
                                   )
        copy_original = True

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
        itm = self.itm

        # Handle the case where we want to use a pre-existing converted image
        pre_existing_im = self.checkBox_useConvertedImage.isChecked() and \
                          convert_im and \
                          self.checkBox_useConvertedImage.isEnabled()
        if pre_existing_im:
            convert_im = False
            input_image = self.converted_im_file
            copy_original = False

        # Star selection
        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        if self.checkBox_globalAlignment.isChecked():
            smoothing = 'high'
        elif self.checkBox_noSmoothing.isChecked():
            smoothing = 'low'
        else:
            smoothing = 'default'

        star_selection = self.groupBox_starSelector.isChecked()
        star_selectiongui = self.radioButton_starSelectorGUI.isChecked()
        if not self.radioButton_regfileStarSelector.isChecked():
            in_file = None
        else:
            in_file = self.lineEdit_regfileStarSelector.text()

        # File writer
        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        file_writer = self.groupBox_fileWriter.isChecked()
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

        # Shift image to ID attitude:
        shift_id_attitude = self.checkBox_id_attitude.isChecked()
        crowded_field = self.radioButton_crowded_id_attitude.isChecked()

        # Rewrite .prc and guiding_selections*.txt ONLY
        if self.checkBox_rewritePRC.isChecked():
            # Open converted FGS file
            data, _ = utils.get_data_and_header(self.converted_im_file)

            # Update array for showing in LogNorm
            data[data <= 0] = 1

            # Open all_found_psfs*.txt (list of all identified segments)
            all_psfs = self.all_found_psfs_file
            all_rows = asc.read(all_psfs)
            x = all_rows['x'].data
            y = all_rows['y'].data

            # Run the select stars GUI to determine the new orientation
            inds = run_SelectStars(data, x, y, 20, masterGUIapp=self.app)

            # Rewrite the id.prc and acq.prc files
            rewrite_prc.rewrite_prc(inds, guider, root, out_dir)

            # Update converted image preview
            self.update_filepreview()
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
                # Initialize the dialog
                self._test_sg_dialog = segment_guiding.SegmentGuidingGUI.SegmentGuidingDialog(
                                       "POF", None, self.program_id, self.observation_num, self.visit_num, log=None
                )
                # Generate the file
                segment_guiding.generate_photometry_override_file(
                    root, self.program_id, self.observation_num, self.visit_num, out_dir=out_dir,
                    parameter_dialog=True, dialog_obj=self._test_sg_dialog, log=LOGGER
                )

            else:
                # Get APT program information from parsed header
                self.parse_header(input_image)

                # Define location of all_found_psfs catalog file
                if self.lineEdit_regfileSegmentGuiding.text() is '':
                    if self.radioButton_shifted.isChecked():
                        segment_infile = self.shifted_all_found_psfs_file
                    else:
                        segment_infile = self.all_found_psfs_file
                else:
                    segment_infile = self.lineEdit_regfileSegmentGuiding.text()

                # Verify that the all_found_psfs*.txt file exists
                if not os.path.exists(segment_infile):
                    raise OSError('Provided segment infile {} not found.'.format(segment_infile))

                # Determine which image to use for the click-to-select GUI in generate_segment_override_file and load it
                if self.radioButton_shifted.isChecked():
                    fgs_filename = self.shifted_im_file
                elif self.radioButton_unshifted.isChecked() and hasattr(self, 'converted_im_file') and \
                        os.path.exists(self.converted_im_file):
                    fgs_filename = self.converted_im_file
                elif self.radioButton_unshifted.isChecked():
                    fgs_filename = input_image
                data, _ = utils.get_data_and_header(fgs_filename)

                # Determine whether to load guiding_selections*.txt or run GUI
                GUI = not self.radioButton_regfileSegmentGuiding.isChecked()
                selected_segs = self.lineEdit_regfileSegmentGuiding.text()

                # Run the tool and generate the file
                # Initialize the dialog
                self._test_sg_dialog = segment_guiding.SegmentGuidingGUI.SegmentGuidingDialog(
                    "SOF", guider, self.program_id, self.observation_num, self.visit_num,
                    ra=self.gs_ra, dec=self.gs_dec,log=None
                )

                # Generate the file
                segment_guiding.generate_segment_override_file(
                    segment_infile, guider, self.program_id, self.observation_num,
                    self.visit_num, ra=self.gs_ra, dec=self.gs_dec,
                    root=root, out_dir=out_dir, selected_segs=selected_segs,
                    click_to_select_gui=GUI, data=data, master_gui_app=self.app,
                    parameter_dialog=True, dialog_obj=self._test_sg_dialog, log=LOGGER
                )

            # Update converted image preview
            self.update_filepreview()
            return
        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

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
                              out_dir=out_dir,
                              convert_im=convert_im,
                              star_selection=star_selection,
                              star_selection_gui=star_selectiongui,
                              file_writer=file_writer,
                              masterGUIapp=self.app,
                              copy_original=copy_original,
                              normalize=normalize,
                              coarse_pointing=coarse_point,
                              jitter_rate_arcsec=jitter_rate_arcsec,
                              itm=itm,
                              shift_id_attitude=shift_id_attitude,
                              crowded_field=crowded_field
                              )

            # Update converted image preview
            self.update_filepreview()

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

        # Derive the root from the filename and assume the default output
        # directory (OUT_PATH)
        if self.radioButton_name_manual.isChecked():
            if self.textEdit_out.toPlainText() == "":
                self.textEdit_out.setEnabled(True)
                self.textEdit_out.setText(OUT_PATH)
            if self.textEdit_out.toPlainText() == OUT_PATH:
                self.textEdit_out.setEnabled(True)

        # Update the example filepath (and converted image preview, if possible)
        self.update_filepreview()

        # Show input image preview
        self.load_input_image_data(filename)

        # Show converted image preview, if possible
        self.update_converted_image_preview()

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

            self.update_filepreview()

    def on_click_infile(self):
        """ Using the Infile Open button (open file) """
        # Determine which infile is being edited
        if self.sender() == self.pushButton_regfileStarSelector:
            to_text = self.lineEdit_regfileStarSelector
        elif self.sender() == self.pushButton_regfileSegmentGuiding:
            to_text = self.lineEdit_regfileSegmentGuiding

        filename = self.open_filename_dialog('In/Reg file', file_type="Input file (*.txt *.incat);;All files (*.*)")
        to_text.setText(filename)
        return filename

    def on_click_root(self):
        """ Using the Set Default Root button"""
        root = utils.make_root(None, self.lineEdit_inputImage.text())
        self.lineEdit_root.setText(root)
        self.update_filepreview()
        return root

    def on_click_bkgdstars(self):
        if self.lineEdit_inputImage.text() == "":
            self.no_inputImage_dialog()
            return

        if not self.buttonGroup_guider.checkedButton():
            self.no_guider_dialog()
            return

        # Enable the textbox
        self.textEdit_backgroundStars.setEnabled(True)

        guider = int(self.buttonGroup_guider.checkedButton().text())

        # If normalization is turned on, read the normalization value & unit
        # and calculate JMag of the guidestar
        if self.checkBox_normalize.isChecked():
            norm_value = float(self.lineEdit_normalize.text())
            norm_unit = self.comboBox_normalize.currentText()
            #TODO - this will have to be changed
            norm_obj = renormalize.NormalizeToCountrate(norm_value, norm_unit, guider)
            fgs_countrate = norm_obj.to_countrate()
            jmag = renormalize.fgs_countrate_to_j_mag(fgs_countrate, guider)
        # If not, determine the FGS counts of the input image
        else:
            input_image = self.lineEdit_inputImage.text()
            data, _ = utils.get_data_and_header(input_image)
            fgs_countrate = np.sum(data[data > np.median(data)])
            jmag = renormalize.fgs_countrate_to_j_mag(fgs_countrate, guider) #We still need this??

        # Eventually, we will want to this to be done with FGS mag
        self.bkgd_stars, method = background_stars.run_background_stars_GUI(guider, jmag, masterGUIapp=self.app)

        method_adverb = {'random': 'randomly',
                         'user-defined': 'as defined by the user',
                         'catalog': 'from a GSC query'}

        if isinstance(self.bkgd_stars, dict) and method is None:
            self.textEdit_backgroundStars.setText('{} background stars added {}'.
                                                  format(len(self.bkgd_stars['x']), method_adverb[method]))

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

    def toggle_convert_im(self):
        # TODO: it's unclear why I (KJB) set this to false. but we want to be able
        # to normalize ITM images so I have commented this out until I better understand it
        # if self.itm:
        #     self.checkBox_normalize.setEnabled(False)
        pass

    def update_segment_guiding_shift(self):
        if self.sender() == self.groupBox_segmentGuiding:
            if self.groupBox_segmentGuiding.isChecked():
                self.radioButton_shifted.setEnabled(os.path.exists(self.shifted_im_file))
            else:
                self.radioButton_shifted.setEnabled(False)

        else:
            if self.radioButton_shifted.isChecked():
                self.lineEdit_regfileSegmentGuiding.setText(self.shifted_guiding_selections_file)
            elif self.radioButton_unshifted.isChecked():
                self.lineEdit_regfileSegmentGuiding.setText(self.guiding_selections_file)

    def update_naming_method(self):
        if self.radioButton_name_manual.isChecked():
            self.stackedWidget.setCurrentIndex(1)
            if self.textEdit_out.toPlainText() == "":
                self.textEdit_out.setEnabled(True)
                self.textEdit_out.setText(OUT_PATH)
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
                self.gs_id, self.apt_guider, self.gs_ra, self.gs_dec = self.query_apt_for_gs(self.program_id, self.observation_num)
            else:
                raise ValueError('Must set both program ID and observation number to use APT')
        elif self.radioButton_name_commissioning.isChecked():
            self.program_id = int(self.lineEdit_commid.text())
            self.observation_num = int(self.comboBox_obs.currentText())
            self.visit_num = 1  # Will we ever have a visit that's not 1?
            self.gs_id, self.apt_guider, self.gs_ra, self.gs_dec = self.query_apt_for_gs(self.program_id, self.observation_num)

        # Check the guider in the APT file matches what's chosen in the GUI
        self.check_guider_against_apt()

        # Update GSID in image normalization
        self.lineEdit_normalize.setText(str(self.gs_id))

    def update_commissioning_name(self):
        # Check which values have been selected already
        valid_practice = self.comboBox_practice.currentText() != '- Select Practice -'
        valid_car = self.comboBox_car.currentText() != '- Select CAR -'

        # When the CAR step is changed...
        if valid_car and self.sender() == self.comboBox_car:
            # Update the observation dropdown box to include the possible observation numbers
            n_obs = int(self.commissioning_dict[self.comboBox_car.currentText().lower()]['observations'])
            self.comboBox_obs.clear()
            self.comboBox_obs.addItem('- Select Obs -')
            for i_obs in range(n_obs):
                self.comboBox_obs.addItem('{:02d}'.format(i_obs + 1))
            for i_obs in np.arange(n_obs, n_obs + 3):
                self.comboBox_obs.addItem('+{:02d}'.format(i_obs + 1))

            # Add/Update the current APT program number
            self.lineEdit_commid.setText(str(self.commissioning_dict[self.comboBox_car.currentText().lower()]['apt']))

        valid_obs = '- Select Obs -' not in self.comboBox_obs.currentText() and self.comboBox_obs.currentText() != ""
        valid_all = valid_practice and valid_car and valid_obs

        # Update which boxes are enabled and disabled accordingly
        self.comboBox_car.setEnabled(valid_practice)
        self.lineEdit_commid.setEnabled(valid_practice)
        self.comboBox_obs.setEnabled(valid_practice & valid_car)
        self.pushButton_commid.setEnabled(valid_practice & valid_car)

        # Update the preview output path
        if valid_all:
            path = os.path.join(SOGS_PATH,
                                self.comboBox_practice.currentText(),
                                self.comboBox_car.currentText().lower().replace('-', ''),
                                'out',
                                'for_obs{:02d}'.format(int(self.comboBox_obs.currentText()))
                                )
            self.textEdit_name_preview.setText(path)
        else:
            self.textEdit_name_preview.setText('')

        # Update file previews
        self.update_filepreview()

        # Update population of guide star information
        if valid_all and any([self.sender() == self.comboBox_car, self.sender() == self.comboBox_obs,
                              self.sender() == self.pushButton_commid]):
            self.update_apt_gs_values()


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

    def open_filename_dialog(self, title, file_type="All Files (*)"):
        """ Dialog box for opening a NIRCam of FGS image"""
        # options = QFileDialog.Options()
        fileName, _ = QFileDialog.getOpenFileName(self, "Open {}".format(title),
                                                  "", file_type)
        # options=options)
        if fileName:
            return fileName

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

            # Load data
            data, _ = utils.get_data_and_header(self.converted_im_file)
            data[data <= 0] = 1

            # Load all_found_psfs*.text
            x, y = [None, None]
            if os.path.exists(self.all_found_psfs_file):
                psf_list = asc.read(self.all_found_psfs_file)
                x = psf_list['x']
                y = psf_list['y']

            # Plot data image and peak locations from all_found_psfs*.txt
            self.canvas_converted.compute_initial_figure(self.canvas_converted.fig, data, x, y)

            # If possible, plot the selected stars in the guiding_selections*.txt
            if os.path.exists(self.guiding_selections_file):
                selected_psf_list = asc.read(self.guiding_selections_file)
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

        # If not, show nothing.
        else:
            # Update textbox showing filepath
            self.textEdit_showingConverted.setText(
                'No converted guider {} image found at {}.'.format(self.buttonGroup_guider.checkedButton().text(),
                                                                   self.converted_im_file))
            self.textEdit_showingConverted.setEnabled(False)

            # Disable the "use converted image" buttons
            self.checkBox_useConvertedImage.setChecked(False)
            self.checkBox_useConvertedImage.setEnabled(False)

            # Disable the "show stars" button
            self.checkBox_showStars.setEnabled(False)

            # Update plot to not show anything
            dummy_img = self.canvas_converted.axes.imshow(
                np.array([[1e4, 1e4], [1e4, 1e4]]), cmap='bone', clim=(1e-1, 1e2)
            )

        return self.canvas_converted.draw()

    def update_shifted_image_preview(self):
        # Are all the necessary fields filled in? If not, don't even try.
        if not self.is_valid_path_defined():
            return

        # Does a shift image exist? If so, show it!
        if os.path.exists(self.shifted_im_file):
            # Prepare to show shifted image
            self.canvas_shifted.axes.set_visible(True)
            self.tabWidget.setCurrentIndex(2)
            self.textEdit_showingShifted.setEnabled(True)

            # Update filepath
            self.textEdit_showingShifted.setText(self.shifted_im_file)

            # Toggle the "use shifted image" buttons
            self.radioButton_shifted.setChecked(True)
            self.lineEdit_regfileSegmentGuiding.setText(self.shifted_guiding_selections_file)

            # Enable the "show stars" button
            self.checkBox_showStars_shifted.setEnabled(True)

            # Load data
            data, _ = utils.get_data_and_header(self.shifted_im_file)
            data[data <= 0] = 1

            # Load all_found_psfs*.text
            x, y = [None, None]
            if os.path.exists(self.shifted_all_found_psfs_file):
                psf_list = asc.read(self.shifted_all_found_psfs_file)
                x = psf_list['x']
                y = psf_list['y']

            # Plot data image and peak locations from all_found_psfs*.txt
            self.canvas_shifted.compute_initial_figure(self.canvas_shifted.fig, data, x, y)

            # If possible, plot the selected stars in the guiding_selections*.txt
            if os.path.exists(self.shifted_guiding_selections_file):
                selected_psf_list = asc.read(self.shifted_guiding_selections_file)
                x_selected = selected_psf_list['x']
                y_selected = selected_psf_list['y']

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

        # If not, show nothing.
        else:
            # Update textbox showing filepath
            self.textEdit_showingShifted.setText(
                'No shifted guider {} image found at {}.'.format(self.buttonGroup_guider.checkedButton().text(),
                                                                   self.shifted_im_file))
            self.textEdit_showingShifted.setEnabled(False)

            # Disable the "use shifted image" buttons
            self.radioButton_unshifted.setChecked(True)
            self.radioButton_shifted.setEnabled(False)

            # Disable the "show stars" button
            self.checkBox_showStars_shifted.setEnabled(False)

            # Update plot to not show anything
            dummy_img = self.canvas_shifted.axes.imshow(
                np.array([[1e4, 1e4], [1e4, 1e4]]), cmap='bone', clim=(1e-1, 1e2)
            )

        return self.canvas_shifted.draw()

    def update_filepreview(self):
        # If either:
        #   1) manual naming is selected and the root, out_dir, and guider have been defined, or
        #   2) commissioning naming is selected and the practice, CAR, and observation have
        #       been selected
        # show an example filepath to a simulated image, and auto-populate the guiding_selections*.txt.text
        # filepath. Also, auto-generate important filenames as class attributes.

        if self.is_valid_path_defined():
            manual_naming = self.radioButton_name_manual.isChecked()
            guider = self.buttonGroup_guider.checkedButton().text()

            # Determine root directory
            if manual_naming:
                root = self.lineEdit_root.text()
                root_dir = os.path.join(self.textEdit_out.toPlainText(), 'out',
                                        self.lineEdit_root.text())
            else:
                root = 'for_obs{:02d}'.format(int(self.comboBox_obs.currentText()))
                root_dir = self.textEdit_name_preview.toPlainText()

            # Set log if not already set (for first file created with MAGIC)
            if self.log is None:
                self.log, self.log_filename = utils.create_logger_from_yaml(__name__, root=root, level='DEBUG')
            # If root is changed, need to create a new log file
            if root != self.log_filename.split('/')[-1].split('masterGUI_')[-1].split('.log')[0]:
                self.log, self.log_filename = utils.create_logger_from_yaml(__name__, root=root, level='DEBUG')

            # Note: maintaining if statements and "old" file names for backwards compatibility.

            # Update guiding selections file path
            guiding_selections_file = os.path.join(
                root_dir, 'guiding_selections_{}_G{}.txt'.format(root, guider)
            )
            guiding_selections_file_old = os.path.join(
                root_dir, '{}_G{}_regfile.txt'.format(root, guider)
            )
            if os.path.exists(guiding_selections_file_old):
                self.guiding_selections_file = guiding_selections_file_old
            else:
                self.guiding_selections_file = guiding_selections_file

            # Update all found PSFs file path
            all_found_psfs_file = os.path.join(
                root_dir, 'all_found_psfs_{}_G{}.txt'.format(root, guider)
            )
            all_found_psfs_file_old = os.path.join(
                root_dir, '{}_G{}_ALLpsfs.txt'.format(root, guider)
            )
            if os.path.exists(all_found_psfs_file_old):
                self.all_found_psfs_file = all_found_psfs_file_old
            else:
                self.all_found_psfs_file = all_found_psfs_file

            # Update converted FGS image filepath
            self.converted_im_file = os.path.join(
                root_dir, 'FGS_imgs', '{}_G{}.fits'.format(root, guider)
            )

            # Update shifted FGS image & catalog filepaths
            self.shifted_im_file = os.path.join(
                root_dir, 'shifted', '{}_G{}.fits'.format(root, guider)
            )
            self.shifted_all_found_psfs_file = os.path.join(
                root_dir, 'shifted', 'all_found_psfs_{}_G{}.txt'.format(root, guider)
            )
            self.shifted_guiding_selections_file = os.path.join(
                root_dir, 'shifted', 'guiding_selections_{}_G{}.txt'.format(root, guider)
            )

            # Update default guiding_selections*.txt paths in GUI
            if os.path.exists(self.guiding_selections_file):
                self.lineEdit_regfileStarSelector.setText(self.guiding_selections_file)
                if not self.radioButton_shifted.isChecked():
                    self.lineEdit_regfileSegmentGuiding.setText(self.guiding_selections_file)
                else:
                    self.lineEdit_regfileSegmentGuiding.setText(self.shifted_guiding_selections_file)
            else:
                self.lineEdit_regfileStarSelector.setText("")
                self.lineEdit_regfileSegmentGuiding.setText("")

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
        observation = observation_list[obs_number - 1]  # indexes from 1
        sr = [x for x in observation.iterchildren() if x.tag.split(namespace_tag)[1] == "SpecialRequirements"][0]

        # Try to pull the Guide Star information
        try:
            gs = [x for x in sr.iterchildren() if x.tag.split(namespace_tag)[1] == "GuideStarID"][0]
        except IndexError:
            self.lineEdit_normalize.setText('')
            LOGGER.error("Master GUI: This observation doesn't have a Guide Star Special Requirement")
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
        data_frame = fgscountrate.query_gsc(gs_id=gs_id)

        # Check there's only 1 line in the GSC with this GS ID
        if len(data_frame) == 1:
            gsc_series = data_frame.iloc[0]
        else:
            self.lineEdit_normalize.setText('')
            LOGGER.error("Master GUI: This Guide Star ID points to multiple lines in catalog")
            raise ValueError("This Guide Star ID points to multiple lines in catalog")

        # Pull RA and DEC
        ra = gsc_series['ra']
        dec = gsc_series['dec']

        LOGGER.info('Master GUI: The Guide Star Catalog have been queried and found RA of {} and DEC of {} '.format(ra, dec))

        return gs_id, guider, ra, dec


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# MAIN FUNCTION
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

def run_MasterGui(root=None, norm_value=12.0, norm_unit="FGS Magnitude", nircam_det=None,
                  nircam=True, smoothing='default', steps=None, in_file=None,
                  bkgd_stars=False, out_dir=OUT_PATH, convert_im=True,
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
