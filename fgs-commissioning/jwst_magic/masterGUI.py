"""GUI interface to JWST MAGIC

The primary interface to JWST MAGIC, from which almost
all of the tool functions can be operated.

Authors
-------
    - Keira Brooks
    - Lauren Chambers

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

import glob
import os
import re
import sys

import matplotlib
import yaml
from PyQt5 import QtCore, uic, QtGui
from PyQt5.QtCore import Qt, pyqtSlot
from PyQt5.QtWidgets import (QApplication, QMainWindow, QMessageBox, QFileDialog,
                             QDialog)

if matplotlib.get_backend() != 'Qt5Agg':
    matplotlib.use('Qt5Agg')  # Make sure that we are using Qt5
from astropy.io import ascii as asc
import numpy as np

from . import run_magic, utils
from jwst_magic import background_stars
from .convert_image import renormalize
from .fsw_file_writer import rewrite_prc
from .segment_guiding import segment_guiding
from .star_selector.SelectStarsGUI import StarClickerMatplotlibCanvas, run_SelectStars

# Define all needed paths
__location__ = os.path.realpath(os.path.join(os.getcwd(), os.path.dirname(__file__)))
PACKAGE_PATH = os.path.dirname(os.path.realpath(__file__))
OUT_PATH = os.path.split(PACKAGE_PATH)[0]  # Location of out/ and logs/ directory
SOGS_PATH = '/data/jwst/wss/guiding/'

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
                 nircam_det=None, nircam=True, global_alignment=False, steps=None,
                 in_file=None, bkgd_stars=False, out_dir=OUT_PATH, convert_im=True,
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

        # Initialize SGT attributes
        self.prognum = None
        self.obsnum = None
        self.visitnum = None

        # Initialize main window object
        QMainWindow.__init__(self)

        # Import .ui file
        uic.loadUi(os.path.join(__location__, 'masterGUI.ui'), self)

        # Set up the custom output stream
        EmittingStream(self.textEdit_log)

        # Create and load GUI session
        self.setWindowTitle('JWST MaGIC')
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
        self.lineEdit_root.editingFinished.connect(self.update_filepreview)
        self.pushButton_root.clicked.connect(self.on_click_root)
        self.pushButton_out.clicked.connect(self.on_click_out)
        self.textEdit_out.installEventFilter(self)

        # Image convertor widgets
        self.pushButton_backgroundStars.clicked.connect(self.on_click_bkgdstars)
        self.horizontalSlider_coarsePointing.sliderReleased.connect(self.on_change_jitter)
        self.lineEdit_coarsePointing.editingFinished.connect(self.on_change_jitter)
        self.checkBox_useConvertedImage.toggled.connect(self.toggle_convert_im)

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

        # # FOR TESTING ONLY
        # if not utils.on_sogs_network():
        #     sogs_dirs = ['wfp_january_2019']

        # If on SOGS, pull out practice names from existing practice directories
        else:
            sogs_search = os.path.join(SOGS_PATH, '*')
            sogs_dirs = [os.path.basename(dir) for dir in glob.glob(sogs_search)]
            for d in ['data', 'processing', 'MAGIC_logs']:
                sogs_dirs.remove(d)

        # Load all OTE cars from commissioning activities YAML file
        commissioning_yaml = os.path.join(__location__, 'data', 'commissioning_activities.yaml')
        with open(commissioning_yaml, encoding="utf-8") as f:
            self.commissioning_dict = yaml.load(f.read())
        cars_list = list(self.commissioning_dict.keys())

        # Use to populate practice and CAR dropdown boxes
        for practice in sogs_dirs:
            self.comboBox_practice.addItem(practice)
        for car in cars_list:
            self.comboBox_car.addItem(car.upper())

        # Connect combo boxes to one another
        self.comboBox_practice.currentIndexChanged.connect(self.update_commissioning_name)
        self.comboBox_car.currentIndexChanged.connect(self.update_commissioning_name)
        self.comboBox_obs.currentIndexChanged.connect(self.update_commissioning_name)


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
        if self.lineEdit_inputImage.text() == "":
            self.no_inputImage_dialog()
            return
        if not self.buttonGroup_guider.checkedButton():
            self.no_guider_dialog()
            return

        # General input
        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        input_image = self.lineEdit_inputImage.text()
        guider = int(self.buttonGroup_guider.checkedButton().text())
        if self.radioButton_name_manual.isChecked():
            root = self.lineEdit_root.text()
            out_dir = self.textEdit_out.toPlainText().rstrip()
        elif self.radioButton_name_commissioning.isChecked():
            root = 'for_obs{}'.format(self.comboBox_obs.currentText())
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
        norm_value = float(self.lineEdit_normalize.text())
        norm_unit = self.comboBox_normalize.currentText()
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
        global_alignment = self.checkBox_globalAlignment.isChecked()
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
        crowded_field =  self.radioButton_crowded_id_attitude.isChecked()

        # Rewrite .prc and regfile.txt ONLY
        if self.checkBox_rewritePRC.isChecked():
            # Open converted FGS file
            data, _ = utils.get_data_and_header(self.converted_im_file)

            # Update array for showing in LogNorm
            data[data <= 0] = 1

            # Open ALLpsfs.txt (list of all identified segments)
            all_psfs = self.all_found_psfs_file
            all_rows = asc.read(all_psfs)
            x = all_rows['x'].data
            y = all_rows['y'].data

            # Run the select stars GUI to determine the new orientation
            inds = run_SelectStars(data, x, y, 20, masterGUIapp=self.app)

            # Rewrite the id.prc and acq.prc files
            rewrite_prc.rewrite_prc(inds, guider, root, out_dir)
            print("** Run Complete **\n\n")

            # Update converted image preview
            self.update_filepreview()
            return

        # Segment guiding
        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        if self.groupBox_segmentGuiding.isChecked():
            # Get APT program information from parsed header
            self.parse_header(input_image)
            program_id = self.prognum
            observation_num = self.obsnum
            visit_num = self.visitnum

            # Check if this is a photometry only override file or segment override file
            if self.radioButton_photometryOverride.isChecked():
                # Generate the file
                segment_guiding.generate_photometry_override_file(
                    root, program_id, observation_num, visit_num, out_dir=out_dir
                )
            else:
                # Define location of ALLpsfs catalog file
                if self.radioButton_shifted.isChecked():
                    segment_infile = self.shifted_all_found_psfs_file
                else:
                    segment_infile = self.all_found_psfs_file

                # Verify that the ALLpsfs.txt file exists
                if not os.path.exists(segment_infile):
                    raise OSError('Provided segment infile {} not found.'.format(segment_infile))

                # Determine which image to use and load it
                if self.radioButton_shifted.isChecked():
                    fgs_filename = self.shifted_im_file
                elif self.radioButton_unshifted.isChecked() and os.path.exists(self.converted_im_file):
                    fgs_filename = self.converted_im_file
                elif self.radioButton_unshifted.isChecked():
                    fgs_filename = input_image
                data, _ = utils.get_data_and_header(fgs_filename)

                # Determine whether to load regfile or run GUI
                GUI = not self.radioButton_regfileSegmentGuiding.isChecked()
                selected_segs = self.lineEdit_regfileStarSelector.text()

                # Run the tool and generate the file
                segment_guiding.generate_segment_override_file(
                    segment_infile, guider, program_id, observation_num, visit_num,
                    root=root, out_dir=out_dir, selected_segs=selected_segs,
                    click_to_select_gui=GUI, data=data, master_gui_app=self.app
                )

            print("** Run Complete **")

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
                              global_alignment=global_alignment,
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
            print("** Run Complete **")

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
        root = utils.make_root(self.root_default, self.lineEdit_inputImage.text())
        # self.lineEdit_root.setText(root)
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
            norm_obj = renormalize.NormalizeToCountrate(norm_value, norm_unit, guider)
            fgs_countrate = norm_obj.to_countrate()
            jmag = renormalize.fgs_countrate_to_j_mag(fgs_countrate, guider)
        # If not, determine the FGS counts of the input image
        else:
            input_image = self.lineEdit_inputImage.text()
            data, _ = utils.get_data_and_header(input_image)
            fgs_countrate = np.sum(data[data > np.median(data)])
            jmag = renormalize.fgs_countrate_to_j_mag(fgs_countrate, guider)

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
        if self.sender == self.checkBox_showStars:
            for line in self.converted_im_circles:
                line[0].set_visible(show)
            self.canvas_converted.peaks.set_visible(show)

            self.canvas_converted.draw()

        elif self.sender == self.checkBox_showStars_shifted:
            for line in self.shifted_im_circles:
                line[0].set_visible(show)
            self.canvas_shifted.peaks.set_visible(show)

            self.canvas_shifted.draw()

    def toggle_convert_im(self):
        if self.itm:
            self.checkBox_normalize.setEnabled(False)

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

    def update_commissioning_name(self):
        # Check which values have been selected already
        valid_practice = self.comboBox_practice.currentText() != '- Select Practice -'
        valid_car = self.comboBox_car.currentText() != '- Select CAR -'
        valid_obs = self.comboBox_obs.currentText() != '- Select Obs -'

        # When the CAR step is changed...
        if valid_car and self.sender() == self.comboBox_car:
            # Update the observation dropdown box to include the possible observation numbers
            obs_list = self.commissioning_dict[self.comboBox_car.currentText().lower()]['observations']
            self.comboBox_obs.clear()
            self.comboBox_obs.addItem('- Select Obs -')
            for i_obs in range(int(obs_list)):
                self.comboBox_obs.addItem('{:02d}'.format(i_obs + 1))
            valid_obs = False

            # Add the current APT program number
            self.label_apt.setText('APT: {}'.format(
                self.commissioning_dict[self.comboBox_car.currentText().lower()]['apt'])
            )

        # Update which boxes are enabled and disabled accordingly
        self.comboBox_car.setEnabled(valid_practice)
        self.comboBox_obs.setEnabled(valid_practice & valid_car)

        # Update the preview output path
        if valid_practice and valid_car and valid_obs:
            path = os.path.join(SOGS_PATH,
                                self.comboBox_practice.currentText(),
                                self.comboBox_car.currentText().lower().replace('-', ''),
                                'out',
                                'for_obs{}'.format(self.comboBox_obs.currentText())
                                )
            self.textEdit_name_preview.setText(path)
        else:
            self.textEdit_name_preview.setText('')

        # Update file previews
        self.update_filepreview()


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

        # Set defaults from parsed header
        SGT_dialog.lineEdit_programNumber.setText(self.prognum)
        SGT_dialog.lineEdit_observationNumber.setText(self.obsnum)
        SGT_dialog.lineEdit_visitNumber.setText(self.visitnum)

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
                self.checkBox_normalize.setEnabled(False)
            else:
                self.itm = False
        except KeyError:
            pass

        # Parse APT program, observation, and visit information
        keywords = ['PROGRAM', 'OBSERVTN', 'VISIT']
        attributes = ['prognum', 'obsnum', 'visitnum']
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
            if self.textEdit_name_preview.toPlainText() == "":
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

            # Load ALLpsfs.text
            x, y = [None, None]
            if os.path.exists(self.all_found_psfs_file):
                psf_list = asc.read(self.all_found_psfs_file)
                x = psf_list['x']
                y = psf_list['y']

            # Plot data image and peak locations from ALLpsfs.txt
            self.canvas_converted.compute_initial_figure(self.canvas_converted.fig, data, x, y)

            # If possible, plot the selected stars in the regfile.txt
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

            # Enable the "show stars" button
            self.checkBox_showStars_shifted.setEnabled(True)

            # Load data
            data, _ = utils.get_data_and_header(self.shifted_im_file)
            data[data <= 0] = 1

            # Load ALLpsfs.text
            x, y = [None, None]
            if os.path.exists(self.shifted_all_found_psfs_file):
                psf_list = asc.read(self.shifted_all_found_psfs_file)
                x = psf_list['x']
                y = psf_list['y']

            # Plot data image and peak locations from ALLpsfs.txt
            self.canvas_shifted.compute_initial_figure(self.canvas_shifted.fig, data, x, y)

            # If possible, plot the selected stars in the regfile.txt
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
        # show an example filepath to a simulated image, and auto-populate the regfile.text
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
                root = 'for_obs{}'.format(self.comboBox_obs.currentText())
                root_dir = self.textEdit_name_preview.toPlainText()

            # Note: maintaining if statements and "old" file names for backwards compatibility.

            # Update guiding selections file path
            guiding_selections_file = os.path.join(root_dir,
                                                    'guiding_selections_{}_G{}.txt'.format(root, guider))
            guiding_selections_file_old = os.path.join(root_dir,
                                                       '{}_G{}_regfile.txt'.format(root, guider))
            if os.path.exists(guiding_selections_file_old):
                self.guiding_selections_file = guiding_selections_file_old
            else:
                self.guiding_selections_file = guiding_selections_file

            # Update all found PSFs file path
            all_found_psfs_file = os.path.join(root_dir,
                                               'all_found_psfs_{}_G{}.txt'.format(root, guider))
            all_found_psfs_file_old = os.path.join(root_dir,
                                                   '{}_G{}_ALLpsfs.txt'.format(root, guider))
            if os.path.exists(all_found_psfs_file_old):
                self.all_found_psfs_file = all_found_psfs_file_old
            else:
                self.all_found_psfs_file = all_found_psfs_file


            # Update converted FGS image filepath
            self.converted_im_file = os.path.join(root_dir, 'FGS_imgs',
                                                  '{}_G{}.fits'.format(root, guider))

            # Update shifted FGS image & catalog filepaths
            self.shifted_im_file = os.path.join(root_dir, 'shifted',
                                                '{}_G{}.fits'.format(root, guider))
            self.shifted_all_found_psfs_file = os.path.join(root_dir, 'shifted',
                                                'all_found_psfs_{}_G{}.txt'.format(root, guider))
            self.shifted_guiding_selections_file = os.path.join(root_dir, 'shifted',
                                                'guiding_selections_{}_G{}.txt'.format(root, guider))

            # Update default regfile paths in GUI
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


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# MAIN FUNCTION
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

def run_MasterGui(root=None, norm_value=12.0, norm_unit="FGS Magnitude", nircam_det=None,
                  nircam=True, global_alignment=False, steps=None, in_file=None,
                  bkgd_stars=False, out_dir=OUT_PATH, convert_im=True,
                  star_selection=True, star_selection_gui=True, file_writer=True,
                  segment_guiding=False, itm=False):
    # RUN GUI
    app = QtCore.QCoreApplication.instance()  # Use existing instance, if there is one
    if app is None:
        app = QApplication(sys.argv)

    ex = MasterGui(root, norm_value, norm_unit, nircam_det, nircam, global_alignment, steps,
                   in_file, bkgd_stars, out_dir, convert_im, star_selection_gui,
                   file_writer, segment_guiding, app=app, itm=itm)
    # #return ex.settings
    out = app.exec()

    return out
