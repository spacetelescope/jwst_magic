"""GUI interface to FGS Commissioning Tools

~Description goes here~

Authors
-------
    - Keira Brooks

Use
---
This GUI can be run in the python shell as such:
    ::
    from jwst_fgs_commissioning_tools.master_gui import masterGUI
    masterGUI.run_MasterGui()

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
FGS Commissioning Tools package, be sure to use the existing instance
of QApplication (access it at QtCore.QCoreApplication.instance()) when
calling the QApplication instance to run a window/dialog/GUI.
"""

import os
import sys
import inspect

from PyQt5.QtWidgets import (QApplication, QWidget, QPushButton, QInputDialog,
                             QLineEdit, QMainWindow, QAction, QMessageBox,
                             QFileDialog, QLabel, QCheckBox, QFrame, QSizePolicy,
                             QHBoxLayout, QSplitter, QTextEdit, QStyleFactory,
                             QSpinBox, QComboBox, QRadioButton, QGridLayout,
                             QGroupBox, qApp, QButtonGroup)
from PyQt5 import QtGui, QtCore
from PyQt5.QtGui import QIcon, QFont
from PyQt5.QtCore import pyqtSlot, Qt, QSettings, QTimer
import matplotlib
if matplotlib.get_backend() != 'Qt5Agg':
    matplotlib.use('Qt5Agg')  # Make sure that we are using Qt5

from jwst_fgs_commissioning_tools import run_fgs_commissioning_tool
from jwst_fgs_commissioning_tools.convert_image import counts_to_jmag
from jwst_fgs_commissioning_tools import utils

MASTERGUI_PATH = os.path.dirname(os.path.realpath(__file__))
PACKAGE_PATH = os.path.split(MASTERGUI_PATH)[0]
OUT_PATH = os.path.split(PACKAGE_PATH)[0]  # Location of out/ and logs/ directory

class MasterGui(QMainWindow):

    def __init__(self, root=None, fgs_counts=None, jmag=11.0,
                 nircam_det=None, nircam=True, global_alignment=False, steps=None,
                 in_file=None, bkgd_stars=False, out_dir=OUT_PATH, convert_im=True,
                 star_selection=True, star_selection_gui=True, file_writer=True,
                 segment_guiding=False, app=None):

        # Initialize main window object (?!?!)
        QMainWindow.__init__(self)

        self.title = 'Master GUI - FGS Commissioning Tools'
        self.image_dim = 1000
        self.app = app

        # Set stages
        self.convert_im = convert_im
        self.star_selection = star_selection
        self.file_writer = file_writer
        self.segment_guiding = segment_guiding

        # Set defaults
        self.root_default = root
        self.fgs_counts_default = fgs_counts
        self.jmag_default = jmag
        self.nircam_det_default = nircam_det
        self.nircam_default = nircam
        self.global_alignment_default = global_alignment
        if not steps:
            self.steps_default = ['ID', 'ACQ1', 'ACQ2', 'TRK']
        self.bkgd_stars_default = bkgd_stars
        self.star_selection_gui_default = star_selection_gui
        self.in_file_default = in_file
        self.out_dir_default = out_dir

        self.initUI()

    def initUI(self):
        ''' Set up GUI window
        Create boxes for each section and cretae the "run" and "exit" buttons
        '''
        self.setWindowTitle(self.title)
        self.main_widget = QWidget(self)
        mainGrid = QGridLayout()  # set grid layout
        self.main_widget.setLayout(mainGrid)
        self.main_widget.setFocus()
        self.setCentralWidget(self.main_widget)

        # Set the font format for titles
        self.title_font = QFont()
        self.title_font.setBold(True)

        ### Build the MAIN section ---------------------------------------------
        inputGroupBox = self.create_input_section()
        mainGrid.addWidget(inputGroupBox, 0, 0, 1, 2)

        # ### Build the CONVERT IMAGE section ----------------------------------
        convertGroupBox = self.create_convert_section()
        mainGrid.addWidget(convertGroupBox, 1, 0, 1, 2)


        # ### Build the STAR SELECTION section ---------------------------------
        starGroupBox = self.create_star_section()
        mainGrid.addWidget(starGroupBox, 2, 0, 1, 2)


        # ### Build the FSW FILE WRITER section --------------------------------
        filewriterGroupBox = self.create_filewriter_section()
        mainGrid.addWidget(filewriterGroupBox, 3, 0, 1, 2)

        # ###  Build the SEGMENT GUIDING section -------------------------------
        segmentGroupBox = self.create_segment_section()
        mainGrid.addWidget(segmentGroupBox, 4, 0, 1, 2)

        # ###  Save out inputs and run tool ------------------------------------
        button_run = QPushButton("Run", self)
        button_run.resize(button_run.minimumSizeHint())
        button_run.clicked.connect(self.run_tool)
        mainGrid.addWidget(button_run, 6, 0, 1, 1)

        button_done = QPushButton("Exit", self)
        button_done.resize(button_done.minimumSizeHint())
        button_done.clicked.connect(self.close_application)
        mainGrid.addWidget(button_done, 6, 1, 1, 1)

        # Show the GUI ---------------------------------------------------------
        self.show()

    def create_input_section(self):
        '''
        Create the first box: Standard inputs

        This section includes:
        *Input Image
        *Guider
        *Root
        *Out directory
        *Non-standard psfs button
        *Background stars button
        '''

        inputGroupBox = QGroupBox('Standard Inputs', self)
        inputGrid = QGridLayout()
        inputGroupBox.setLayout(inputGrid)

        # Input Image
        inputGrid.addWidget(QLabel('Input Image (required)', self), 1, 0)
        self.textbox_input_im = QLineEdit(self)
        self.textbox_input_im.setMinimumSize(400, 20)
        inputGrid.addWidget(self.textbox_input_im, 2, 0)

        self.button_input_image = QPushButton('Open', self)
        self.button_input_image.setToolTip('Open the input image')
        self.button_input_image.clicked.connect(self.on_click_input)
        inputGrid.addWidget(self.button_input_image, 2, 1)

        # Guider
        inputGrid.addWidget(QLabel('Guider (required)', self), 3, 0)
        self.cb_guider = QComboBox(self)
        self.cb_guider.addItem("-Select-")
        self.cb_guider.addItem("1")
        self.cb_guider.addItem("2")
        self.cb_guider.activated.connect(self.update_filepreview)
        inputGrid.addWidget(self.cb_guider, 4, 0, 1, 2)

        # Root
        inputGrid.addWidget(QLabel('Root', self), 5, 0)
        self.textbox_root = QLineEdit(self)
        self.textbox_root.editingFinished.connect(self.update_filepreview)
        inputGrid.addWidget(self.textbox_root, 6, 0, 1, 2)

        # Out directory
        inputGrid.addWidget(QLabel('Out Directory - the directory where files are saved', self), 7, 0)
        self.textbox_outdir = QLineEdit(self)
        self.textbox_outdir.setText(self.out_dir_default)
        self.textbox_outdir.setEnabled(False)
        inputGrid.addWidget(self.textbox_outdir, 8, 0)
        self.button_input_image = QPushButton('Open', self)
        self.button_input_image.setToolTip('Open the out directory')
        self.button_input_image.clicked.connect(self.on_click_out)
        inputGrid.addWidget(self.button_input_image, 8, 1)

        # Filepath preview
        inputGrid.addWidget(QLabel('Example output filepath', self), 9, 0)
        self.textbox_filepreview = QTextEdit(self)
        self.textbox_filepreview.setStyleSheet("QTextEdit { background-color : #ECECEC }");
        self.textbox_filepreview.setMaximumSize(2000, 40)
        inputGrid.addWidget(self.textbox_filepreview, 10, 0, 1, 2)

        ## Global Alignment flag
        self.cb_ga = QCheckBox('Non-standard PSFs (i.e. Global alignment, image array)', self)
        self.cb_ga.setCheckable(True)
        if self.global_alignment_default:
            self.cb_ga.toggle()
        inputGrid.addWidget(self.cb_ga, 11, 0)

        # Background stars?
        self.cb_bkgd_stars = QCheckBox('Background stars', self)
        self.cb_bkgd_stars.setCheckable(True)
        if self.bkgd_stars_default:
            self.cb_bkgd_stars.toggle()
        inputGrid.addWidget(self.cb_bkgd_stars, 12, 0)

        return inputGroupBox

    def create_convert_section(self):
        '''
        Create the second box: Convert Image
        Take a FGS or NIRCam image and convert and renormalize it to create a
        RAW FGS image

        This section includes:
        *NIRCam vs FGS Image
        *NIRCAm Detector (only for if this is a NIRCam image)
        *FGS Counts vs J mag w/ place to inlude counts or J mag for renormalization
        '''
        convertGroupBox = QGroupBox('Convert Input Image to an FGS Raw Image', self)
        convertGrid = QGridLayout()
        convertGroupBox.setLayout(convertGrid)

        # # If we want to convert image
        self.cb_convert_im = QCheckBox('Convert Image', self)
        self.cb_convert_im.setCheckable(True)
        self.cb_convert_im.setFont(self.title_font)
        convertGrid.addWidget(self.cb_convert_im, 1, 0)
        if self.convert_im:
            self.cb_convert_im.toggle()

        instrument_group = QButtonGroup(self.main_widget) # Number group
        self.rb_nircam = QRadioButton("NIRCam Image")
        instrument_group.addButton(self.rb_nircam)
        if self.nircam_default:
            self.rb_nircam.setChecked(True)
        self.rb_nircam.toggled.connect(lambda: self.btnstate(self.rb_nircam))
        convertGrid.addWidget(self.rb_nircam, 2, 0)

        self.rb_fgs = QRadioButton("FGS Image")
        instrument_group.addButton(self.rb_fgs)
        if not self.nircam_default:
            self.rb_fgs.setChecked(True)
            self.cb_nircam_det  = None
        self.rb_fgs.toggled.connect(lambda: self.btnstate(self.rb_fgs))
        convertGrid.addWidget(self.rb_fgs, 2, 1)

        # Set NIRCam detector
        convertGrid.addWidget(QLabel('NIRCam Detector:', self), 3, 0)
        self.cb_nircam_det = QComboBox(self)
        self.cb_nircam_det.addItem("-Select-")
        self.cb_nircam_det.addItem("A1")
        self.cb_nircam_det.addItem("A2")
        self.cb_nircam_det.addItem("A3")
        self.cb_nircam_det.addItem("A4")
        self.cb_nircam_det.addItem("A5 - long wave")
        self.cb_nircam_det.addItem("B1")
        self.cb_nircam_det.addItem("B2")
        self.cb_nircam_det.addItem("B3")
        self.cb_nircam_det.addItem("B4")
        self.cb_nircam_det.addItem("B5 - long wave")

        if self.nircam_det_default:
            index = self.cb_nircam_det.findText(self.nircam_det_default, Qt.MatchFixedString)
            if index >= 0:
                self.cb_nircam_det.setCurrentIndex(index)

        convertGrid.addWidget(self.cb_nircam_det, 3, 1)

        #Normalize by FGS counts or JMAG
        magnitude_group = QButtonGroup(self.main_widget)
        self.rb_fgs_counts = QRadioButton("FGS Counts")
        magnitude_group.addButton(self.rb_fgs_counts)
        self.rb_fgs_counts.toggled.connect(lambda: self.btnstate(self.rb_fgs_counts))
        convertGrid.addWidget(self.rb_fgs_counts, 4, 0)

        self.rb_jmag = QRadioButton("J mag")
        magnitude_group.addButton(self.rb_jmag)
        self.rb_jmag.toggled.connect(lambda: self.btnstate(self.rb_jmag))
        convertGrid.addWidget(self.rb_jmag, 4, 1)

        convertGrid.addWidget(QLabel('Value:', self), 4, 3)
        self.textbox_value = QLineEdit(self)
        if self.jmag_default:
            self.rb_jmag.setChecked(True)
            self.textbox_value.setText(str(self.jmag_default))
        elif self.fgs_counts_default:
            self.rb_fgs_counts.setChecked(True)
            self.textbox_value.setText(str(self.fgs_counts_default))
        elif not self.jmag_default and not self.fgs_counts_default:
            self.rb_jmag.setChecked(True)
            self.textbox_value.setText(str(self.jmag_default))

        convertGrid.addWidget(self.textbox_value, 4, 4)

        # Set toggling
        self.cb_convert_im.toggled.connect(self.on_check_convertimage)

        return convertGroupBox

    def create_star_section(self):
        '''
        Create the third box: Do Star Selection
        User decides if star selection will happen from pre-made file ("in/reg file")
        or by opening the Star Selection GUI

        This section includes:
        *Star Selection GUI - User chooses own guide and reference stars
        *In/Reg file - file giving x,y positions and count rates of guide and
            reference stars (in that order)
        '''
        # If we want to use the star slection GUI
        starGroupBox = QGroupBox('Guide and Reference Star Selection', self)
        starGrid = QGridLayout()
        starGroupBox.setLayout(starGrid)

        self.cb_star_selection = QCheckBox('Complete Star Selection', self)
        self.cb_star_selection.setCheckable(True)
        self.cb_star_selection.setFont(self.title_font)
        starGrid.addWidget(self.cb_star_selection, 1, 0)
        if self.cb_star_selection:
            self.cb_star_selection.toggle()

        ## Run the star selection GUI
        # Add checkbox
        self.cb_star_selection_gui = QCheckBox('Star Selection GUI', self)
        self.cb_star_selection_gui.setCheckable(True)
        self.cb_star_selection_gui.toggled.connect(self.on_check_starselectgui)
        starGrid.addWidget(self.cb_star_selection_gui, 2, 0)

        ## Pass in an in or reg file
        # Add checkbox
        self.cb_infile = QCheckBox('In/Reg file', self)
        self.cb_infile.setCheckable(True)
        self.cb_infile.toggled.connect(self.on_check_infile)
        starGrid.addWidget(self.cb_infile, 3, 0)
        # Add textbox
        self.textbox_infile = QLineEdit(self)
        self.textbox_infile.setMinimumSize(200, 20)
        starGrid.addWidget(self.textbox_infile, 3, 2)
        # Add 'Open' button
        self.button_input_infile = QPushButton('Open', self)
        self.button_input_infile.setToolTip('Open the in_file directory')
        self.button_input_infile.clicked.connect(self.on_click_infile)
        starGrid.addWidget(self.button_input_infile, 3, 5)

        self.cb_star_selection.toggled.connect(self.on_check_starselect)

        #self.cb_star_selection.toggled.connect(self.on_check_starselectgui)
        # Check the box that is set to default
        if self.in_file_default:
            self.cb_infile.setChecked(True)
            self.textbox_infile.setEnabled(True)
            self.textbox_infile.setText(self.in_file_default)
            self.button_input_infile.setEnabled(True)
            self.cb_star_selection_gui.setChecked(False)
        else:
            self.cb_star_selection_gui.setChecked(True)
            self.cb_infile.setChecked(False)
            self.textbox_infile.setEnabled(False)
            self.button_input_infile.setEnabled(False)


        return starGroupBox

    def create_filewriter_section(self):
        '''
        Create the fourth box: FSW file writer
        Create all files necessary for running this image through flight software
        simulators

        This section includes:
        *ID/ACQ1/ACQ2/TRK/FG - steps to create files for
            *ID - Identification
            *ACQ1 - Acquisition 1
            *ACQ2 - Acquisition 2
            *TRK - Track
            *FG - Fine guide
        '''
        # Use the FSW file writer
        filewriterGroupBox = QGroupBox('Output FSW files', self)
        filewriterGrid = QGridLayout()
        filewriterGroupBox.setLayout(filewriterGrid)

        # Write the files out
        self.cb_file_writer = QCheckBox('FSW File Writer', self)
        self.cb_file_writer.setCheckable(True)
        self.cb_file_writer.toggled.connect(self.on_check_filewriter)

        self.cb_file_writer.setFont(self.title_font)
        filewriterGrid.addWidget(self.cb_file_writer, 1, 0)
        # Which steps to write out
        self.cb_id = QCheckBox('ID', self)
        # if 'ID' in steps:
        #     self.cb_id.
        filewriterGrid.addWidget(self.cb_id, 2, 1)
        self.cb_acq1 = QCheckBox('ACQ1', self)
        filewriterGrid.addWidget(self.cb_acq1, 2, 2)
        self.cb_acq2 = QCheckBox('ACQ2', self)
        filewriterGrid.addWidget(self.cb_acq2, 2, 3)
        self.cb_trk = QCheckBox('TRK', self)
        filewriterGrid.addWidget(self.cb_trk, 2, 4)
        self.cb_fg = QCheckBox('FG', self)
        filewriterGrid.addWidget(self.cb_fg, 2, 5)

        if self.file_writer:
            self.cb_file_writer.toggle()

        return filewriterGroupBox

    def create_segment_section(self):
        '''
        Create the fifth box: Segment Guiding
        Will pop open separate GUI for setting the parameters for running the
        Segment Guiding Tool and creating the segment override file.

        '''
        # Do segment guiding
        segmentGroupBox = QGroupBox('Segment Guiding', self)
        segmentGrid = QGridLayout()
        segmentGroupBox.setLayout(segmentGrid)

        self.cb_segment_guiding = QCheckBox('Segment Guiding', self)
        self.cb_segment_guiding.setCheckable(True)
        if self.segment_guiding:
            self.cb_segment_guiding.toggle()

        self.cb_segment_guiding.setFont(self.title_font)
        segmentGrid.addWidget(self.cb_segment_guiding, 1, 0)

        self.cb_segment_guiding.toggled.connect(self.on_check_segment_guiding)

        return segmentGroupBox


    def close_application(self):
        ''' Close the window'''
        self.close()
        self.app.quit()

    def close_dialog(self):
        ''' Needed to make sure somethign happens when a button is clicked, but
        does nothing.'''
        pass

    def run_tool(self):
        '''
        Takes inputs provided by user and runs run_fgs_commissioning_tool
        '''
        # Required
        if self.textbox_input_im.text() == "":
            no_input_im_dialog = QMessageBox()
            no_input_im_dialog.setText('No input image' + ' ' * 50)
            no_input_im_dialog.setInformativeText('The tool will not be able to continue. Please pass in an input image.')
            no_input_im_dialog.setStandardButtons(QMessageBox.Ok)
            no_input_im_dialog.buttonClicked.connect(self.close_dialog)
            no_input_im_dialog.exec()
            return

        if str(self.cb_guider.currentText()) == '-Select-':
            no_guider_dialog = QMessageBox()
            no_guider_dialog.setText('No guider choosen' + ' ' * 50)
            no_guider_dialog.setInformativeText('The tool will not be able to continue. Please select Guider 1 or 2.')
            no_guider_dialog.setStandardButtons(QMessageBox.Ok)
            no_guider_dialog.buttonClicked.connect(self.close_dialog)
            no_guider_dialog.exec()
            return

        input_image = self.textbox_input_im.text()
        guider = int(self.cb_guider.currentText())

        root = self.textbox_root.text()
        out_dir = self.textbox_outdir.text()

        bkgd_stars = self.cb_bkgd_stars.isChecked()
        global_alignment = self.cb_ga.isChecked()

        #Convert image
        convert_im = self.cb_convert_im.isChecked()
        nircam = self.rb_nircam.isChecked()
        nircam_det = str(self.cb_nircam_det.currentText())
        if self.rb_fgs_counts.isChecked():
            fgs_counts = float(self.textbox_value.text())
        else:
            fgs_counts = None
        if self.rb_jmag.isChecked():
            jmag = float(self.textbox_value.text())
        else:
            jmag = None

        #Star selection
        star_selection = self.cb_star_selection.isChecked()
        star_selectiongui = self.cb_star_selection_gui.isChecked()
        if not self.cb_infile.isChecked():
            in_file = None
        else:
            in_file = self.textbox_infile.text()

        #File writer
        file_writer = self.cb_file_writer.isChecked()
        steps = []
        if self.cb_id.isChecked():
            steps.append('ID')
        if self.cb_acq1.isChecked():
            steps.append('ACQ1')
        if self.cb_acq2.isChecked():
            steps.append('ACQ2')
        if self.cb_trk.isChecked():
            steps.append('LOSTRK')
        if self.cb_fg.isChecked():
            steps.append('FG')


        segment_guiding = self.cb_segment_guiding.isChecked()

        if convert_im or star_selection or file_writer:
            run_fgs_commissioning_tool.run_all(input_image, guider, root,
                                               fgs_counts, jmag, nircam_det,
                                               nircam, global_alignment,
                                               steps, in_file, bkgd_stars,
                                               out_dir, convert_im, star_selection,
                                               star_selectiongui, file_writer, self.app)
            print("** Run Complete **")
        if segment_guiding:
            pass

    # What to do when specific buttons are pressed -----------------------------
    def on_click_input(self):
        ''' Using the Input Image Open button (open file) '''
        filename = self.open_filename_dialog("NIRCam or FGS image")
        self.textbox_input_im.setText(filename)
        if self.textbox_root.text() == "":
            root = utils.make_root(self.root_default, self.textbox_input_im.text())
            self.textbox_root.setText(root)
            if self.textbox_outdir.text() == "":
                self.textbox_outdir.setEnabled(True)
                self.textbox_outdir.setText(OUT_PATH)
            if self.textbox_outdir.text() == OUT_PATH:
                self.textbox_outdir.setEnabled(True)

        self.update_filepreview()
        return filename

    def on_click_out(self):
        ''' Using the Out Dir Open button (open directory) '''
        dirname = self.open_dirname_dialog()
        self.textbox_outdir.setEnabled(True)
        self.textbox_outdir.setText(dirname)
        self.update_filepreview()
        return dirname

    def on_click_infile(self):
        ''' Using the Infile Open button (open file) '''
        filename = self.open_filename_dialog('In/Reg file')
        self.textbox_infile.setText(filename)
        return filename

    # What to do when specific check boxes are checked -------------------------
    def on_check_convertimage(self, is_toggle):
        ''' Checking the convert image box - defaults set here'''
        if is_toggle:
            self.rb_nircam.setEnabled(True)
            self.rb_fgs.setEnabled(True)
            self.cb_nircam_det.setEnabled(True)
            self.rb_fgs_counts.setEnabled(True)
            self.rb_jmag.setEnabled(True)
            self.textbox_value.setEnabled(True)

        else:
            self.rb_nircam.setEnabled(False)
            self.rb_fgs.setEnabled(False)
            self.cb_nircam_det.setEnabled(False)
            self.rb_fgs_counts.setEnabled(False)
            self.rb_jmag.setEnabled(False)
            self.textbox_value.setEnabled(False)


        if not self.rb_nircam.isChecked():
            self.rb_nircam.toggle()
        if not self.rb_jmag.isChecked():
            self.rb_jmag.toggle()

    def on_check_starselect(self, is_toggle):
        ''' Checking that the star selection step occurs box '''
        if is_toggle:
            self.cb_star_selection_gui.setEnabled(True)
            self.cb_infile.setEnabled(True)
            self.textbox_infile.setEnabled(True)
            self.button_input_infile.setEnabled(True)

        else:
            self.cb_star_selection_gui.setEnabled(False)
            self.cb_infile.setEnabled(False)
            self.textbox_infile.setEnabled(False)
            self.button_input_infile.setEnabled(False)


    def on_check_starselectgui(self, is_toggle):
        ''' Checking the star selection GUI box '''
        if is_toggle:
            self.textbox_infile.setEnabled(False)
            self.button_input_infile.setEnabled(False)
            self.cb_infile.setChecked(False)
        else:
            self.textbox_infile.setEnabled(True)
            self.button_input_infile.setEnabled(True)
            self.cb_infile.setChecked(True)

    def on_check_infile(self, is_toggle):
        ''' Checking the star selection infile box '''
        if is_toggle:
            self.textbox_infile.setEnabled(True)
            self.button_input_infile.setEnabled(True)
            if self.cb_star_selection_gui.isChecked():
                self.cb_star_selection_gui.toggle()
        else:
            if not self.cb_star_selection_gui.isChecked():
                self.cb_star_selection_gui.toggle()

    def on_check_filewriter(self, is_toggle):
        ''' Checking the filewriter box - defaults set here'''
        idstep = bool('ID' in self.steps_default)
        acq1step = bool('ACQ1' in self.steps_default)
        acq2step = bool('ACQ2' in self.steps_default)
        trkstep = bool('TRK' in self.steps_default)
        fgstep = bool('FG' in self.steps_default)

        if is_toggle:
            self.cb_id.setEnabled(True)
            self.cb_acq1.setEnabled(True)
            self.cb_acq2.setEnabled(True)
            self.cb_trk.setEnabled(True)
            self.cb_fg.setEnabled(True)
            #Set defaults
            if not self.cb_id.isChecked() == idstep:
                self.cb_id.toggle()
            if not self.cb_acq1.isChecked() == acq1step:
                self.cb_acq1.toggle()
            if not self.cb_acq2.isChecked() == acq2step:
                self.cb_acq2.toggle()
            if not self.cb_trk.isChecked() == trkstep:
                self.cb_trk.toggle()
            if not self.cb_fg.isChecked() == fgstep:
                self.cb_fg.toggle()
        else:
            self.cb_id.setEnabled(False)
            self.cb_acq1.setEnabled(False)
            self.cb_acq2.setEnabled(False)
            self.cb_trk.setEnabled(False)
            self.cb_fg.setEnabled(False)

    def on_check_segment_guiding(self, is_toggle):
        ''' Checking the segment guiding box. Turn off all other steps in GUI.'''
        if is_toggle:
            if self.cb_convert_im.isChecked():
                self.cb_convert_im.toggle()
            self.cb_convert_im.setEnabled(False)
            if self.cb_file_writer.isChecked():
                self.cb_file_writer.toggle()
            self.cb_file_writer.setEnabled(False)
            if self.cb_star_selection.isChecked():
                self.cb_star_selection.toggle()
            self.cb_star_selection.setEnabled(False)
        else:
            if not self.cb_convert_im.isChecked():
                self.cb_convert_im.toggle()
            self.cb_convert_im.setEnabled(True)
            if not self.cb_file_writer.isChecked():
                self.cb_file_writer.toggle()
            self.cb_file_writer.setEnabled(True)
            if not self.cb_star_selection.isChecked():
                self.cb_star_selection.toggle()
            self.cb_star_selection.setEnabled(True)


    # What to do with specific buttons ------------------------------------------
    def btnstate(self, button):
        ''' Simple button checking code. Only needs to do something specific
        for choosing an FGS image vs NIRCam image.'''
        if button.text() == "FGS Image":
            if button.isChecked():
                self.cb_nircam_det.setEnabled(False)
            else:
                self.cb_nircam_det.setEnabled(True)
        if button.text() == "NIRCam Image":
            if button.isChecked():
                self.cb_nircam_det.setEnabled(True)
            else:
                self.cb_nircam_det.setEnabled(False)

    # Dialog Boxes -------------------------------------------------------------
    def open_filename_dialog(self, title):
        ''' Dialog box for opening a NIRCam of FGS image'''
        options = QFileDialog.Options()
        options |= QFileDialog.DontUseNativeDialog
        fileName, _ = QFileDialog.getOpenFileName(self, "Open {}".format(title),
                                                  "", "All Files (*);;Python Files (*.py)",
                                                  options=options)
        if fileName:
            return fileName

    def open_dirname_dialog(self):
        ''' Dialog box for opening a directory'''
        options = QFileDialog.Options()
        options |= QFileDialog.DontUseNativeDialog
        dirName = str(QFileDialog.getExistingDirectory(self, "Select Directory"))
        if dirName:
            return dirName

    # Helper Functions ---------------------------------------------------------
    def update_filepreview(self):
        if self.textbox_outdir.text() != "" and self.textbox_root.text() != ""\
           and self.cb_guider.currentText() != "-Select-":
            self.textbox_filepreview.setText(os.path.join(self.textbox_outdir.text(),
                                             'out', self.textbox_root.text(),
                                             '{}_G{}_strips.fits'.format(self.textbox_root.text(),
                                                                         self.cb_guider.currentText())))


#===============================================================================
def run_MasterGui(root=None, fgs_counts=None, jmag=11.0, nircam_det=None,
                  nircam=True, global_alignment=False, steps=None, in_file=None,
                  bkgd_stars=False, out_dir=OUT_PATH, convert_im=True,
                  star_selection=True, star_selection_gui=True, file_writer=True,
                  segment_guiding=False):
    # RUN GUI
    app = QtCore.QCoreApplication.instance()  # Use existing instance, if there is one
    if app is None:
        app = QApplication(sys.argv)

    ex = MasterGui(root, fgs_counts, jmag, nircam_det, nircam, global_alignment, steps,
                   in_file, bkgd_stars, out_dir, convert_im, star_selection_gui,
                   file_writer, segment_guiding, app=app)
    # #return ex.settings
    out = app.exec()

    return out


if __name__ == '__main__':
    out = run_MasterGui()
    sys.exit(out)
