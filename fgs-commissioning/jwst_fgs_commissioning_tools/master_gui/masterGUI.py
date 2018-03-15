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

import sys
import inspect
from PyQt5.QtWidgets import (QApplication, QWidget, QPushButton, QInputDialog,
                             QLineEdit, QMainWindow, QAction, QMessageBox,
                             QFileDialog, QLabel, QCheckBox, QFrame, QSizePolicy,
                             QHBoxLayout, QSplitter, QTextEdit, QStyleFactory,
                             QSpinBox, QComboBox, QRadioButton, QGridLayout,
                             QGroupBox, qApp)
from PyQt5 import QtGui
from PyQt5.QtGui import QIcon, QFont
from PyQt5.QtCore import pyqtSlot, Qt, QSettings, QTimer
import matplotlib
if matplotlib.get_backend() != 'Qt5Agg':
    matplotlib.use('Qt5Agg')  # Make sure that we are using Qt5


from jwst_fgs_commissioning_tools import run_fgs_commissioning_tool
from jwst_fgs_commissioning_tools.convert_image import counts_to_jmag


class MasterGui(QWidget):

    def __init__(self, root=None, fgs_counts=None, jmag=None,
                 nircam_det=None, nircam=True, global_alignment=False, steps=None,
                 in_file=None, bkgd_stars=False, out_dir=None, convert_im=True,
                 star_selection_gui=True, file_writer=True, segment_guiding=False,
                 parent=None, app=None):
        super().__init__()
        self.title = 'Master GUI - FGS Commissioning Tools'
        self.image_dim = 1000
        self.convert_im = convert_im
        self.star_selection_gui = star_selection_gui
        self.file_writer = file_writer
        self.segment_guiding = segment_guiding
        self.app = app
        self.initUI(root, fgs_counts, jmag, nircam_det, nircam, global_alignment,
                    steps, in_file, bkgd_stars, out_dir, parent)

    def initUI(self, root, fgs_counts, jmag, nircam_det, nircam, global_alignment,
               steps, in_file, bkgd_stars, out_dir, parent):
        # Set default values
        self.jmag_default = 11.0

        # Set the font format for titles
        title_font = QFont()
        title_font.setBold(True)

        # Set window attributes
        QMainWindow.__init__(self, parent)
        self.setWindowTitle(self.title)
        self.main_widget = QWidget(self)
        mainGrid = QGridLayout()  # set grid layout
        self.main_widget.setLayout(mainGrid)
        self.main_widget.setFocus()

        ### Build the MAIN section ---------------------------------------------
        inputGroupBox = QGroupBox('', self)
        inputGrid = QGridLayout()
        inputGroupBox.setLayout(inputGrid)
        # Input Image
        inputGrid.addWidget(QLabel('Input Image', self), 1, 0)
        self.textbox_input_im = QLineEdit(self)
        self.textbox_input_im.setMinimumSize(200, 20)
        inputGrid.addWidget(self.textbox_input_im, 1, 2)

        self.button_input_image = QPushButton('Open', self)
        self.button_input_image.setToolTip('Open the input image')
        self.button_input_image.clicked.connect(self.on_click_input)
        inputGrid.addWidget(self.button_input_image, 1, 5)

        # Guider
        inputGrid.addWidget(QLabel('Guider', self), 2, 0)
        self.textbox_guider = QLineEdit(self)
        inputGrid.addWidget(self.textbox_guider, 2, 2)

        # Root
        inputGrid.addWidget(QLabel('Root', self), 3, 0)
        self.textbox_root = QLineEdit(self)
        inputGrid.addWidget(self.textbox_root, 3, 2)

        # Out directory
        inputGrid.addWidget(QLabel('Out Directory', self), 4, 0)
        self.textbox_outdir = QLineEdit(self)
        inputGrid.addWidget(self.textbox_outdir, 4, 2)
        self.button_input_image = QPushButton('Open', self)
        self.button_input_image.setToolTip('Open the out directory')
        self.button_input_image.clicked.connect(self.on_click_out)
        inputGrid.addWidget(self.button_input_image, 4, 5)

        # Background stars?
        self.cb_bkgd_stars = QCheckBox('Background stars', self)
        self.cb_bkgd_stars.setCheckable(True)
        if bkgd_stars:
            self.cb_bkgd_stars.toggle()
        inputGrid.addWidget(self.cb_bkgd_stars, 5, 0)

        mainGrid.addWidget(inputGroupBox, 1, 0, 5, 3)

        # ### Build the CONVERT IMAGE section ----------------------------------
        convertGroupBox = QGroupBox('Convert Input Image to an FGS Raw Image', self)
        convertGrid = QGridLayout()
        convertGroupBox.setLayout(convertGrid)

        # # If we want to convert image
        self.cb_convert_im = QCheckBox('Convert Image', self)
        self.cb_convert_im.setCheckable(True)
        self.cb_convert_im.setFont(title_font)
        convertGrid.addWidget(self.cb_convert_im, 1, 0)
        if self.convert_im:
            self.cb_convert_im.toggle()

        # Set to be a NIRCam or FGS image
        self.cb_nircam = QCheckBox('NIRCam Image', self)
        self.cb_nircam.setCheckable(True)
        if nircam:
            self.cb_nircam.toggle()
        convertGrid.addWidget(self.cb_nircam, 2, 0)

        self.cb_fgs = QCheckBox('FGS Image', self)
        self.cb_fgs.setCheckable(True)
        if not nircam:
            self.cb_fgs.toggle()
        convertGrid.addWidget(self.cb_fgs, 2, 2)

        # Set NIRCam detector
        convertGrid.addWidget(QLabel('NIRCam Detector:', self), 3, 0)
        self.textbox_nircam_det = QLineEdit(self)
        convertGrid.addWidget(self.textbox_nircam_det, 3, 2)

        #Normalize by FGS counts or JMAG
        # Add checkbox
        self.cb_fgs_counts = QCheckBox('FGS Counts', self)
        convertGrid.addWidget(self.cb_fgs_counts, 4, 0)
        # Add checkbox
        self.cb_jmag = QCheckBox('J mag', self)
        convertGrid.addWidget(self.cb_jmag, 4, 2)
        # Add textbox
        convertGrid.addWidget(QLabel('Value:', self), 4, 3)
        self.textbox_value = QLineEdit(self)
        if jmag:
            self.textbox_value.setText(str(jmag))
        elif not fgs_counts:
            self.textbox_value.setText(str(self.jmag_default))
        convertGrid.addWidget(self.textbox_value, 4, 4)
        # Set toggling for FGS counts and Jmag
        self.cb_fgs_counts.toggled.connect(self.on_check_fgscounts)
        self.cb_jmag.toggled.connect(self.on_check_jmag)
        if fgs_counts:
            self.cb_fgs_counts.toggle()
        else:
            self.cb_jmag.toggle()

        ## Global Alighment flat
        self.cb_ga = QCheckBox('Global Alignment', self)
        self.cb_ga.setCheckable(True)
        if global_alignment:
            self.cb_ga.toggle()
        convertGrid.addWidget(self.cb_ga, 5, 0)

        # Set toggling
        self.cb_convert_im.toggled.connect(self.on_check_convertimage)
        mainGrid.addWidget(convertGroupBox, 7, 0, 3, 3)

        # ### Build the STAR SELECTION section ---------------------------------
        # If we want to use the star slection GUI
        starGroupBox = QGroupBox('Guide and Reference Star Selection', self)
        starGrid = QGridLayout()
        starGroupBox.setLayout(starGrid)

        ## Run the star selection GUI
        # Add checkbox
        self.cb_star_selection = QCheckBox('Star Selection GUI', self)
        self.cb_star_selection.setCheckable(True)
        self.cb_star_selection.toggled.connect(self.on_check_starselect)
        starGrid.addWidget(self.cb_star_selection, 1, 0)

        ## Pass in an in or reg file
        # Add checkbox
        self.cb_infile = QCheckBox('In/Reg file', self)
        self.cb_infile.setCheckable(True)
        self.cb_infile.toggled.connect(self.on_check_infile)
        starGrid.addWidget(self.cb_infile, 2, 0)
        # Add textbox
        self.textbox_infile = QLineEdit(self)
        self.textbox_infile.setMinimumSize(200, 20)
        starGrid.addWidget(self.textbox_infile, 2, 2)
        # Add 'Open' button
        self.button_input_infile = QPushButton('Open', self)
        self.button_input_infile.setToolTip('Open the in_file directory')
        self.button_input_infile.clicked.connect(self.on_click_infile)
        starGrid.addWidget(self.button_input_infile, 2, 5)
        # Check the box that is set to default
        if in_file:
            self.cb_infile.toggle()
        else:
            self.cb_star_selection.toggle()

        mainGrid.addWidget(starGroupBox, 12, 0, 5, 3)


        # ### Build the FSW FILE WRITER section --------------------------------
        # Use the FSW file writer
        filewriterGroupBox = QGroupBox('Output FSW files', self)
        filewriterGrid = QGridLayout()
        filewriterGroupBox.setLayout(filewriterGrid)

        # Write the files out
        self.cb_file_writer = QCheckBox('FSW File Writer', self)
        self.cb_file_writer.setCheckable(True)
        self.cb_file_writer.toggled.connect(self.on_check_filewriter)

        self.cb_file_writer.setFont(title_font)
        filewriterGrid.addWidget(self.cb_file_writer, 1, 0)
        # Which steps to write out
        self.cb_id = QCheckBox('ID', self)
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

        mainGrid.addWidget(filewriterGroupBox, 17, 0, 5, 3)

        # ###  Build the SEGMENT GUIDING section -------------------------------
        # Do segment guiding
        segmentGroupBox = QGroupBox('Do Segment Guiding', self)
        segmentGrid = QGridLayout()
        segmentGroupBox.setLayout(segmentGrid)

        self.cb_segment_guiding = QCheckBox('Segment Guiding', self)
        self.cb_segment_guiding.setCheckable(True)
        self.cb_segment_guiding.toggle()
        self.cb_segment_guiding.setFont(title_font)
        segmentGrid.addWidget(self.cb_segment_guiding, 1, 0)

        mainGrid.addWidget(segmentGroupBox, 22, 0, 5, 3)

        # ###  Save out inputs and run tool ------------------------------------
        button_run = QPushButton("Run", self)
        button_run.resize(button_run.minimumSizeHint())
        button_run.move(self.image_dim*.6, self.image_dim*.75)

        button_done = QPushButton("Done", self)
        button_done.resize(button_done.minimumSizeHint())
        button_done.move(self.image_dim*.75, self.image_dim*.75)

        button_run.clicked.connect(self.run_tool)
        button_done.clicked.connect(self.close_application)
        self.show()


    def close_application(self):
        sys.exit()


    def run_tool(self):
        # Required
        input_image = self.textbox_input_im.text()
        guider = self.textbox_guider.text()
        root = self.textbox_root.text()
        out_dir = self.textbox_outdir.text()
        bkgd_stars = self.cb_bkgd_stars.isChecked()

        #Convert image
        convert_im = self.cb_convert_im.isChecked()
        nircam = self.cb_nircam.isChecked()
        nircam_det = self.textbox_nircam_det.text()
        if self.cb_fgs_counts.isChecked():
            fgs_counts = float(self.textbox_value.text())
        else:
            fgs_counts = None
        if self.cb_jmag.isChecked():
            jmag = float(self.textbox_value.text())
        else:
            jmag = None
        global_alignment = self.cb_ga.isChecked()

        #Star selection
        star_selection = self.cb_star_selection.isChecked()

        #File writer
        file_writer = self.cb_file_writer.isChecked()
        steps = []
        if self.cb_id.isChecked():
            steps.append('ID')
        elif self.cb_acq1.isChecked():
            steps.append('ACQ1')
        elif self.cb_acq2.isChecked():
            steps.append('ACQ2')
        elif self.cb_trk.isChecked():
            steps.append('LOSTRK')
        elif self.cb_fg.isChecked():
            steps.append('FG')

        if not self.cb_infile.isChecked():
            in_file = None
        else:
            in_file = self.textbox_infile.text()


        segment_guiding = self.cb_segment_guiding.isChecked()

        if convert_im or star_selection or file_writer:
            run_fgs_commissioning_tool.run_all(input_image, guider, root,
                                               fgs_counts, jmag, nircam_det,
                                               nircam, global_alignment,
                                               steps, in_file, bkgd_stars,
                                               out_dir, convert_im, star_selection,
                                               self.app)
        if segment_guiding:
            pass


    @pyqtSlot()
    def on_click_input(self):
        ''' Using the Input Image Open button (open file) '''
        filename = self.openFileNameDialog()
        self.textbox_input_im.setText(filename)
        return filename

    @pyqtSlot()
    def on_click_out(self):
        ''' Using the Out Dir Open button (open directory) '''
        dirname = self.openDirNameDialog()
        self.textbox_outdir.setText(dirname)
        return dirname

    @pyqtSlot()
    def on_click_infile(self):
        ''' Using the Infile Open button (open file) '''
        filename = self.openFileNameDialog()
        self.textbox_infile.setText(filename)
        return filename

    def on_check_convertimage(self, is_toggle):
        ''' Checking the convert image box - defaults set here'''
        if is_toggle:
            self.cb_nircam.setEnabled(True)
            self.cb_fgs.setEnabled(True)
            self.textbox_nircam_det.setEnabled(True)
            self.cb_fgs_counts.setEnabled(True)
            self.cb_jmag.setEnabled(True)
            self.textbox_value.setEnabled(True)
            self.cb_ga.setEnabled(True)

        else:
            self.cb_nircam.setEnabled(False)
            self.cb_fgs.setEnabled(False)
            self.textbox_nircam_det.setEnabled(False)
            self.cb_fgs_counts.setEnabled(False)
            self.cb_jmag.setEnabled(False)
            self.textbox_value.setEnabled(False)
            self.cb_ga.setEnabled(False)

        if not self.cb_nircam.isChecked():
            self.cb_nircam.toggle()
        if not self.cb_jmag.isChecked():
            self.cb_jmag.toggle()

    def on_check_nircam(self, is_toggle):
        ''' Checking the NIRCam box in Convert Image - turn off FGS'''
        if is_toggle:
            if self.cb_fgs.isChecked():
                self.cb_fgs.toggle()

    def on_check_fgs(self, is_toggle):
        ''' Checking the FGS box in Convert Image - turn off NIRCam'''
        if is_toggle:
            if self.cb_nircam.isChecked():
                self.cb_nircam.toggle()    

    def on_check_fgscounts(self, is_toggle):
        ''' Checking the FGS counts box in Convert Image - turn off J mag box'''
        if is_toggle:
            if self.cb_jmag.isChecked():
                self.cb_jmag.toggle()
            fgs_counts_default = counts_to_jmag.jmag_to_fgs_counts(self.jmag_default,
                                                                   self.textbox_guider.text())
            self.textbox_value.setText(str(int(fgs_counts_default)))

    def on_check_jmag(self, is_toggle):
        ''' Checking the J mag box in Convert Image - turn off FGS Counts box'''
        if is_toggle:
            if self.cb_fgs_counts.isChecked():
                self.cb_fgs_counts.toggle()
            # Set Value back to default
            self.textbox_value.setText(str(self.jmag_default))

    def on_check_starselect(self, is_toggle):
        ''' Checking the star selection GUI box '''
        if is_toggle:
            self.textbox_infile.setEnabled(False)
            self.button_input_infile.setEnabled(False)
            if self.cb_infile.isChecked():
                self.cb_infile.toggle()
        else:
            self.textbox_infile.setEnabled(True)
            self.button_input_infile.setEnabled(True)

    def on_check_infile(self, is_toggle):
        ''' Checking the star selection infile box '''
        if is_toggle:
            self.textbox_infile.setEnabled(False)
            self.button_input_infile.setEnabled(False)
            if self.cb_star_selection.isChecked():
                self.cb_star_selection.toggle()

    def on_check_filewriter(self, is_toggle):
        ''' Checking the filewriter box - defaults set here'''
        if is_toggle:
            self.cb_id.setEnabled(True)
            self.cb_acq1.setEnabled(True)
            self.cb_acq2.setEnabled(True)
            self.cb_trk.setEnabled(True)
            self.cb_fg.setEnabled(True)
            #Set defaults
            if not self.cb_id.isChecked():
                self.cb_id.toggle()
            if not self.cb_acq1.isChecked():
                self.cb_acq1.toggle()
            if not self.cb_acq2.isChecked():
                self.cb_acq2.toggle()
            if not self.cb_trk.isChecked():
                self.cb_trk.toggle()
        else:
            self.cb_id.setEnabled(False)
            self.cb_acq1.setEnabled(False)
            self.cb_acq2.setEnabled(False)
            self.cb_trk.setEnabled(False)
            self.cb_fg.setEnabled(False)

    def openFileNameDialog(self):
        options = QFileDialog.Options()
        options |= QFileDialog.DontUseNativeDialog
        fileName, _ = QFileDialog.getOpenFileName(self,"Open NIRCam or FGS image",
                                                  "", "All Files (*);;Python Files (*.py)",
                                                  options=options)
        if fileName:
            return fileName

    def openDirNameDialog(self):
        options = QFileDialog.Options()
        options |= QFileDialog.DontUseNativeDialog
        dirName = str(QFileDialog.getExistingDirectory(self,"Select Directory"))
        if dirName:
            return dirName


#===============================================================================
def run_MasterGui(root=None, fgs_counts=None, jmag=None, nircam_det=None,
                  nircam=True, global_alignment=False, steps=None, in_file=None,
                  bkgd_stars=False, out_dir=None, convert_im=True,
                  star_selection_gui=True, file_writer=True, segment_guiding=False):
    app = QApplication(sys.argv)
    ex = MasterGui(root, fgs_counts, jmag, nircam_det, nircam, global_alignment, steps,
              in_file, bkgd_stars, out_dir, convert_im, star_selection_gui,
              file_writer, segment_guiding, app=app)
    sys.exit(app.exec_())
    #return ex.settings

if __name__ == '__main__':
    run_MasterGui()
    # app = QApplication(sys.argv)
    # ex = MasterGui()
    # sys.exit(app.exec_())
