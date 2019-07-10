"""Interactive star selector GUI for segment guiding.

Builds a GUI with PyQt5 that prompts the user to click on a data image
to select guide and/or reference stars, which are then saved as segment
guiding override commands and written as segment override files in
``segment_guiding.py``. The dynamic matplotlib plot that is users
interact with is imported from ``SelectStarsGUI.py``.

Authors
-------
    - Lauren Chambers

Use
---
    This module can be used as such:
    ::
        from jwst_magic.segment_guiding import SegmentGuidingGUI
        SegmentGuidingGUI.run_segment_override_gui(data, x, y, dist)

    Required arguments:
        ``data`` - image data (2048 x 2048)
        ``x`` - list of x-coordinates of identified PSFs
        ``y`` - list of y-coordinates of identified PSFs
        ``dist`` - minimum distance between identified PSFs; maximum
            distance from a star the user can click to select that star

    Optional arguments:
        ``print_output`` - enable output to the terminal
        ``selected_segs`` - filepath containing locations and count
            rates of segments selected as the guide and reference stars
        ``masterGUIapp`` - qApplication instance of parent GUI

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
JWST MaGIC package, be sure to use the existing instance
of QApplication (access it at QtCore.QCoreApplication.instance()) when
calling the QApplication instance to run a window/dialog/GUI.
"""

# Standard Library Imports
from __future__ import unicode_literals
import glob
import logging
import re
import sys
import os

# Third Party Imports
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.io import ascii as asc
import numpy as np
from PyQt5 import QtCore, uic
from PyQt5.QtWidgets import (QApplication, QDialog, QMessageBox, QTableWidgetItem, QWidget)
from PyQt5.QtGui import QIcon

# Local Imports
from ..star_selector.SelectStarsGUI import StarSelectorWindow
from ..utils.coordinate_transforms import nrca3pixel_offset_to_v2v3_offset

# Paths
__location__ = os.path.realpath(os.path.join(os.getcwd(), os.path.dirname(__file__)))


class SegmentGuidingWindow(StarSelectorWindow, QDialog):
    def __init__(self, data, x, y, dist, qApp, in_master_GUI,
                 print_output=False, selected_segs=None):
        """Initializes class; calls initUI() method to set up user interface.

        Parameters
        ----------
        data : 2-D numpy array
            Image data (2048 x 2048)
        x : list
            List of x-coordinates of identified PSFs
        y : list
            List of y-coordinates of identified PSFs
        dist : int
            Minimum distance between identified PSFs; maximum distance from a star the user
            can click to select that star
        qApp : qApplication
            qApplication instance of parent GUI
        in_master_GUI : bool
            Denotes if the GUI is being launched as a dialog box from
            the master GUI
        print_output : bool, optional
            Flag enabling output to the terminal
        selected_segs : str, optional
            Filepath containing locations and count rates of segments
            selected as the guide and reference stars
        """
        # Initialize attributes
        self.n_orientations = 0
        self.segNum = None
        self.center = None
        self.selected_segs = selected_segs

        # Initialize dialog object
        QDialog.__init__(self, modal=True)

        # Import .ui file
        uic.loadUi(os.path.join(__location__, 'SegmentGuidingGUI.ui'), self)

        # Initialize the StarSelector Window that is inherited
        StarSelectorWindow.__init__(self, data, x, y, dist, qApp, in_master_GUI,
                                    print_output=print_output, in_SGT_GUI=True)

        # Add the star selector GUI to the segment guiding window
        starSelectorCentralWidget = self.central_widget
        self.frame_selectStarsGUI.layout().addWidget(starSelectorCentralWidget)

        # Create and load segment guiding GUI session
        self.show()
        self.setWindowTitle('JWST MaGIC - Segment Guiding Command Generator')
        self.adjust_screen_size_segment_guiding()
        self.define_SegmentGuidingGUI_connections()
        self.show()

        # Add the number of stars/segments to the pointing center dropdown menu
        for i in range(len(x)):
            self.comboBox_segmentCenter.addItem('{}'.format(i+1))

        # Load stars from guiding_selections*.txt, if it exists
        if self.selected_segs is not None and os.path.exists(self.selected_segs):
            self.load_orientation()
        self.selected_segs = None

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # GUI CONSTRUCTION
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    def adjust_screen_size_segment_guiding(self):
        """ Adjust the GUI sizes for laptop screens
        """
        # Determine screen size
        screen_size = self.qApp.desktop().screenGeometry()
        width, height = screen_size.width(), screen_size.height()

        # Adjust the scroll window size
        minimum_width = 1110 + 400
        self.scrollArea_segmentGuiding.setMinimumWidth(minimum_width)

        if width - 200 < self.scrollArea_segmentGuiding.minimumWidth():
            # Window is too wide
            self.scrollArea_segmentGuiding.setMinimumWidth(width - 200)
        if height - 200 < self.scrollArea_segmentGuiding.minimumHeight():
            # Window is too tall
            self.scrollArea_segmentGuiding.setMinimumHeight(height - 200)

        # Adjust the widget sizes
        if height < 1100:
            self.tableWidget_commands.setMinimumSize(220, 250)
            self.groupBox_commands.setMinimumSize(280, 380 - 50)
            self.frame_segmentGuiding.layout().setSpacing(20)
        else:
            self.scrollArea_segmentGuiding.setMinimumSize(minimum_width + 200, 950 + 200)

    def define_SegmentGuidingGUI_connections(self):
        """Connect widgets' signals to the appropriate methods.
        """
        # Main dialog window widgets
        self.pushButton_done.clicked.connect(self.quit_segment_guiding)  # Redefine

        # General segment guiding widgets
        self.pushButton_save.clicked.connect(self.save_orientation_to_list)

        # Command widgets
        self.pushButton_loadCommand.clicked.connect(self.load_orientation)
        self.pushButton_deleteCommand.clicked.connect(self.delete_orientation)
        self.toolButton_moveUp.clicked.connect(self.move_orientation)
        self.toolButton_moveDown.clicked.connect(self.move_orientation)

        # Override center widgets
        self.checkBox_meanCenter.toggled.connect(self.update_center_mean)
        self.comboBox_segmentCenter.activated.connect(self.update_center_seg)

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # WIDGET CONNECTIONS
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    def save_orientation_to_list(self):
        """Save the currently selected segment orientation to the list
        of override commands.
        """
        if not self.inds:
            no_stars_selected_dialog = QMessageBox()
            no_stars_selected_dialog.setText('No stars selected' + ' ' * 50)
            no_stars_selected_dialog.setInformativeText(
                'Please select at least one star to command.'
            )
            no_stars_selected_dialog.setStandardButtons(QMessageBox.Ok)
            no_stars_selected_dialog.exec()
            return

        if not self.tableWidget_commands.isEnabled():
            self.tableWidget_commands.setEnabled(True)
            self.tableWidget_commands.removeRow(0)

        n_commands = self.tableWidget_commands.rowCount()

        self.reorder_indices()

        orientation_summary = ', '.join([str(i + 1) for i in self.inds])
        item = QTableWidgetItem(orientation_summary)
        self.tableWidget_commands.insertRow(n_commands)
        self.tableWidget_commands.setItem(n_commands, 0, item)
        self.clear_selected_stars()
        self.n_orientations += 1

    def load_orientation(self):
        """Plot the segments in the currently selected override command
        on the matplotlib canvas

        Parameters
        ----------
        selected_segs : str, optional
            Filepath to guiding_selections*.txt file with list of locations and
            countrates for the selected segments (guide and reference stars).
        """
        # Determine what are the indices of the stars to load

        # Read them from a guiding_selections*.txt
        if self.selected_segs is not None:
            guiding_selections_locations = asc.read(self.selected_segs)
            x_reg = guiding_selections_locations['x']
            y_reg = guiding_selections_locations['y']
            selected_indices = []
            for x, y in zip(x_reg, y_reg):
                for i in range(len(self.x)):
                    if self.x[i] == x and self.y[i] == y:
                        selected_indices.append(i + 1)
            guiding_selections_file_message = ' from guiding_selections*.txt'

            # If the guiding_selections*.txt doesn't match the all_found_psfs*.txt locations,
            # don't try to load anything.
            if not selected_indices:
                return

        # Read them from a table
        elif self.selected_segs is None:
            # Remove the stars that are currently shown
            self.clear_selected_stars()

            # Read the table
            selected_row = self.tableWidget_commands.selectedItems()[0].row()
            selected_orientation = self.tableWidget_commands.item(selected_row, 0).text()
            selected_indices = [int(s) for s in selected_orientation.split(', ')]
            guiding_selections_file_message = ''

        # Update the table and canvas
        for i_row, ind in enumerate(selected_indices):
            ind = ind - 1
            # Add rows to the table
            gs = i_row == 0
            if gs:
                self.gs_ind = 0
            else:
                self.tableWidget_selectedStars.insertRow(self.tableWidget_selectedStars.rowCount())
            # Populate row with data
            countrate = np.sum(self.data[int(self.y[ind] - 1):int(self.y[ind] + 2),
                                         int(self.x[ind] - 1):int(self.x[ind] + 2)])
            values = ['', str(ind + 1), str(int(self.x[ind])),
                      str(int(self.y[ind])), str(int(countrate))]
            for i_col, value in enumerate(values):
                if gs and i_col == 0:
                    icon = QIcon(os.path.join(__location__, '../star_selector/gs_icon.png'))
                    item = QTableWidgetItem(icon, '')
                else:
                    item = QTableWidgetItem(value)
                item.setFlags(item.flags() & ~QtCore.Qt.ItemIsEditable)
                self.tableWidget_selectedStars.setItem(i_row, i_col, item)

            # Load matplotlib canvas
            if gs:
                c = 'yellow'
            else:
                c = 'darkorange'
            self.circles += self.canvas.axes.plot(self.x[ind], self.y[ind], 'o',
                                                  ms=25, mfc='none', mec=c,
                                                  mew=2, lw=0)

        self.canvas.draw()

        # Load index list
        self.inds = [i - 1 for i in selected_indices]

        # Send message to output
        remstar_string = 'Loaded orientation{}: {}'.\
            format(guiding_selections_file_message, ', '.join([str(i) for i in selected_indices]))
        self.textEdit_output.setHtml(remstar_string + "<br>" + self.textEdit_output.toHtml())
        if self.print_output:
            print(remstar_string)

    def delete_orientation(self):
        """Delete the selected override command.
        """
        selected_orientation_index = self.tableWidget_commands.selectedItems()[0].row()
        self.tableWidget_commands.removeRow(selected_orientation_index)
        self.n_orientations -= 1

    def move_orientation(self):
        """Move the selected override command up or down within the
        list of commands.
        """
        # Remove the item
        selected_orientation_index = self.tableWidget_commands.selectedItems()[0].row()
        item = self.tableWidget_commands.takeItem(selected_orientation_index, 0)
        self.tableWidget_commands.removeRow(selected_orientation_index)

        # Move it up or down (but don't move if at the end of the list)
        if self.sender() == self.toolButton_moveUp:
            new_index = max(0, selected_orientation_index - 1)
        elif self.sender() == self.toolButton_moveDown:
            new_index = min(self.tableWidget_commands.rowCount(), selected_orientation_index + 1)
        self.tableWidget_commands.insertRow(new_index)
        self.tableWidget_commands.setItem(new_index, 0, item)

        # Maintain which item is selected
        self.tableWidget_commands.setCurrentItem(item)

    def update_center_seg(self):
        """Use the location of a specific segment as the pointing center.
        """
        # Uncheck "use segment center" box
        self.checkBox_meanCenter.setChecked(False)

        # Remove old center
        if self.center:
            self.canvas.axes.lines.remove(self.center[0])
            self.center = None
        self.segNum = None

        if self.comboBox_segmentCenter.currentText() != "-Select Segment-":
            # Replace with new center
            i_seg_center = int(self.comboBox_segmentCenter.currentText()) - 1
            self.center = self.canvas.axes.plot(self.x[i_seg_center],
                                                self.y[i_seg_center], 'x', ms=20,
                                                alpha=0.8, mfc='red',
                                                mec='red', mew=5, lw=0)

            self.segNum = int(self.comboBox_segmentCenter.currentText())

        self.canvas.draw()

    def update_center_mean(self, use_mean_as_center):
        """Use the mean of the array as the pointing center.
        """
        # Remove old center
        if self.center:
            self.canvas.axes.lines.remove(self.center[0])
            self.center = None
        self.segNum = None

        if use_mean_as_center:
            self.comboBox_segmentCenter.setCurrentIndex(0)

            # Calculate center of array
            x_mean = np.average(self.x)
            y_mean = np.average(self.y)

            # Plot mean location of array on canvas
            self.center = self.canvas.axes.plot(x_mean, y_mean, 'x', ms=20, alpha=0.8,
                                                mfc='red', mec='red', mew=5, lw=0)

            self.segNum = 0

        self.canvas.draw()

    def quit_segment_guiding(self):
        """Ensures that all necessary values have been defined and
        closes the segment guiding window.
        """
        # If the center segment number hasn't been set, don't quit.
        if self.segNum is None:
            no_segNum_selected_dialog = QMessageBox()
            no_segNum_selected_dialog.setText('No center segment number' + ' ' * 50)
            no_segNum_selected_dialog.setInformativeText(
                'The center of override pointing (segNum) has not been defined.'
                ' Please define before quitting.'
            )
            no_segNum_selected_dialog.setStandardButtons(QMessageBox.Ok)
            no_segNum_selected_dialog.exec()
            return

        self.answer = True
        # If the user selected stars but didn't explicitly save them as a
        # command, do it for them
        if self.inds != [] and self.n_orientations == 0:
            self.save_orientation_to_list()

        # If the user really didn't choose any stars, ask if they really want to quit.
        elif self.inds == [] and self.n_orientations == 0:
            no_stars_selected_dialog = QMessageBox()
            no_stars_selected_dialog.setText('No stars selected' + ' ' * 50)
            no_stars_selected_dialog.setInformativeText(
                'The tool will not be able to continue. Do you want to quit anyway?'
            )
            no_stars_selected_dialog.setStandardButtons(QMessageBox.Yes | QMessageBox.No)
            no_stars_selected_dialog.exec()
        # If they do, then quit.
        if self.answer:
            # If not being called from the master GUI, exit the whole application
            if not self.in_master_GUI:
                self.qApp.exit(0)  # Works only with self.close() after; same as qApp.quit()

            # Close the star selector dialog window
            self.close()


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
    observation_num : optional, int
        Observation number
    visit_num : optional, int
        Visit number

    Returns
    -------
    tup
        Tuple containing the following arguments: (guide_star_params_dict,
        program_id, observation_num, visit_num, threshold_factor,
        countrate_factor)
    """
    def __init__(self, override_type, guider, program_id, observation_num, visit_num, log=None):
        # Initialize attributes
        self.override_type = override_type
        self.guider = guider
        self.program_id = program_id
        self.observation_num = observation_num
        self.visit_num = visit_num

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
                v2_offset, v3_offset = nrca3pixel_offset_to_v2v3_offset(x_offset,
                                                                        y_offset)
                self.log.info(
                    'Segment Guiding: Applying boresight offset of {}, {} arcsec (Converted from {}, {} pixels)'.
                        format(v2_offset, v3_offset, x_offset, y_offset)
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


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# MAIN FUNCTION
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


def run_segment_override_gui(data, x, y, dist, print_output=False,
                             selected_segs=None, masterGUIapp=None):
    """Constructs a PyQt GUI to allow interactive user selection of
    guide and reference stars.

    Parameters
    ----------
    data : 2-D numpy array
        Image data (2048 x 2048)
    x : list
        List of x-coordinates of identified PSFs
    y : list
        List of y-coordinates of identified PSFs
    dist : int
        Minimum distance between identified PSFs; maximum distance from a star the user
        can click to select that star
    print_output : bool, optional
        Flag enabling output to the terminal
    selected_segs : str, optional
        Filepath containing locations and count rates of segments
        selected as the guide and reference stars
    masterGUIapp : qApplication, optional
        qApplication instance of parent GUI

    Returns
    -------
    inds : list
        List of indices of positions of selected stars
    segNum : int
        Number of the segment used as the center of the pointing (the
        mean of the array was used if segNum is 0)
    """

    # Alter data to accurately display null or negative values in log scale plot
    data[data <= 0] = 1

    # RUN GUI
    if masterGUIapp:
        qApp = masterGUIapp
        in_master_GUI = True
    else:
        qApp = QtCore.QCoreApplication.instance()
        if qApp is None:
            qApp = QApplication(sys.argv)
        in_master_GUI = False

    window = SegmentGuidingWindow(data=data, x=x, y=y, dist=dist, qApp=qApp,
                                  in_master_GUI=in_master_GUI,
                                  print_output=print_output,
                                  selected_segs=selected_segs)

    # Begin interactive session; pauses until window.exit() is called
    if masterGUIapp:
        window.exec_()
    else:
        qApp.exec_()

    # Save indices of selected stars to pass
    inds = []
    for i in range(window.n_orientations):
        orientation = window.tableWidget_commands.item(i, 0).text()
        selected_indices = [int(s) for s in orientation.split(', ')]
        inds.append(selected_indices)

    # Save index of center segment (pointing)
    segNum = window.segNum

    return inds, segNum
