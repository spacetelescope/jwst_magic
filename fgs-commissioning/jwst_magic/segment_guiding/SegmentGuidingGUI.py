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
import sys
import os

# Third Party Imports
import numpy as np
from astropy.io import ascii as asc
from PyQt5 import QtCore, uic
from PyQt5.QtWidgets import (QApplication, QDialog, QMessageBox, QTableWidgetItem)
from PyQt5.QtGui import QIcon

# Local Imports
from ..star_selector.SelectStarsGUI import StarSelectorWindow

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
