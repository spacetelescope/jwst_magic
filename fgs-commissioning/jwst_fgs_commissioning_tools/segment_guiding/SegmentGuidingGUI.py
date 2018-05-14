"""Interactive star selector GUI for segment guiding.

Authors
-------
    - Lauren Chambers


References
----------
For matplotlib canvas example, see:
    https://matplotlib.org/examples/user_interfaces/embedding_in_qt5.html

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

# Standard Library
from __future__ import unicode_literals
import sys
import os

# Third Party
import numpy as np
from PyQt5 import QtCore, uic
from PyQt5.QtWidgets import (QApplication, QDialog, QMessageBox, QTableWidgetItem)
from PyQt5.QtGui import QIcon

from ..star_selector.SelectStarsGUI import StarSelectorWindow

__location__ = os.path.realpath(os.path.join(os.getcwd(), os.path.dirname(__file__)))


class SegmentGuidingWindow(StarSelectorWindow, QDialog):
    def __init__(self, data, x, y, dist, qApp, in_master_GUI,
                 print_output=False):
        '''Defines attributes; calls initUI() method to set up user interface.'''
        # Initialize attributes
        self.n_orientations = 0
        self.segNum = None
        self.center = None

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

        # Create and load sgement guiding GUI session
        self.setWindowTitle('FGS Commissioning Tools - Segment Guiding Command Generator')
        self.define_SegmentGuidingGUI_connections()
        self.show()

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # GUI CONSTRUCTION
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    def define_SegmentGuidingGUI_connections(self):
        # Main dialog window widgets
        self.pushButton_done.clicked.connect(self.fileQuit)  # Redefine

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
        if self.inds == []:
            no_stars_selected_dialog = QMessageBox()
            no_stars_selected_dialog.setText('No stars selected' + ' ' * 50)
            no_stars_selected_dialog.setInformativeText('Please select at least one star to command.')
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
        # Remove the stars that are currently shown
        self.clear_selected_stars()

        # Determine what are the indices of the stars to load
        selected_row = self.tableWidget_commands.selectedItems()[0].row()
        selected_orientation = self.tableWidget_commands.item(selected_row, 0).text()
        selected_indices = [int(s) for s in selected_orientation.split(', ')]

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
            values = ['', str(ind + 1), str(int(self.x[ind])), str(int(self.y[ind])), str(int(countrate))]
            for i_col, value in enumerate(values):
                if gs and i_col == 0:
                    item = QTableWidgetItem(QIcon(os.path.join(__location__, '../star_selector/gs_icon.png')), '')
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
        remstar_string = 'Loaded orientation: {}'.format(', '.join([str(i + 1) for i in self.inds]))
        self.textEdit_output.setHtml(remstar_string + "<br>" + self.textEdit_output.toHtml())
        if self.print_output:
            print(remstar_string)

    def delete_orientation(self):
        selected_orientation_index = self.tableWidget_commands.selectedItems()[0].row()
        self.tableWidget_commands.removeRow(selected_orientation_index)
        self.n_orientations -= 1

    def move_orientation(self):
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
        '''Use the location of a specific segment as the pointing center
        '''
        # Uncheck "use segment center" box
        self.checkBox_meanCenter.setChecked(False)

        # Remove old center
        if self.center:
            self.canvas.axes.lines.remove(self.center[0])
        self.segNum = None

        if self.comboBox_segmentCenter.currentText() != "-Select Segment-":
            # Replace with new center
            i_seg_center = int(self.comboBox_segmentCenter.currentText()) - 1
            self.center = self.canvas.axes.plot(self.x[i_seg_center],
                                                self.y[i_seg_center], 'x', ms=20,
                                                alpha=0.8, mfc='white',
                                                mec='white', mew=5, lw=0)

            self.segNum = int(self.comboBox_segmentCenter.currentText())

        self.canvas.draw()

    def update_center_mean(self, use_mean_as_center):
        '''Use the mean of the array as the pointing center
        '''
        # Remove old center
        if self.center:
            self.canvas.axes.lines.remove(self.center[0])
        self.segNum = None

        if use_mean_as_center:
            self.comboBox_segmentCenter.setCurrentIndex(0)

            # Calculate center of array
            x_mean = np.average(self.x)
            y_mean = np.average(self.y)

            # Plot mean location of array on canvas
            self.center = self.canvas.axes.plot(x_mean, y_mean, 'x', ms=20, alpha=0.8,
                                                mfc='white', mec='white', mew=5, lw=0)

            self.segNum = 0

        self.canvas.draw()

    def fileQuit(self):
        '''Closes the star selector window'''
        # If the center segment number hasn't been set, don't quit.
        if self.segNum is None:
            no_segNum_selected_dialog = QMessageBox()
            no_segNum_selected_dialog.setText('No center segment number' + ' ' * 50)
            no_segNum_selected_dialog.setInformativeText('The center of override pointing (segNum) has not been defined. Please define before quitting.')
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
            no_stars_selected_dialog.setInformativeText('The tool will not be able to continue. Do you want to quit anyway?')
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


def run_SelectSegmentOverride(data, x, y, dist, print_output=False, masterGUIapp=None):
    '''Calls a PyQt GUI to allow interactive user selection of guide and reference stars.

    Params
    ------
    data (2D array)
        Opened FITS image data (2048 x 2048)
    x (list)
        List of x-coordinates of identified PSFs
    y (list)
        List of y-coordinates of identified PSFs
    dist (int)
        Minimum distance between identified PSFs; maximum distance from a star the user
        can click to select that star
    print_output (boolean)
        Flag enabling output to the terminal

    Returns
    -------
    inds (list)
        List of indices of positions of selected stars
    '''

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
                                  print_output=print_output)

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
