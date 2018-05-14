"""Interactive star selector GUI.

Builds a GUI with PyQt5 that prompts the user to click on a data image
to select guide and/or reference stars, which are then saved as output.
To help the user select desirable stars, the GUI includes a dynamic
matplotlib plot that displays the shape of the signal under the cursor.
The user can also adjust the image boundaries and stretch. Once
multiple stars have been selected, the user can choose to delete a
selected star, or convert a reference star into the guide star.
Built using template found at matplotlib.org (see references).

Authors
-------
    - Lauren Chambers

Use
---
This GUI can be run in the python shell or as a module, as such:
    ::
    from jwst_fgs_commissioning_tools.star_selector import SelectStarsGUI
    inds = SelectStarsGUI.run_SelectStars(data_array, x_list, y_list, dist)

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
import matplotlib
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure
from matplotlib.colors import LogNorm
import matplotlib.pyplot as plt
from PyQt5 import QtCore, uic
from PyQt5.QtWidgets import (QApplication, QDialog, QMessageBox, QSizePolicy,
                             QTableWidgetItem, QWidget, QVBoxLayout)
from PyQt5.QtCore import pyqtSlot, QSize
from PyQt5.QtGui import QIcon

matplotlib.rcParams['font.family'] = 'serif'
matplotlib.rcParams['font.weight'] = 'light'
matplotlib.rcParams['mathtext.bf'] = 'serif:normal'

__location__ = os.path.realpath(os.path.join(os.getcwd(), os.path.dirname(__file__)))

class StarClickerMatplotlibCanvas(FigureCanvas):
    """Creates a matplotlib canvas as a PyQt widget to plot an FGS image.

    Initializes and draws a matplotlib canvas which plots the following:
    the provided FGS image and the locations of the PSFs as identified by
    photutils.find_peaks.

    Attributes
    ----------
    data (2D array)
        Opened FITS image data (should be 2048 x 2048)
    fig (matplotlib Figure)
        Figure to be used in the GUI; one subplot is on self.axes
    axes (matplotlib Axes)
        Axes to be used in the GUI


    Methods
    -------
    compute_initial_figure
        Plots image of data (with colorbar) and locations of identified PSFs
    zoom_to_fit
        Resets axes limits to (0, 2048)
    """

    def __init__(self, parent=None, width=None, height=None, dpi=100, data=None,
                 x=None, y=None, left=None, right=None, bottom=0.05, top=0.95,
                 profile=False):

        self.data = data
        self.x = x
        self.y = y

        if width and height:
            self.fig = Figure(figsize=(width, height), dpi=dpi)
        else:
            self.fig = Figure(dpi=dpi)
        self.axes = self.fig.add_subplot(111)
        self.fig.subplots_adjust(top=top, left=left, right=right,
                                 bottom=bottom)
        self.cbar = []
        self.peaks = []

        self.profile = profile

        FigureCanvas.__init__(self, self.fig)
        self.setParent(parent)
        self.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Expanding)

    def compute_initial_figure(self, fig, data, x, y,
                               xlabel='Raw FGS X (–V3) [pixels]',
                               ylabel='Raw FGS Y (V2) [pixels]'):
        # Plot data
        self.fitsplot = self.axes.imshow(data, cmap='bone', interpolation='nearest',
                                         clim=(0.1, 1e3), norm=LogNorm())
        if self.peaks != []:
            self.peaks.remove()
        self.peaks = self.axes.scatter(x, y, c='r', marker='+')

        # Update colorbar
        if self.cbar != []:
            self.cbar.remove()
        self.cbar = self.fig.colorbar(self.fitsplot, ax=self.axes, fraction=0.046, pad=0.04, norm=LogNorm())

        # Add axis labels
        self.axes.set_xlabel(xlabel)
        self.axes.set_ylabel(ylabel)
        self.draw()

    def init_profile(self):
        profile_min = max(np.min(self.data) / 5, 1e-1)
        self.axes.set_ylim(profile_min, 5 * np.max(self.data))
        self.axes.set_yscale('log')
        self.axes.set_ylabel('Counts')
        self.axes.set_xlabel('X Pixels')

        self.countrate_label = self.axes.text(0.02, 0.9, '3 x 3 countrate:',
                                              transform=self.axes.transAxes)
        self.draw()

    def zoom_to_fit(self):
        self.axes.set_xlim(0, np.shape(self.data)[0])
        self.axes.set_ylim(np.shape(self.data)[1], 0)
        self.draw()

    def zoom_to_crop(self):
        # Calculate initial axis limits to show all sources
        xarray, yarray = [self.y, self.x]
        x_mid = (min(xarray) + max(xarray)) / 2
        y_mid = (min(yarray) + max(yarray)) / 2
        x_range = max(xarray) - min(xarray)
        y_range = max(yarray) - min(yarray)
        ax_range = max(x_range, y_range)  # Choose the larger of the dimensions
        ax_range += 100  # Make sure not to clip off the edge of border PSFS

        self.axes.set_ylim(min(2048, x_mid + ax_range / 2),
                           max(0, x_mid - ax_range / 2))
        self.axes.set_xlim(max(0, y_mid - ax_range / 2),
                           min(2048, y_mid + ax_range / 2))
        self.draw()

    def sizeHint(self):
        if self.profile == True:
            return QSize(200, 250)
        else:
            return QSize(800, 800)

    def hasHeightForWidth(self):
        return self.profile == False


class StarSelectorWindow(QDialog):
    """Interactive PyQt GUI window used to select stars from an FGS image.


    Main Widgets
    ------------
    matplotlib canvas
        Interactive plot of FITS data and locations of PSFs where users select PSFs
    cursor position
        Reports the position of the cursor within the matplotlib canvas in real-time
    output log
        Gives user feedback about selected stars and warnings
    axes limits
        Allows user to change the limits of the FITS data to zoom in or out of the image
    selected stars form
        Shows the positions of user-selected stars in the order they were chosen and
        allows the user to change the selected guide star
    done
        Exits the program


    Attributes
    ----------
    qApp (QApplication)
        Instance of QApplication used to start interactive session
    print_output (boolean)
        Flag to print star selection output to terminal
    title (str)
        Title of application window
    image_dim (int)
        Minimum size of matplotlib canvas, in pixels
    data (2D array)
        Opened FITS image data (should be 2048 x 2048)
    x (list)
        List of x-coordinates of indentified PSFs
    y (list)
        List of y-coordinates of indentified PSFs
    _ind (int)
        Index of the position of the most recently selected star
    inds (list)
        List of indices of the positions of all selected stars
    epsilon (int)
        Furthest distance from a star, in pixels, that a user can click to
        select a star


    Methods
    -------
    initUI
        Sets up user interface with all widgets listed above
    update_axes
        Changes the axes limits of the matplotlib canvas to the current values
        of the axis limit textboxes
    update_textboxes
        Changes the values of the axis limit textboxes to the current axis
        limits of the matplotlib canvas
    update_cursor_position
        Updates the cursor position textbox when the cursor moves within the
        matplotlib axis
    button_press_callback
        Determines if a button press is within the matplotlib axis; if so calls
        get_ind_under_point
    get_ind_under_point
        If user clicks within one epsilon of an identified star, appends that
        star's position to inds; notifies user if no star is nearby or if
        selected star is already selected
    make_setGuideStar
        Updates the order of inds (i.e. the selected guide star) if user changes
        the guide star and updates matplotlib axis accordingly
    fileQuit
        Closes the application
    """
    def __init__(self, data, x, y, dist, qApp, in_master_GUI,
                 print_output=False, in_SGT_GUI=False):
        '''Defines attributes; calls initUI() method to set up user interface.'''
        # Initialize runtime attributes
        self.qApp = qApp
        self.print_output = print_output
        self.in_master_GUI = in_master_GUI

        # Initialize construction attributes
        self.image_dim = 800

        # Initialize data attributes
        self.data = data
        self.x = x
        self.y = y
        self.epsilon = dist
        self._ind = None
        self.inds = []
        self.inds_of_inds = []
        self.n_stars_max = 11
        self.gs_ind = None
        self.currentRow = -1
        self.circles = []

        # Import .ui file
        # (It is imported as a widget, rather than a QDialog window, so that it
        # can also be imported into the SegmentGuidingGUI module)
        self.central_widget = QWidget()
        uic.loadUi(os.path.join(__location__, 'SelectStarsGUI.ui'), self.central_widget)
        self.__dict__.update(self.central_widget.__dict__)

        # Initialize dialog object and add imported UI
        if not in_SGT_GUI:
            QDialog.__init__(self, modal=True)
            self.setLayout(QVBoxLayout())
            self.layout().addWidget(self.central_widget)

        # Create and load GUI session
        self.setWindowTitle('FGS Commissioning Tools - Guide and Reference Star Selector')
        self.init_matplotlib()
        self.define_StarSelectionGUI_connections()
        self.show()

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # GUI CONSTRUCTION
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    def init_matplotlib(self):
        '''Set up the two matplotlib canvases that will preview the
        input image and converted image in the "Image Preview" section.
        '''

        # Connect main matplotlib canvas and add to layout
        self.canvas = StarClickerMatplotlibCanvas(
            parent=self.frame_canvas,  data=self.data, x=self.x, dpi=100,
            y=self.y, left=0.1, right=0.93)
        self.canvas.compute_initial_figure(self.canvas.fig, self.data, self.x,
                                           self.y)
        self.canvas.zoom_to_crop()
        self.update_textboxes()
        self.frame_canvas.layout().insertWidget(0, self.canvas)

        # Connect profile matplotlib canvas and add to layout
        self.canvas_profile = StarClickerMatplotlibCanvas(
            parent=self.frame_profile, data=self.data, dpi=100, profile=True,
            left=0.2, bottom=0.15, right=0.95, top=0.98)
        self.canvas_profile.init_profile()
        self.frame_profile.layout().insertWidget(0, self.canvas_profile)

    def define_StarSelectionGUI_connections(self):
        # Main dialog widgets
        self.pushButton_done.clicked.connect(self.fileQuit)
        self.pushButton_cancel.clicked.connect(self.cancel)

        # Main matplotlib canvas widget
        #   Update cursor position and pixel value under cursor
        self.canvas.mpl_connect('motion_notify_event', self.update_cursor_position)
        self.canvas.mpl_connect('motion_notify_event', self.update_profile)
        #   Star selection!
        self.canvas.mpl_connect('button_press_event', self.button_press_callback)

        # Colorbar widgets
        self.pushButton_colorbar.clicked.connect(self.update_cbar)

        # Axis limits widgets
        self.pushButton_updateAxes.clicked.connect(self.update_axes)
        self.pushButton_zoomFit.clicked.connect(self.canvas.zoom_to_fit)
        self.pushButton_zoomFit.clicked.connect(self.update_textboxes)
        self.pushButton_cropToData.clicked.connect(self.canvas.zoom_to_crop)
        self.pushButton_cropToData.clicked.connect(self.update_textboxes)

        # Star selection widgets
        self.pushButton_makeGuideStar.clicked.connect(self.set_guidestar)
        self.pushButton_deleteStar.clicked.connect(self.remove_star)
        self.pushButton_clearStars.clicked.connect(self.clear_selected_stars)
        self.pushButton_deleteStar.installEventFilter(self)
        self.pushButton_clearStars.installEventFilter(self)
        self.tableWidget_selectedStars.viewport().installEventFilter(self)

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # WIDGET CONNECTIONS
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    @pyqtSlot()
    def eventFilter(self, source, event):
        '''Event filter for "Delete Star" button, "Clear Selection"
        button, and selected stars table.
        '''
        if hasattr(self, 'pushButton_deleteStar') and \
           hasattr(self, 'tableWidget_selectedStars') and \
           hasattr(self, 'pushButton_clearStars'):
            # Parse out the cursor moving in or out of a star deletion button,
            # and updating the matplotlib axis with a red highlighted circle
            # accordingly
            if source == self.pushButton_deleteStar:
                if event.type() == QtCore.QEvent.Enter and len(self.inds) > 0:
                    # Determine index of star corresponding to button
                    star_ind = self.tableWidget_selectedStars.currentRow()

                    # Re-draw selected star as red
                    self.circles[star_ind].set_markeredgecolor('red')
                    self.canvas.draw()
                    return True

                if event.type() == QtCore.QEvent.Leave and len(self.inds) > 0:
                    # Determine index of star corresponding to button
                    star_ind = self.tableWidget_selectedStars.currentRow()

                    # Re-draw selected star as original color
                    if star_ind == self.gs_ind:
                        self.circles[star_ind].set_markeredgecolor('yellow')
                    else:
                        self.circles[star_ind].set_markeredgecolor('darkorange')

                    self.canvas.draw()
                    return True

            # Parse out the cursor moving in or out of the "clear selection"
            # button, and updating the matplotlib axis with red highlighted
            # circles accordingly
            elif source == self.pushButton_clearStars:
                n_stars = len(self.inds)
                if event.type() == QtCore.QEvent.Enter and n_stars > 0:
                    # Re-draw all stars as red
                    for i in range(n_stars):
                        self.circles[i].set_markeredgecolor('red')
                    self.canvas.draw()
                    return True

                if event.type() == QtCore.QEvent.Leave and n_stars > 0:
                    # Re-draw all stars as original color
                    for i in range(n_stars):
                        # Set the guide star back to yellow
                        if i == self.gs_ind:
                            self.circles[i].set_markeredgecolor('yellow')
                        # Set the reference stars back to orange
                        else:
                            self.circles[i].set_markeredgecolor('darkorange')
                    self.canvas.draw()
                    return True

            # Parse any mouse clicks within the table widget, and update the
            # matplotlib axis with a blue highlighted circle to represent the
            # rows that is selected
            elif source == self.tableWidget_selectedStars.viewport() and \
                 event.type() == QtCore.QEvent.MouseButtonPress:
                # Record the previous row that was selected
                previousRow = self.currentRow

                # Determine where the click happened
                if self.tableWidget_selectedStars.itemAt(event.pos()) is None:
                    currentRow = -1
                else:
                    currentRow = self.tableWidget_selectedStars.itemAt(event.pos()).row()

                # Re-draw the newly selected star as blue
                if currentRow >= 0:
                    self.circles[currentRow].set_markeredgecolor('cornflowerblue')

                # Update the table to show the selection
                self.tableWidget_selectedStars.setCurrentCell(currentRow, 0)

                # Put the previously selected star back to normal
                if previousRow == currentRow:
                    pass
                elif previousRow == self.gs_ind:
                    self.circles[previousRow].set_markeredgecolor('yellow')
                elif previousRow >= 0:
                    self.circles[previousRow].set_markeredgecolor('darkorange')

                # Update current row to be the newly clicked one
                self.currentRow = currentRow

                self.canvas.draw()
                return True

        return False

    def update_axes(self):
        '''Changes the axes limits of the matplotlib canvas to the current
        values of the axis limit textboxes
        '''
        self.canvas.axes.set_xlim(float(self.lineEdit_xmin.text()),
                                  float(self.lineEdit_xmax.text()))
        self.canvas.axes.set_ylim(float(self.lineEdit_ymax.text()),
                                  float(self.lineEdit_ymin.text()))
        self.canvas.fig.subplots_adjust(top=.95, left=.07)
        self.canvas.draw()

    def update_cbar(self):
        '''Changes the limits of the colorbar to the current values of the
        colorbar limit textboxes
        '''
        self.canvas.fitsplot.set_clim(float(self.lineEdit_vmin.text()),
                                      float(self.lineEdit_vmax.text()))
        self.canvas.draw()

    def update_textboxes(self):
        '''Changes the values of the axis limit textboxes to the current axis
        limits of the matplotlib canvas'''
        self.lineEdit_xmin.setText(str(self.canvas.axes.get_xlim()[0]))
        self.lineEdit_xmax.setText(str(self.canvas.axes.get_xlim()[1]))
        self.lineEdit_ymax.setText(str(self.canvas.axes.get_ylim()[0]))
        self.lineEdit_ymin.setText(str(self.canvas.axes.get_ylim()[1]))

        self.lineEdit_vmin.setText(str(self.canvas.fitsplot.get_clim()[0]))
        self.lineEdit_vmax.setText(str(self.canvas.fitsplot.get_clim()[1]))

    def update_cursor_position(self, event):
        '''Updates the cursor position textbox when the cursor moves within the
        matplotlib axis'''
        if event.inaxes:
            self.lineEdit_cursorPosition.setText('({:.0f}, {:.0f})'.format(event.xdata,
                                                                           event.ydata))
            self.lineEdit_cursorValue.setText('{:.0f}'.format(self.data[int(event.ydata),
                                                                        int(event.xdata)]))

    def update_profile(self, event):
        '''Updates the profile plot when the cursor moves within the matplotlib axis'''
        if event.inaxes:
            # Determine profile under the cursor
            x = np.arange(max(0, np.floor(event.xdata - self.epsilon)),
                          min(2047, np.ceil(event.xdata + self.epsilon))).astype(int)
            y = self.data[int(event.ydata), x]

            # Make previous lines white so only the current profile is visible
            for line in self.canvas_profile.axes.lines:
                line.set_color('white')

            # Plot the new line and update axis limits
            self.canvas_profile.axes.plot(x, y, c='cornflowerblue')
            self.canvas_profile.axes.set_xlim(int(np.floor(event.xdata - self.epsilon)),
                                       int(np.ceil(event.xdata + self.epsilon)))

            # Update countrate indicator
            countrate = np.sum(self.data[int(event.ydata) - 1:int(event.ydata) + 2,
                                         int(event.xdata) - 1:int(event.xdata) + 2])
            self.canvas_profile.countrate_label.set_text('3 x 3 countrate: {:.0f}'.format(countrate))

            self.canvas_profile.draw()

    def button_press_callback(self, event):
        '''Whenever the mouse is clicked, determines if the click is
        within the matplotlib axis; if so calls get_ind_under_point'''
        if not event.inaxes:
            return
        if event.button != 1:
            return
        self._ind = self.get_ind_under_point(event)

    def set_guidestar(self):
        '''
        Change order of self.inds to reflect the new guide star
        Update list of original index order
        Update GUI plotted circles
        '''
        # Determine index of old guide star
        ind_of_OLD_GS = self.gs_ind

        # Re-draw old guide star as regular star
        self.circles[ind_of_OLD_GS].set_markeredgecolor('darkorange')

        # Determine index of new guide star
        guide_ind = self.tableWidget_selectedStars.currentRow()
        self.gs_ind = guide_ind

        # Re-draw new guide star as guide star
        self.circles[guide_ind].set_markeredgecolor('yellow')

        # Move the guide star icon
        self.tableWidget_selectedStars.setItem(ind_of_OLD_GS, 0, QTableWidgetItem(''))
        self.tableWidget_selectedStars.setItem(guide_ind, 0, QTableWidgetItem(QIcon(os.path.join(__location__, 'gs_icon.png')), ''))

        self.canvas.draw()

    def remove_star(self):
        '''
        Remove star and update inds list
        '''

        # Determine index of star being removed
        star_ind = self.tableWidget_selectedStars.currentRow()

        # If user tries to delete guide star, don't let them
        if star_ind == self.gs_ind:
            redText = "<span style=\" color:#ff0000;\" >"
            redText += 'The guide star cannot be deleted. Please choose \
                        another star or change the guide star.'
            redText += "</span><br>"
            self.textEdit_output.setHtml(redText + self.textEdit_output.toHtml())
            return

        # Un-draw circle
        l_rem = self.circles.pop(star_ind)
        del l_rem

        # Send deletion message to output
        ind_of_ind = int(self.tableWidget_selectedStars.item(star_ind, 1).text()) - 1
        remstar_string = 'Deleted star at: x={:.1f}, y={:.1f}'.format(self.x[ind_of_ind],
                                                                      self.y[ind_of_ind])
        redText = "<span style=\" color:#ff0000;\" >" + \
                  remstar_string + "</span><br>"
        self.textEdit_output.setHtml(redText + self.textEdit_output.toHtml())
        if self.print_output:
            print(remstar_string)

        # Update inds list
        self.inds.pop(star_ind)

        # Update guide star index
        if star_ind < self.gs_ind:
            self.gs_ind -= 1

        # Update table
        self.tableWidget_selectedStars.removeRow(star_ind)

        self.canvas.draw()

    def clear_selected_stars(self):
        # Reset table
        self.tableWidget_selectedStars.setRowCount(0)

        # Add back empty row
        self.tableWidget_selectedStars.insertRow(0)

        # Reset matplotlib canvas
        for circle in self.circles:
            self.canvas.axes.lines.remove(circle)
        self.canvas.draw()

        # Reset index list, guide star index, and circle list
        self.inds = []
        self.gs_ind = None
        self.circles = []

        # Send deletion message to output
        remstar_string = 'Cleared all selected stars.'
        redText = "<span style=\" color:#ff0000;\" >" + \
                  remstar_string + "</span><br>"
        self.textEdit_output.setHtml(redText + self.textEdit_output.toHtml())
        if self.print_output:
            print(remstar_string)

    def fileQuit(self):
        '''Closes the star selector window'''

        self.answer = True
        # If the user didn't choose any stars, ask if they really want to quit.
        if self.inds == []:
            no_stars_selected_dialog = QMessageBox()
            no_stars_selected_dialog.setText('No stars selected' + ' ' * 50)
            no_stars_selected_dialog.setInformativeText('The tool will not be able to continue. Do you want to quit anyway?')
            no_stars_selected_dialog.setStandardButtons(QMessageBox.Yes | QMessageBox.No)
            no_stars_selected_dialog.buttonClicked.connect(self.nostars_dialog)
            no_stars_selected_dialog.exec()

        # If they do, then quit.
        if self.answer:
            # Save out the indices, applying guide star information
            self.reorder_indices()

            # If not being called from the master GUI, exit the whole application
            if not self.in_master_GUI:
                self.qApp.exit(0)  # Works only with self.close() after; same as qApp.quit()

            # Close the star selector dialog window
            self.close()

    def cancel(self):
        '''Closes the star selector window and clears indices'''

        # Clear the indices (i.e. don't save user selections)
        self.inds = []

        # If not being called from the master GUI, exit the whole application
        if not self.in_master_GUI:
            self.qApp.exit(0)  # Works only with self.close() after; same as qApp.quit()

        # Close the star selector dialog window
        self.close()

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # HELPER FUNCTIONS
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    def get_ind_under_point(self, event):
        '''If user clicks within one epsilon of an identified star, appends that
        star's position to inds; notifies user if no star is nearby or if
        selected star is already selected'''

        dist = np.sqrt((self.x - event.xdata)**2 + (self.y - event.ydata)**2)
        indseq = np.nonzero(np.equal(dist, np.amin(dist)))[0]
        ind = indseq[0]
        n_stars = len(self.inds)

        # If the user has already selected the maximum number of stars
        if len(self.inds) == self.n_stars_max:
            if self.print_output:
                print('Maximum number of stars, {}, already selected.'.format(self.n_stars_max))

            redText = "<span style=\" color:#ff0000;\" >"
            redText += 'Maximum number of stars, {}, already selected.'.format(self.n_stars_max)
            redText += "</span><br>"

            self.textEdit_output.setHtml(redText + self.textEdit_output.toHtml())
            return

        # If there is no star within ``dist`` of the mouse click
        elif dist[ind] >= self.epsilon:
            if self.print_output:
                print('No star within {} pixels. No star selected.'.format(self.epsilon))

            redText = "<span style=\" color:#ff0000;\" >"
            redText += 'No star within {} pixels. No star selected.'.format(self.epsilon)
            redText += "</span><br>"

            self.textEdit_output.setHtml(redText + self.textEdit_output.toHtml())

        # If the user clicked on a star that already exists
        elif ind in self.inds:
            if self.print_output:
                print('Star already selected, please choose another star')

            redText = "<span style=\" color:#ff0000;\" >"
            redText += 'Star already selected, please choose another star'
            redText += "</span><br>"

            self.textEdit_output.setHtml(redText + self.textEdit_output.toHtml())

        # Otherwise, actually add the star to the list!
        else:
            gs = len(self.inds) == 0
            if gs:  # First star selected
                c = 'yellow'
                newstar_string = '1 star selected: x={:.1f}, y={:.1f}'.format(self.x[ind],
                                                                              self.y[ind])
                self.gs_ind = 0
            else:  # >= second star selected
                c = 'darkorange'
                newstar_string = '{} stars selected: x={:.1f}, y={:.1f}'.format(len(self.inds) + 1,
                                                                                self.x[ind],
                                                                                self.y[ind])

            # Update log
            if self.print_output:
                print(newstar_string)
            self.textEdit_output.setHtml(newstar_string + '<br>' + self.textEdit_output.toHtml())

            # Add to list
            # Add empty row
            if not gs:
                self.tableWidget_selectedStars.insertRow(self.tableWidget_selectedStars.rowCount())
            n_rows = self.tableWidget_selectedStars.rowCount()
            # Set value of row
            countrate = np.sum(self.data[int(self.y[ind] - 1):int(self.y[ind] + 2),
                                         int(self.x[ind] - 1):int(self.x[ind] + 2)])
            # Populate row with data
            values = ['', str(ind + 1), str(int(self.x[ind])), str(int(self.y[ind])), str(int(countrate))]
            for i_col, value in enumerate(values):
                if gs and i_col == 0:
                    item = QTableWidgetItem(QIcon(os.path.join(__location__, 'gs_icon.png')), '')
                else:
                    item = QTableWidgetItem(value)
                item.setFlags(item.flags() & ~QtCore.Qt.ItemIsEditable)
                self.tableWidget_selectedStars.setItem(n_rows - 1, i_col, item)

            # Plot the new star on the canvas
            self.circles += self.canvas.axes.plot(self.x[ind], self.y[ind], 'o', ms=25,
                                                 mfc='none', mec=c, mew=2, lw=0)
            self.canvas.draw()

            # Add to the list of indices
            self.inds.append(ind)

        return ind

    def reorder_indices(self):
        print(self.gs_ind)
        gs_ID = self.inds[self.gs_ind]
        print(self.inds)
        self.inds.pop(self.gs_ind)
        inds_final = [gs_ID] + self.inds
        self.inds = inds_final
        print(self.inds)

    def nostars_dialog(self, button):
        if 'No' in button.text():
            self.answer = False

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# MAIN FUNCTION
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

def run_SelectStars(data, x, y, dist, print_output=False, masterGUIapp=None):
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

    # RUN GUI
    if masterGUIapp:
        qApp = masterGUIapp
        in_master_GUI = True
    else:
        qApp = QtCore.QCoreApplication.instance()
        if qApp is None:
            qApp = QApplication(sys.argv)
        in_master_GUI = False

    window = StarSelectorWindow(data=data, x=x, y=y, dist=dist, qApp=qApp,
                                in_master_GUI=in_master_GUI,
                                print_output=print_output)

    try:
        plt.get_current_fig_manager().window.raise_()  # Bring window to front
    except AttributeError:
        pass

    if masterGUIapp:
        window.exec_()  # Begin interactive session; pauses until window.exit() is called
    else:
        qApp.exec_()
    inds = window.inds

    return inds
