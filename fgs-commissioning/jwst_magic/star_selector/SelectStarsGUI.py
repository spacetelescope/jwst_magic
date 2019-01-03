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
    from jwst_magic.star_selector import SelectStarsGUI
    inds = SelectStarsGUI.run_SelectStars(data, x, y, dist)

    Required arguments:
        ``data`` - image data (2048 x 2048)
        ``x`` - list of x-coordinates of identified PSFs
        ``y`` - list of y-coordinates of identified PSFs
        ``dist`` - minimum distance between identified PSFs; maximum
            distance from a star the user can click to select that star

    Optional arguments:
        ``print_output`` - enable output to the terminal
        ``masterGUIapp`` - qApplication instance of parent GUI

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
import matplotlib
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure
from matplotlib.colors import LogNorm
import matplotlib.pyplot as plt
from PyQt5 import QtCore, uic
from PyQt5.QtWidgets import (QApplication, QDialog, QMessageBox, QSizePolicy,
                             QTableWidgetItem, QWidget, QVBoxLayout, QScrollArea,
                             QFrame)
from PyQt5.QtCore import pyqtSlot, QSize
from PyQt5.QtGui import QIcon

# Local Imports
from ..match_to_wss import MatchToWss

# Adjust matplotlib parameters
matplotlib.rcParams['font.family'] = 'serif'
matplotlib.rcParams['font.weight'] = 'light'
matplotlib.rcParams['mathtext.bf'] = 'serif:normal'

# Paths
__location__ = os.path.realpath(os.path.join(os.getcwd(), os.path.dirname(__file__)))


class StarClickerMatplotlibCanvas(FigureCanvas):

    def __init__(self, parent=None, width=None, height=None, dpi=100, data=None,
                 x=None, y=None, left=None, right=None, bottom=0.05, top=0.95,
                 profile=False):
        """Creates a matplotlib canvas as a PyQt widget to plot an FGS image.

        Initializes and draws a matplotlib canvas which plots the following
        the provided FGS image and the locations of the PSFs as identified by
        photutils.find_peaks.

        Parameters
        ----------
        parent : QWidget, optional
            Parent widget of the matplotlib canvas
        width : int, optional
            Width of the matplotlib axis
        height : int, optional
            Height of the matplotlib axis
        dpi : int, optional
            DPI of the matplotlib figure
        data : 2-D numpy array, optional
            Image data
        x : list
            List of x-coordinates of all identified PSFs
        y : list
            List of y-coordinates of all identified PSFs
        left : float, optional
            Left boundary of matplotlib figure
        right : float, optional
            Right boundary of matplotlib figure
        bottom : float, optional
            Bottom boundary of matplotlib figure
        top : float, optional
            Top boundary of matplotlib figure
        profile : bool, optional
            Denotes if the canvas will show profile data (True) or image
            data (False)
        """
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
        """Plot the figure axes and colorbar with given image data and
        coordinates.

        Parameters
        ----------
        fig : matplotlib Figure
            Figure to plot to
        data : 2-D numpy array
            Image data
        x : list
            List of x-coordinates of all identified PSFs
        y : list
            List of y-coordinates of all identified PSFs
        xlabel : str, optional
            X axis label
        ylabel : str, optional
            Y axis label
        """
        # Plot data
        self.fitsplot = self.axes.imshow(data, cmap='bone', interpolation='nearest',
                                         clim=(0.1, 1e3), norm=LogNorm())
        if self.peaks:
            self.peaks.remove()
        self.peaks = self.axes.scatter(x, y, c='r', marker='+')

        # Update colorbar
        if self.cbar:
            self.cbar.remove()
        self.cbar = self.fig.colorbar(self.fitsplot, ax=self.axes,
                                      fraction=0.046, pad=0.04, norm=LogNorm())

        # Add axis labels
        self.axes.set_xlabel(xlabel)
        self.axes.set_ylabel(ylabel)
        self.draw()

    def init_profile(self):
        """Initialize the PSF profile figure.
        """
        profile_min = max(np.min(self.data) / 5, 1e-1)
        self.axes.set_ylim(profile_min, 5 * np.max(self.data))
        self.axes.set_yscale('log')
        self.axes.set_ylabel('Counts')
        self.axes.set_xlabel('X Pixels')

        self.countrate_label = self.axes.text(0.02, 0.9, '3 x 3 countrate:',
                                              transform=self.axes.transAxes)
        self.draw()

    def zoom_to_fit(self):
        """Update the axes limits to fit the entire image.
        """
        self.axes.set_xlim(0, np.shape(self.data)[0])
        self.axes.set_ylim(np.shape(self.data)[1], 0)
        self.draw()

    def zoom_to_crop(self):
        """Update the axes limits to crop to just the source locations
        """
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
        """Reimplementation of the sizeHint method for the matplotlib
        canvas widget. Determines what the size hint should be based on
        if the canvas is plotting a profile or an image.

        Returns
        -------
        QSize
            Size of widget
        """
        if self.profile is True:
            return QSize(200, 250)
        else:
            return QSize(800, 800)

    def hasHeightForWidth(self):
        """Reimplementation of the hasHeightForWidth method for the
        matplotlib canvas widget. Ensures the widget is square if the
        canvas is plotting an image.

        Returns
        -------
        QSize
            Size of widget
        """
        return not self.profile


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
        List of x-coordinates of identified PSFs
    y (list)
        List of y-coordinates of identified PSFs
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
    quit
        Closes the application
    """
    def __init__(self, data, x, y, dist, qApp, in_master_GUI,
                 print_output=False, in_SGT_GUI=False):
        """Initializes class; sets up user interface.

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
        in_SGT_GUI : bool, optional
            Denotes if this instance of the StarSelectorWindow was
            created within ``SegmentGuidingGUI.py``
        """
        # Initialize runtime attributes
        self.qApp = qApp
        self.print_output = print_output
        self.in_master_GUI = in_master_GUI
        self.in_SGT_GUI = in_SGT_GUI

        # Initialize construction attributes
        self.image_dim = 800

        # Initialize data attributes
        self.data = data
        self.x = x
        self.y = y
        self.epsilon = dist
        self._ind = None
        self.inds = [] # FIXME! What is this?
        self.n_stars_max = 11
        self.gs_ind = None
        self.current_row = -1
        self.circles = []

        # Get WSS seg numbers
        #FIXME!
        coords = list(map(list, zip(x, y)))
        wss_obj = MatchToWss(coords)
        self.wss_dict = wss_obj.dictionary
        # FIXME
        # I think I might want to set these here so that you don't need to select segs
        # in order to have this attached to the class
        # self.inds will match order of x, y unless otherwise states??
        # Or inds are selected stars and self.segnums are the associated segment numbers
        self.segnums = []
        for coordx, coordy in zip(self.x, self.y):
            self.segnums.append(next((name for name, seg in self.wss_dict.items() if seg['coords'] == (coordx, coordy)), None))
        print(self.segnums)
        # Import .ui file
        # (It is imported as a widget, rather than a QDialog window, so that it
        # can also be imported into the SegmentGuidingGUI module)
        self.central_widget = QWidget()
        uic.loadUi(os.path.join(__location__, 'SelectStarsGUI.ui'), self.central_widget)
        self.__dict__.update(self.central_widget.__dict__)

        # Initialize dialog object and add imported UI
        if not in_SGT_GUI:
            QDialog.__init__(self, modal=True)
            # Set dialog window layout
            self.setLayout(QVBoxLayout())

            # Create scroll area and add it to the dialog window
            self.scrollArea_selectStars = QScrollArea()
            self.layout().addWidget(self.scrollArea_selectStars)

            # Modify scroll area
            self.scrollArea_selectStars.setMinimumSize(1110, 950)
            self.scrollArea_selectStars.setFrameShape(QFrame.NoFrame)
            self.scrollArea_selectStars.setWidgetResizable(True)

            # Add central widget to scroll area
            self.scrollArea_selectStars.setWidget(self.central_widget)

        # Create and load GUI session
        self.setWindowTitle('JWST MaGIC - Guide and Reference Star Selector')
        canvas_left, profile_bottom = self.adjust_screen_size_select_stars()
        self.init_matplotlib(canvas_left, profile_bottom)
        self.define_StarSelectionGUI_connections()
        self.show()

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # GUI CONSTRUCTION
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    def adjust_screen_size_select_stars(self):
        """ Adjust the GUI sizes for laptop screens
        """

        # Determine screen size
        screen_size = self.qApp.desktop().screenGeometry()
        width, height = screen_size.width(), screen_size.height()

        # Adjust the matplotlib frame sizes
        if height < 1100:
            canvas_left = 0.13
            profile_bottom = 0.2
        else:
            self.frame_canvas.setMinimumSize(800, 0)
            self.tableWidget_selectedStars.setMinimumSize(420, 300)
            self.groupBox_selectedStars.setMinimumSize(440, 380)
            self.frame_profile.setMinimumSize(0, 300)
            if not self.in_SGT_GUI:
                self.scrollArea_selectStars.setMinimumSize(1110 + 200, 950 + 200)
            canvas_left = 0.1
            profile_bottom = 0.15

        # If the window is standalone, adjust the scroll window size
        if not self.in_SGT_GUI:
            # Window is too wide
            if width - 200 < self.scrollArea_selectStars.minimumWidth():
                # print(
                #     'Star selection window is too wide (screen size: {}); resizing to {}'.
                #     format(width, width - 200)
                # )
                self.scrollArea_selectStars.setMinimumWidth(width - 200)
            # Window is too tall
            if height - 200 < self.scrollArea_selectStars.minimumHeight():
                # print(
                #     'Star selection window is too tall (screen size: {}); resizing to {}'.
                #     format(height, height - 200)
                # )
                self.scrollArea_selectStars.setMinimumHeight(height - 200)

        return canvas_left, profile_bottom

    def init_matplotlib(self, canvas_left, profile_bottom):
        """Set up the two matplotlib canvases that will preview the
        input image and converted image in the "Image Preview" section.

        Parameters
        ----------
        canvas_left : float
            Left boundary of matplotlib figure
        profile_bottom : float
            Bottom boundary of matplotlib figure
        """

        # Connect main matplotlib canvas and add to layout
        self.canvas = StarClickerMatplotlibCanvas(
            parent=self.frame_canvas,  data=self.data, x=self.x, dpi=100,
            y=self.y, left=canvas_left, right=0.93)
        self.canvas.compute_initial_figure(self.canvas.fig, self.data, self.x,
                                           self.y)
        self.canvas.zoom_to_crop()
        self.update_textboxes()
        self.frame_canvas.layout().insertWidget(0, self.canvas)

        # Connect profile matplotlib canvas and add to layout
        self.canvas_profile = StarClickerMatplotlibCanvas(
            parent=self.frame_profile, data=self.data, dpi=100, profile=True,
            left=0.2, bottom=profile_bottom, right=0.95, top=0.98)
        self.canvas_profile.init_profile()
        self.frame_profile.layout().insertWidget(0, self.canvas_profile)

    def define_StarSelectionGUI_connections(self):
        """Connect widgets' signals to the appropriate methods.
        """
        # If not in SGT, define the main GUI widgets
        if not self.in_SGT_GUI:
            self.pushButton_done.clicked.connect(self.quit_star_selection)
        self.pushButton_cancel.clicked.connect(self.cancel)

        # Main matplotlib canvas widget
        #   Update cursor position and pixel value under cursor
        self.canvas.mpl_connect('motion_notify_event', self.update_cursor_position)
        self.canvas.mpl_connect('motion_notify_event', self.update_profile)
        self.canvas.mpl_connect('motion_notify_event', self.show_segment_id)
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
        """Event filter for "Delete Star" button, "Clear Selection"
        button, and selected stars table.

        Parameters
        ----------
        source : QWidget
            Widget that signalled to the event filter
        event : QEvent
            The event that occurred within the source widget

        Returns
        -------
        bool
            Denotes an action occurred within the event filter
        """
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
                    if star_ind != -1:
                        self.circles[star_ind].set_markeredgecolor('red')
                        self.canvas.draw()
                    return True

                if event.type() == QtCore.QEvent.Leave and len(self.inds) > 0:
                    # Determine index of star corresponding to button
                    star_ind = self.tableWidget_selectedStars.currentRow()

                    # Re-draw selected star as selected color
                    if star_ind != -1:
                        self.circles[star_ind].set_markeredgecolor('cornflowerblue')

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
            elif (source == self.tableWidget_selectedStars.viewport()) and \
                 (event.type() == QtCore.QEvent.MouseButtonPress):
                # Record the previous row that was selected
                previous_row = self.current_row

                # Determine where the click happened
                if self.tableWidget_selectedStars.itemAt(event.pos()) is None:
                    current_row = -1
                else:
                    current_row = self.tableWidget_selectedStars.itemAt(event.pos()).row()

                # Re-draw the newly selected star as blue
                if current_row >= 0:
                    self.circles[current_row].set_markeredgecolor('cornflowerblue')

                # Update the table to show the selection
                self.tableWidget_selectedStars.setCurrentCell(current_row, 0)

                # Put the previously selected star back to normal
                if previous_row == current_row:
                    pass
                elif previous_row == self.gs_ind:
                    self.circles[previous_row].set_markeredgecolor('yellow')
                elif previous_row >= 0 and len(self.inds) > 0:
                    # Except for the case where all the stars have been deleted
                    self.circles[previous_row].set_markeredgecolor('darkorange')

                # Update current row to be the newly clicked one
                self.current_row = current_row

                self.canvas.draw()
                return True

        return False

    def update_axes(self):
        """Changes the axes limits of the matplotlib canvas to the current
        values of the axis limit textboxes.
        """
        self.canvas.axes.set_xlim(float(self.lineEdit_xmin.text()),
                                  float(self.lineEdit_xmax.text()))
        self.canvas.axes.set_ylim(float(self.lineEdit_ymax.text()),
                                  float(self.lineEdit_ymin.text()))
        self.canvas.fig.subplots_adjust(top=.95, left=.07)
        self.canvas.draw()

    def update_cbar(self):
        """Changes the limits of the colorbar to the current values of the
        colorbar limit textboxes.
        """
        self.canvas.fitsplot.set_clim(float(self.lineEdit_vmin.text()),
                                      float(self.lineEdit_vmax.text()))
        self.canvas.draw()

    def update_textboxes(self):
        """Changes the values of the axis limit textboxes to the current axis
        limits of the matplotlib canvas.
        """
        self.lineEdit_xmin.setText(str(self.canvas.axes.get_xlim()[0]))
        self.lineEdit_xmax.setText(str(self.canvas.axes.get_xlim()[1]))
        self.lineEdit_ymax.setText(str(self.canvas.axes.get_ylim()[0]))
        self.lineEdit_ymin.setText(str(self.canvas.axes.get_ylim()[1]))

        self.lineEdit_vmin.setText(str(self.canvas.fitsplot.get_clim()[0]))
        self.lineEdit_vmax.setText(str(self.canvas.fitsplot.get_clim()[1]))

    def update_cursor_position(self, event):
        """Updates the cursor position textbox when the cursor moves within the
        matplotlib axis.

        Parameters
        ----------
        event : QEvent
            The event that occurred within the signalling widget
        """
        if event.inaxes:
            self.lineEdit_cursorPosition.setText('({:.0f}, {:.0f})'.format(event.xdata,
                                                                           event.ydata))
            self.lineEdit_cursorValue.setText('{:.0f}'.format(self.data[int(event.ydata),
                                                                        int(event.xdata)]))

    def update_profile(self, event):
        """Updates the profile plot when the cursor moves within the
        matplotlib axis.

        Parameters
        ----------
        event : QEvent
            The event that occurred within the signalling widget
        """
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
            self.canvas_profile.countrate_label.set_text(
                '3 x 3 countrate: {:.0f}'.
                format(countrate)
            )

            self.canvas_profile.draw()

    def show_segment_id(self, event):
        """Show the segment ID of a segment when the cursor mouses over
        a segment on the matplotlib axis.

        Parameters
        ----------
        event : QEvent
            The event that occurred within the signalling widget
        """
        if event.inaxes:
            # Determine if the cursor is close to a segment
            dist = np.sqrt((self.x - event.xdata)**2 + (self.y - event.ydata)**2)
            indseq = np.nonzero(np.equal(dist, np.amin(dist)))[0]
            ind = indseq[0]
            #FIXME! - add real numbering
            segnum = next((name for name, seg in self.wss_dict.items() if seg['coords'] == (self.x[ind], self.y[ind])), None)

            if dist[ind] <= self.epsilon:
                # Update the tool tip to show the segment ID
                #self.canvas.setToolTip(str(ind + 1))
                self.canvas.setToolTip(str(segnum))

            else:
                # Update the tool tip to be blank
                self.canvas.setToolTip(str(''))

    def button_press_callback(self, event):
        """Whenever the mouse is clicked, determines if the click is
        within the matplotlib axis; if so calls get_ind_under_point.

        Parameters
        ----------
        event : QEvent
            The event that occurred within the signalling widget
        """
        if not event.inaxes:
            return
        if event.button != 1:
            return

        self._ind = self.get_ind_under_point(event)


    def set_guidestar(self):
        """Set the currently selected star to be the guide star.

        Change order of self.inds to reflect the new guide star, update
        list of original index order, and update GUI plotted circles
        """
        # Determine index of old guide star
        ind_of_old_guide_star = self.gs_ind

        # Re-draw old guide star as regular star
        self.circles[ind_of_old_guide_star].set_markeredgecolor('darkorange')

        # Determine index of new guide star
        guide_ind = self.tableWidget_selectedStars.currentRow()
        self.gs_ind = guide_ind

        # Re-draw new guide star as guide star
        self.circles[guide_ind].set_markeredgecolor('yellow')

        # Move the guide star icon
        self.tableWidget_selectedStars.setItem(ind_of_old_guide_star, 0,
                                               QTableWidgetItem(''))
        icon = QIcon(os.path.join(__location__, 'gs_icon.png'))
        self.tableWidget_selectedStars.setItem(guide_ind, 0,
                                               QTableWidgetItem(icon, ''))

        self.canvas.draw()

    def remove_star(self):
        """Remove the currently selected star and update inds list.
        """

        # Determine index of star being removed
        star_ind = self.tableWidget_selectedStars.currentRow()

        # Un-draw circle
        delete_circle = self.circles.pop(star_ind)
        self.canvas.axes.lines.remove(delete_circle)

        # Send deletion message to output
        ind_of_ind = int(self.tableWidget_selectedStars.item(star_ind, 1).text()) - 1
        remstar_string = 'Deleted star at: x={:.1f}, y={:.1f}'.format(self.x[ind_of_ind],
                                                                      self.y[ind_of_ind])
        red_text = "<span style=\" color:#ff0000;\" >" + \
                   remstar_string + "</span><br>"
        self.textEdit_output.setHtml(red_text + self.textEdit_output.toHtml())
        if self.print_output:
            print(remstar_string)

        # Update inds list
        self.inds.pop(star_ind)

        # Update guide star index
        gs_deleted, self.gs_ind = self.calc_new_gs_ind(star_ind)

        # Update table
        self.tableWidget_selectedStars.removeRow(star_ind)
        if self.gs_ind is None:
            # Add back empty row
            self.tableWidget_selectedStars.insertRow(0)
        # Update which row is highlighted:
        if self.gs_ind is not None:
            self.current_row = self.tableWidget_selectedStars.currentRow()
            self.circles[self.current_row].set_markeredgecolor('cornflowerblue')

        # Update guide star
        if gs_deleted and self.gs_ind is not None:
            # Re-draw new guide star as guide star
            self.circles[self.gs_ind].set_markeredgecolor('yellow')

            # Move the guide star icon
            icon = QIcon(os.path.join(__location__, 'gs_icon.png'))
            self.tableWidget_selectedStars.setItem(self.gs_ind, 0,
                                                   QTableWidgetItem(icon, ''))

        self.canvas.draw()

    def clear_selected_stars(self):
        """Remove all stars from the table.
        """

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
        red_text = "<span style=\" color:#ff0000;\" >" + \
                   remstar_string + "</span><br>"
        self.textEdit_output.setHtml(red_text + self.textEdit_output.toHtml())
        if self.print_output:
            print(remstar_string)

    def quit_star_selection(self):
        """Closes the star selector window.
        """

        self.answer = True
        # If the user didn't choose any stars, ask if they really want to quit.
        if not self.inds:
            no_stars_selected_dialog = QMessageBox()
            no_stars_selected_dialog.setText('No stars selected' + ' ' * 50)
            no_stars_selected_dialog.setInformativeText(
                'The tool will not be able to continue. Do you want to quit anyway?'
            )
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
        """Closes the star selector window and clears indices.
        """

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
        """Get the index of the closest star to where a click occurred.

        If user clicks within one epsilon of an identified star, appends that
        star's position to inds; notifies user if no star is nearby or if
        selected star is already selected.

        Parameters
        ----------
        event : QEvent
            The event that occurred within the signalling widget

        Returns
        -------
        ind : int
            The index of the star closest to a mouseButtonPress event

        """

        dist = np.sqrt((self.x - event.xdata)**2 + (self.y - event.ydata)**2)
        # TODO - this is probably where we want to check the number of the segment
        #FIXME!self.wss_dict
        # Index of the coord in the list of coords
        indseq = np.nonzero(np.equal(dist, np.amin(dist)))[0]
        ind = indseq[0] #Index of the x,y that is closest to where the user clicked
        #WSS segumber
        if len(dist) <= 16:
            LOGGER.warning("Do not have enough PSFs to match to WSS numberings.")
            segnum = indseq[0]
        else:
            segnum = next((name for name, seg in self.wss_dict.items() if seg['coords'] == (self.x[ind], self.y[ind])), None)

        # If the user has already selected the maximum number of stars
        if len(self.inds) == self.n_stars_max:
            if self.print_output:
                print('Maximum number of stars, {}, already selected.'.format(self.n_stars_max))

            red_text = "<span style=\" color:#ff0000;\" >"
            red_text += 'Maximum number of stars, {}, already selected.'.format(self.n_stars_max)
            red_text += "</span><br>"

            self.textEdit_output.setHtml(red_text + self.textEdit_output.toHtml())
            return

        # If there is no star within ``dist`` of the mouse click
        elif dist[ind] >= self.epsilon:
            if self.print_output:
                print('No star within {} pixels. No star selected.'.format(self.epsilon))

            red_text = "<span style=\" color:#ff0000;\" >"
            red_text += 'No star within {} pixels. No star selected.'.format(self.epsilon)
            red_text += "</span><br>"

            self.textEdit_output.setHtml(red_text + self.textEdit_output.toHtml())

        # If the user clicked on a star that already exists
        elif ind in self.inds:
            if self.print_output:
                print('Star already selected, please choose another star')

            red_text = "<span style=\" color:#ff0000;\" >"
            red_text += 'Star already selected, please choose another star'
            red_text += "</span><br>"

            self.textEdit_output.setHtml(red_text + self.textEdit_output.toHtml())

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
            values = ['', str(ind + 1), str(int(self.x[ind])),
                      str(int(self.y[ind])), str(int(countrate))]
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
            #self.inds.append(ind)
            self.inds.append(segnum)
            #FIXME!
            #self.inds.append(segnum)

        return segnum #ind

    def reorder_indices(self):
        """Re-order self.inds such that the desired guide star is the
        first element on the list.
        """

        # Determine which ID is the guide star
        guide_star_id = self.inds[self.gs_ind]
        # FIXME - rearrange segnums... this code should actually be sufficient

        # Remove it from the list
        self.inds.pop(self.gs_ind)

        # Rewrite self.inds
        inds_final = [guide_star_id] + self.inds
        self.inds = inds_final

    def nostars_dialog(self, button):
        """Process dialog box asking if user wants to close the window
        without selecting stars.

        Parameters
        ----------
        button : QPushButton
            The button that was pushed in the dialog box.
        """
        if 'No' in button.text():
            self.answer = False

    def calc_new_gs_ind(self, star_ind):
        """Determine the new guide star index after a star was removed.

        Parameters
        ----------
        star_ind : int
            Index of the deleted star

        Returns
        -------
        bool
            Denotes if the guide star was deleted
        int
            The new guide star index
        """
        # If a star earlier than the guide star was deleted
        if star_ind < self.gs_ind:
            return False, self.gs_ind - 1

        # If the guide star was deleted
        elif star_ind == self.gs_ind:
            if len(self.inds) == 0:
                return True, None

            elif star_ind == 0:
                # Make the first star the guide star
                return True, 0
            else:
                # Make the star_ind - 1 star the guide star
                return True, star_ind - 1

        # If a star after the guide star was deleted
        elif star_ind > self.gs_ind:
            return False, self.gs_ind

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# MAIN FUNCTION
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


def run_SelectStars(data, x, y, dist, print_output=False, masterGUIapp=None):
    """Calls a PyQt GUI to allow interactive user selection of guide and
    reference stars.

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
    masterGUIapp : qApplication, optional
        qApplication instance of parent GUI

    Returns
    -------
    inds : list
        List of indices of positions of selected stars
    """

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
    print(inds)
    return inds
