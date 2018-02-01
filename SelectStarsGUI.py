# Built using template at:
# https://matplotlib.org/examples/user_interfaces/embedding_in_qt5.html

# Standard Library
from __future__ import unicode_literals
import sys

# Third Party
import numpy as np
import matplotlib
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure
from matplotlib.colors import LogNorm
import matplotlib.pyplot as plt
from PyQt5 import QtCore, QtWidgets
from PyQt5.QtWidgets import QMainWindow, QWidget, QVBoxLayout, QApplication, \
                            QSizePolicy, QFormLayout, QPushButton, QLineEdit, \
                            QLabel, QTextEdit, QRadioButton, QGroupBox, \
                            QGridLayout
from PyQt5.QtCore import pyqtSlot
from astropy.io import fits

matplotlib.rcParams['font.family'] = 'serif'
matplotlib.rcParams['font.weight'] = 'light'
matplotlib.rcParams['mathtext.bf'] = 'serif:normal'


class MyMplCanvas(FigureCanvas):
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

    def __init__(self, parent=None, width=15, height=15, dpi=100, data=None,
                 x=None, y=None, left=None, right=None, bottom=0.05):

        self.data = data
        self.x = x
        self.y = y

        self.fig = Figure(figsize=(width, height), dpi=dpi)
        self.axes = self.fig.add_subplot(111)
        self.fig.subplots_adjust(top=.95, left=left, right=right,
                                 bottom=bottom)

        FigureCanvas.__init__(self, self.fig)
        self.setParent(parent)

        # self.compute_initial_figure(self.fig, data, x, y)

        # FigureCanvas.setSizePolicy(self,
        #                            QtWidgets.QSizePolicy.Expanding,
        #                            QtWidgets.QSizePolicy.Expanding)
        # FigureCanvas.updateGeometry(self)

    def compute_initial_figure(self, fig, data, x, y):
        self.fitsplot = self.axes.imshow(data, cmap='bone', interpolation='nearest',
                                         clim=(0.1, 1e3), norm=LogNorm())
        self.axes.scatter(x, y, c='r', marker='+')
        self.fig.colorbar(self.fitsplot, ax=self.axes, fraction=0.046, pad=0.04)
        self.draw()

    def init_profile(self):
        # x_range = range(2048)
        # self.y_slices = [self.data[i] for i in range(2048)]
        # self.lines = [self.axes.plot(x_range, y, c='white') for y in self.y_slices]
        self.axes.set_ylim(np.min(self.data) / 5, 5 * np.max(self.data))
        self.axes.set_yscale('log')
        self.axes.set_ylabel('Counts')
        self.axes.set_xlabel('X Pixels')
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

        self.axes.set_ylim(min(2048, x_mid + ax_range/2),
                           max(0, x_mid - ax_range/2))
        self.axes.set_xlim(max(0, y_mid - ax_range/2),
                           min(2048, y_mid + ax_range/2))
        self.draw()


class ApplicationWindow(QtWidgets.QMainWindow):
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
    def __init__(self, data=None, x=None, y=None, dist=None, qApp=None,
                 print_output=False):
        '''Defines attributes; calls initUI() method to set up user interface.'''
        self.qApp = qApp
        self.print_output = print_output

        self.title = 'PyQt5 matplotlib example - pythonspot.com'
        self.image_dim = 800

        self.data = data
        self.x = x
        self.y = y

        self._ind = None
        self.inds = []
        self.inds_of_inds = []
        self.epsilon = dist

        # Initialize user interface
        self.initUI()

    def initUI(self):
        '''Sets up interactive graphical user interface.
        '''

        # Set window attributes
        QMainWindow.__init__(self)
        # self.setAttribute(QtCore.Qt.WA_DeleteOnClose)
        self.setWindowTitle("FGS Guide & Reference Star Selector")
        self.main_widget = QWidget(self)
        mainGrid = QGridLayout()  # set grid layout
        self.main_widget.setLayout(mainGrid)
        self.main_widget.setFocus()
        self.setCentralWidget(self.main_widget)

        # Add plot - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        sc = MyMplCanvas(self.main_widget, width=5, height=4, dpi=100,
                         data=self.data, x=self.x, y=self.y, left=0.1, right=0.9)
        self.canvas = sc
        self.canvas.compute_initial_figure(self.canvas.fig, self.data, self.x,
                                           self.y)
        self.canvas.zoom_to_crop()
        mainGrid.addWidget(sc, 0, 0, 4, 2)

        mainGrid.setColumnMinimumWidth(0, self.image_dim - 150)

        # Show current cursor position
        self.cursor_label = QLabel('Cursor Position:', self,
                                   alignment=QtCore.Qt.AlignRight)
        mainGrid.addWidget(self.cursor_label, 4, 0)

        self.cursor_textbox = QLineEdit(self)
        self.cursor_textbox.setPlaceholderText('Move cursor into axes')
        self.cursor_textbox.setMinimumSize(200, 20)
        mainGrid.addWidget(self.cursor_textbox, 4, 1)

        self.canvas.mpl_connect('motion_notify_event',
                                self.update_cursor_position)

        # Show value under cursor
        self.pixel_label = QLabel('Pixel Value:', self,
                                   alignment=QtCore.Qt.AlignRight)
        mainGrid.addWidget(self.pixel_label, 5, 0)

        self.pixel_textbox = QLineEdit(self)
        self.pixel_textbox.setPlaceholderText('Move cursor into axes')
        self.pixel_textbox.setMinimumSize(200, 20)
        mainGrid.addWidget(self.pixel_textbox, 5, 1)

        # Update cursor position and pixel value under cursor
        self.canvas.mpl_connect('motion_notify_event',
                                self.update_cursor_position)

        # Star selection!
        self.canvas.mpl_connect('button_press_event',
                                self.button_press_callback)

        # Add first column  - - - - - - - - - - - - - - - - - - - - - - - - - -

        # Make Group Box to put the label in, because it will misbehave
        instructionsGroupBox = QGroupBox('Instructions', self)
        vBox = QVBoxLayout()

        self.instructions = QLabel('''Click as near to the center of the star \
        as possible. The first star that is choosen will be the <span style=\
        "background-color: #FFFF00">guide star</span>. All additional stars that are \
        clicked on are the <span style="background-color: #FFA500">reference stars\
        </span>. Star locations will appear in the list at right as they are selected; \
        the guide star can be re-selected using the radio buttons. Errors will be shown\
        in the output box below.''', self)
        # self.instructions.setSizePolicy(QSizePolicy.Preferred, QSizePolicy.Maximum)

        self.instructions.setWordWrap(True)

        vBox.addWidget(self.instructions)
        instructionsGroupBox.setLayout(vBox)
        instructionsGroupBox.setSizePolicy(QSizePolicy.Preferred,
                                           QSizePolicy.Maximum)
        mainGrid.addWidget(instructionsGroupBox, 0, 2, 1, 3,
                           alignment=QtCore.Qt.AlignTop)

        # Log to update
        self.log_textbox = QTextEdit(self)
        self.log_textbox.setPlaceholderText('No stars selected.')
        # self.log_textbox.setSizePolicy(QSizePolicy.Preferred, QSizePolicy.MinimumExpanding)
        self.log_textbox.setMinimumSize(10, 300)
        mainGrid.setRowMinimumHeight(1, 320)
        mainGrid.addWidget(self.log_textbox, 1, 2, 1, 2,
                           alignment=QtCore.Qt.AlignTop)

        # mainGrid.setRowMinimumHeight(1, self.image_dim*.2)

        # Plot slice of profile under cursor
        prof = MyMplCanvas(self.main_widget, width=4, height=3, dpi=100,
                           data=self.data, left=0.2, right=0.95, bottom=0.2)
        prof.init_profile()
        self.profile = prof
        self.canvas.mpl_connect('motion_notify_event', self.update_profile)
        mainGrid.setRowMinimumHeight(2, self.image_dim * .25)
        # mainGrid.setColumnMinimumWidth(2, self.image_dim * .5)
        mainGrid.addWidget(self.profile, 2, 2, 1, 2,
                           alignment=QtCore.Qt.AlignVCenter)

        # Create colorbar-updating section
        cbarGroupBox = QGroupBox('Colorbar Limits', self)
        cbarGrid = QGridLayout()
        cbarGroupBox.setLayout(cbarGrid)

        cbarGrid.addWidget(QLabel('Min Value:  ', self), 0, 0)
        cbarGrid.addWidget(QLabel('Max Value:  ', self), 1, 0)

        self.vmin_textbox = QLineEdit(str(sc.fitsplot.get_clim()[0]), self)
        cbarGrid.addWidget(self.vmin_textbox, 0, 1)

        self.vmax_textbox = QLineEdit(str(sc.fitsplot.get_clim()[1]), self)
        cbarGrid.addWidget(self.vmax_textbox, 1, 1)

        self.cbarlims_button = QPushButton('Update colorbar limits', self)
        self.cbarlims_button.clicked.connect(self.update_cbar)   # connect button to function on_click
        cbarGrid.addWidget(self.cbarlims_button, 2, 0, 1, 2)

        mainGrid.addWidget(cbarGroupBox, 3, 3, 3, 1)

        # Create axis-updating section
        axGroupBox = QGroupBox('Axes Limits', self)
        axGrid = QGridLayout()
        axGroupBox.setLayout(axGrid)

        axGrid.addWidget(QLabel('X: ( ', self), 1, 0)
        axGrid.addWidget(QLabel(',', self), 1, 2)
        axGrid.addWidget(QLabel(' )', self), 1, 4)

        axGrid.addWidget(QLabel('Y: ( ', self), 2, 0)
        axGrid.addWidget(QLabel(',', self), 2, 2)
        axGrid.addWidget(QLabel(' )', self), 2, 4)

        self.x1_textbox = QLineEdit(str(sc.axes.get_xlim()[0]), self)
        axGrid.addWidget(self.x1_textbox, 1, 1)

        self.x2_textbox = QLineEdit(str(sc.axes.get_xlim()[1]), self)
        axGrid.addWidget(self.x2_textbox, 1, 3)

        self.y1_textbox = QLineEdit(str(sc.axes.get_ylim()[1]), self)
        axGrid.addWidget(self.y1_textbox, 2, 1)

        self.y2_textbox = QLineEdit(str(sc.axes.get_ylim()[0]), self)
        axGrid.addWidget(self.y2_textbox, 2, 3)

        self.axlims_button = QPushButton('Update axis limits', self)
        self.axlims_button.clicked.connect(self.update_axes)   # connect button to function on_click
        axGrid.addWidget(self.axlims_button, 3, 0, 1, 5)

        # Add "Zoom Fit" button
        self.zoom_button = QPushButton('Zoom Fit', self)
        self.zoom_button.setToolTip('Zoom to encompass all data')
        self.zoom_button.clicked.connect(sc.zoom_to_fit)
        self.zoom_button.clicked.connect(self.update_textboxes)
        axGrid.addWidget(self.zoom_button, 4, 0, 1, 5)

        # Add "Crop to Data" button
        self.crop_button = QPushButton('Crop to Data', self)
        self.crop_button.setToolTip('Zoom crop to data')
        self.crop_button.clicked.connect(sc.zoom_to_crop)
        self.crop_button.clicked.connect(self.update_textboxes)
        axGrid.addWidget(self.crop_button, 5, 0, 1, 5)

        mainGrid.addWidget(axGroupBox, 3, 2, 3, 1)

        # Add second column  - - - - - - - - - - - - - - - - - - - - - - - - - -

        starsGroupBox = QGroupBox(self)
        col2Grid = QGridLayout()

        # Show selected stars
        guide_label = QLabel('Guide\nStar?', self)
        guide_label.setMinimumSize(30, 30)
        guide_label.setSizePolicy(QSizePolicy.Preferred, QSizePolicy.Maximum)
        col2Grid.addWidget(guide_label, 0, 0)
        col2Grid.addWidget(QLabel('Star Position', self), 0, 1)
        col2Grid.setColumnMinimumWidth(0, 50)
        col2Grid.setColumnMinimumWidth(1, 120)
        self.nStars = 11

        self.star_positions = []
        self.guidestar_buttons = []
        for n in range(self.nStars):
            # Create radio buttons
            guidestar_button = QRadioButton(starsGroupBox)
            if n == 0: guidestar_button.setChecked(True)
            guidestar_button.released.connect(self.make_setguidestar(guidestar_button))
            if n == 0: guidestar_button.setEnabled(False)
            if n != 0: guidestar_button.setCheckable(False)

            # Create star position textboxes
            star_position = QLineEdit(starsGroupBox)
            if n == 0: star_position.setPlaceholderText('No stars selected')

            # Create 'remove' buttons
            remove_button = QPushButton('Delete', starsGroupBox)
            remove_button.clicked.connect(self.make_removestar(remove_button))
            remove_button.installEventFilter(self)
            remove_button.setEnabled(False)

            # Add to starsGroupBox grid layout
            col2Grid.addWidget(guidestar_button, n + 1, 0)
            col2Grid.addWidget(star_position, n + 1, 1)
            col2Grid.addWidget(remove_button, n + 1, 2)

        self.star_positions = [child for child in starsGroupBox.findChildren(QLineEdit)]
        self.guidestar_buttons = [child for child in starsGroupBox.findChildren(QRadioButton)]
        self.remove_buttons = [child for child in starsGroupBox.findChildren(QPushButton)]

        starsGroupBox.setLayout(col2Grid)
        starsGroupBox.setSizePolicy(QSizePolicy.Preferred, QSizePolicy.Maximum)
        mainGrid.addWidget(starsGroupBox, 1, 4, 2, 1)

        # Add "Done" button
        self.done_button = QPushButton('Done', self)
        self.done_button.setToolTip('Close the window')
        self.done_button.clicked.connect(self.fileQuit)
        self.done_button.setMinimumSize(150, 50)
        mainGrid.addWidget(self.done_button, 3, 4, 1, 1,
                           alignment=QtCore.Qt.AlignBottom)

        # Add "Cancel" button
        self.cancel_button = QPushButton('Cancel', self)
        self.cancel_button.setToolTip('Close the window and discard changes')
        self.cancel_button.clicked.connect(self.cancel)
        self.cancel_button.setMinimumSize(150, 50)
        mainGrid.addWidget(self.cancel_button, 4, 4, 2, 1,
                           alignment=QtCore.Qt.AlignVCenter)

    @pyqtSlot()
    def eventFilter(self, remove_button, event):
        '''Parse out the cursor moving in or out of a star deletion button, and
        updating the matplotlib axis with a red highlighted circle accordingly'''
        if event.type() == QtCore.QEvent.Enter:
            if remove_button.isEnabled() == True:
                # Determine index of star corresponding to button
                star_ind = self.remove_buttons.index(remove_button)
                ind_of_star_ind = self.inds_of_inds.index(star_ind)

                # Re-draw selected star as red
                self.canvas.axes.lines[ind_of_star_ind].set_markeredgecolor('red')

                self.canvas.draw()
            return True

        if event.type() == QtCore.QEvent.Leave:
            if remove_button.isEnabled() == True:
                # Determine index of star corresponding to button
                star_ind = self.remove_buttons.index(remove_button)
                ind_of_star_ind = self.inds_of_inds.index(star_ind)

                # Re-draw selected star as original color
                if ind_of_star_ind == 0:
                    self.canvas.axes.lines[ind_of_star_ind].set_markeredgecolor('yellow')
                else:
                    self.canvas.axes.lines[ind_of_star_ind].set_markeredgecolor('darkorange')

                self.canvas.draw()
            return True

        return False

    def update_axes(self):
        '''Changes the axes limits of the matplotlib canvas to the current
        values of the axis limit textboxes
        '''
        self.canvas.axes.set_xlim(float(self.x1_textbox.text()),
                                  float(self.x2_textbox.text()))
        self.canvas.axes.set_ylim(float(self.y2_textbox.text()),
                                  float(self.y1_textbox.text()))
        self.canvas.fig.subplots_adjust(top=.95, left=.07)
        self.canvas.draw()

    def update_cbar(self):
        '''Changes the limits of the colorbar to the current values of the
        colorbar limit textboxes
        '''
        self.canvas.fitsplot.set_clim(float(self.vmin_textbox.text()),
                                      float(self.vmax_textbox.text()))
        self.canvas.draw()

    def update_textboxes(self):
        '''Changes the values of the axis limit textboxes to the current axis
        limits of the matplotlib canvas'''
        self.x1_textbox.setText(str(self.canvas.axes.get_xlim()[0]))
        self.x2_textbox.setText(str(self.canvas.axes.get_xlim()[1]))
        self.y1_textbox.setText(str(self.canvas.axes.get_ylim()[1]))
        self.y2_textbox.setText(str(self.canvas.axes.get_ylim()[0]))

    def update_cursor_position(self, event):
        '''Updates the cursor position textbox when the cursor moves within the
        matplotlib axis'''
        if event.inaxes:
            self.cursor_textbox.setText('({:.0f}, {:.0f})'.format(event.xdata,
                                                                  event.ydata))
            self.pixel_textbox.setText('{:.2f}'.format(self.data[int(event.ydata),
                                                                 int(event.xdata)]))

    def update_profile(self, event):
        '''Updates the profile plot when the cursor moves within the matplotlib axis'''
        if event.inaxes:
            x = np.arange(max(0, np.floor(event.xdata - self.epsilon)),
                          min(2047, np.ceil(event.xdata + self.epsilon))).astype(int)
            y = self.data[int(event.ydata), x]

            # Make previous lines white so only the current profile is visible
            for line in self.profile.axes.lines:
                line.set_color('white')

            # self.profile_line[0].set_color('white')
            # self.profile_line = self.profile.lines[int(event.ydata)]
            # self.profile_line[0].set_color('cornflowerblue')
            # self.profile_line[0].set_zorder(2050)
            self.profile.axes.plot(x, y, c='cornflowerblue')
            self.profile.axes.set_xlim(int(np.floor(event.xdata - self.epsilon)),
                                       int(np.ceil(event.xdata + self.epsilon)))
            self.profile.draw()

    def button_press_callback(self, event):
        '''Whenever a mouse button is pressed, determines if a button press is
        within the matplotlib axis; if so calls get_ind_under_point'''
        if not event.inaxes:
            return
        if event.button != 1:
            return
        self._ind = self.get_ind_under_point(event)

    def get_ind_under_point(self, event):
        '''If user clicks within one epsilon of an identified star, appends that
        star's position to inds; notifies user if no star is nearby or if
        selected star is already selected'''

        d = np.sqrt((self.x - event.xdata)**2 + (self.y - event.ydata)**2)
        # i = np.unravel_index(np.nanargmin(d), d.shape)
        indseq = np.nonzero(np.equal(d, np.amin(d)))[0]
        ind = indseq[0]

        if len(self.inds) == self.nStars:
            if self.print_output:
                print('Maximum number of stars, {}, already selected.'.format(self.nStars))

            redText = "<br/><span style=\" color:#ff0000;\" >"
            redText += 'Maximum number of stars, {}, already selected.'.format(self.nStars)
            redText += "</span>"

            self.log_textbox.setHtml(redText + self.log_textbox.toHtml())
            return

        if d[ind] >= self.epsilon:
            if self.print_output:
                print('No star within {} pixels. No star selected.'.format(self.epsilon))

            redText = "<br/><span style=\" color:#ff0000;\" >"
            redText += 'No star within {} pixels. No star selected.'.format(self.epsilon)
            redText += "</span>"

            self.log_textbox.setHtml(redText + self.log_textbox.toHtml())

        elif ind in self.inds:
            if self.print_output:
                print('Star already selected, please choose another star')

            redText = "<br/><span style=\" color:#ff0000;\" >"
            redText += 'Star already selected, please choose another star'
            redText += "</span>"

            self.log_textbox.setHtml(redText + self.log_textbox.toHtml())

        else:
            if len(self.inds) == 0:  # First star selected
                c = 'yellow'
                newstar_string = '1 star selected: x={:.1f}, y={:.1f}'.format(self.x[ind],
                                                                              self.y[ind])
                self.guidestar_buttons[0].setEnabled(True)
            else:  # >= second star selected
                c = 'darkorange'
                newstar_string = '{} stars selected: x={:.1f}, y={:.1f}'.format(len(self.inds) + 1,
                                                                                self.x[ind],
                                                                                self.y[ind])

            if self.print_output:
                print(newstar_string)

            self.log_textbox.setHtml('<br/>' + newstar_string + self.log_textbox.toHtml())
            self.star_positions[len(self.inds)].setText('x={:.1f}, y={:.1f}'.format(self.x[ind],
                                                                                    self.y[ind]))
            self.guidestar_buttons[len(self.inds)].setCheckable(True)
            self.remove_buttons[len(self.inds)].setEnabled(True)
            self.inds_of_inds.append(len(self.inds))

            self.canvas.axes.plot(self.x[ind], self.y[ind], 'o', ms=25,
                                  mfc='none', mec=c, mew=2, lw=0)
            self.canvas.draw()
            self.inds.append(ind)

        return ind

    def make_setguidestar(self, button):
        # Have to define second function inside to return a function to the toggle.connect method
        def setguidestar():
            '''
            Change order of self.inds to reflect the new guide star
            Update list of original index order
            Update GUI plotted circles
            '''

            # Determine index of old guide star
            ind_of_OLD_guide = self.inds_of_inds[0]

            # Re-draw old guide star as regular star
            self.canvas.axes.lines[ind_of_OLD_guide].set_markeredgecolor('darkorange')

            # Determine index of new guide star
            guide_ind = self.guidestar_buttons.index(button)
            ind_of_guide_ind = self.inds_of_inds.index(guide_ind)

            self.inds.insert(0, self.inds.pop(ind_of_guide_ind))
            self.inds_of_inds.insert(0, self.inds_of_inds.pop(ind_of_guide_ind))

            # Re-draw new guide star as guide star
            self.canvas.axes.lines[guide_ind].set_markeredgecolor('yellow')

            self.canvas.draw()

        return setguidestar

    def make_removestar(self, button):
        # Have to define second function inside to return a function to the toggle.connect method
        def removestar():
            '''
            Remove star and update inds list
            '''

            # Determine index of star being removed
            star_ind = self.remove_buttons.index(button)
            ind_of_star_ind = self.inds_of_inds.index(star_ind)

            # If user tries to delete guide star, don't let them
            if self.guidestar_buttons[star_ind].isChecked():
                redText = "<br/><span style=\" color:#ff0000;\" >"
                redText += 'The guide star cannot be deleted. Please choose \
                            another star or change the guide star.'
                redText += "</span>"
                self.log_textbox.setHtml(redText + self.log_textbox.toHtml())
                return

            # Un-draw circle
            l_rem = self.canvas.axes.lines.pop(star_ind)
            del l_rem

            # Send deletion message to output
            remstar_string = 'Deleted star at: x={:.1f}, y={:.1f}'.format(self.x[self.inds[ind_of_star_ind]],
                                                                          self.y[self.inds[ind_of_star_ind]])
            redText = "<br/><span style=\" color:#ff0000;\" >" + \
                      remstar_string + "</span>"
            self.log_textbox.setHtml(redText + self.log_textbox.toHtml())
            if self.print_output:
                print(remstar_string)

            # Update inds list
            self.inds.pop(ind_of_star_ind)
            removed_ind = len(self.inds)

            # Update inds of inds list
            self.inds_of_inds.pop(ind_of_star_ind)
            for i in range(len(self.inds_of_inds)):
                if self.inds_of_inds[i] > star_ind:
                    self.inds_of_inds[i] -= 1

            # Update star position textboxes
            for n in np.arange(star_ind, removed_ind):
                self.star_positions[n].setText(self.star_positions[n + 1].text())
            self.star_positions[removed_ind].setText('')

            # Update guide star radio buttons
            ind_of_guide = self.inds_of_inds[0]
            self.guidestar_buttons[removed_ind].setCheckable(False)
            self.guidestar_buttons[ind_of_guide].setChecked(True)

            # Update remove buttons
            self.remove_buttons[removed_ind].setEnabled(False)

            self.canvas.draw()

        return removestar


    def fileQuit(self):
        '''Closes the application'''
        # if self.inds == []:

        self.qApp.exit(0)  # Works only with self.close() after; same as qApp.quit()
        self.close()

    def cancel(self):
        '''Closes the application and clears indices'''
        self.inds = []
        self.qApp.exit(0)  # Works only with self.close() after; same as qApp.quit()
        self.close()

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

def run_SelectStars(data, x, y, dist, print_output=False):
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
    # plt.close()
    qApp = QApplication(sys.argv)
    aw = ApplicationWindow(data=data, x=x, y=y, dist=dist, qApp=qApp,
                           print_output=print_output)
    aw.show()

    try:
        plt.get_current_fig_manager().window.raise_()  # Bring window to front
    except AttributeError:
        pass

    inds = aw.inds
    qApp.exec_()  # Begin interactive session; pauses until qApp.exit() is called

    return inds


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

if __name__ == '__main__':
    data = fits.open('../data/LMCfootprint_80/LMCfootprint_80.31512449419834_-70.45278790693078.fits')[0].data
    data[data == 0] = 0.1  # Adjust so 0 shows up as such with a logNorm colorbar
    gauss_sigma = 5
    dist = 16

    inds = run_SelectStars(data, gauss_sigma, dist)
