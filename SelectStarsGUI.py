# Built using template at https://matplotlib.org/examples/user_interfaces/embedding_in_qt5.html

#Standard Library
from __future__ import unicode_literals
import sys
import os
import random

# Third Party
import numpy as np
import matplotlib
#matplotlib.use('Qt5Agg') # Make sure that we are using Qt5
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure
from matplotlib.colors import LogNorm
import matplotlib.pyplot as plt
from PyQt5 import QtCore, QtWidgets, QtGui
from PyQt5.QtWidgets import QMainWindow, QWidget, QVBoxLayout, QApplication, QSizePolicy, QFormLayout, \
                            QPushButton, QLineEdit, QLabel, QTextEdit, QRadioButton, QGroupBox, QGridLayout
from PyQt5.QtCore import pyqtSlot
from astropy.io import fits
from astropy.stats import sigma_clipped_stats
from scipy import ndimage
from photutils import find_peaks


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

    def __init__(self, parent=None, width=15, height=15, dpi=100, data=None, x=None, y=None):

        self.data = data

        self.fig = Figure(figsize=(width, height), dpi=dpi)
        self.axes = self.fig.add_subplot(111)
        self.fig.subplots_adjust(top=.95, left = .07)

        FigureCanvas.__init__(self, self.fig)
        self.setParent(parent)

        self.compute_initial_figure(self.fig, data, x, y)
        self.zoom_to_fit()

        # FigureCanvas.setSizePolicy(self,
        #                            QtWidgets.QSizePolicy.Expanding,
        #                            QtWidgets.QSizePolicy.Expanding)
        # FigureCanvas.updateGeometry(self)

    def compute_initial_figure(self, fig, data, x, y):
        fitsplot = self.axes.imshow(data, cmap='bone', interpolation='nearest' , clim=(0.1,100), norm=LogNorm())
        self.axes.scatter(x, y, c='r', marker='+')
        self.fig.colorbar(fitsplot, ax=self.axes, fraction=0.046, pad=0.04)

    def zoom_to_fit(self):
        self.axes.set_xlim(0, np.shape(self.data)[0])
        self.axes.set_ylim(np.shape(self.data)[1], 0)
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
    printOutput (boolean)
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
        Furthest distance from a star, in pixels, that a user can click to select a star


    Methods
    -------
    initUI
        Sets up user interface with all widgets listed above
    update_axes
        Changes the axes limits of the matplotlib canvas to the current values of the axis limit textboxes
    update_textboxes
        Changes the values of the axis limit textboxes to the current axis limits of the matplotlib canvas
    move_in_axes
        Updates the cursor position textbox when the cursor moves within the matplotlib axis
    button_press_callback
        Determines if a button press is within the matplotlib axis; if so calls get_ind_under_point
    get_ind_under_point
        If user clicks within one epsilon of an identified star, appends that star's position to inds;
        notifies user if no star is nearby or if selected star is already selected
    make_setGuideStar
        Updates the order of inds (i.e. the selected guide star) if user changes the guide star and
        updates matplotlib axis accordingly
    fileQuit
        Closes the application
    """
    def __init__(self, data=None, x=None, y=None, dist = None, qApp=None, printOutput=False):
        '''Defines attributes; calls initUI() method to set up user interface.'''
        self.qApp = qApp
        self.printOutput = printOutput

        self.title = 'PyQt5 matplotlib example - pythonspot.com'
        self.image_dim = 800

        self.data = data
        self.x = x
        self.y = y

        self._ind = None
        self.inds = []
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

        mainGrid = QGridLayout() # set grid layout
        self.main_widget.setLayout(mainGrid)
        self.main_widget.setFocus()
        self.setCentralWidget(self.main_widget)

        # Add plot - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        sc = MyMplCanvas(self.main_widget, width=5, height=4, dpi=100, data=self.data, x=self.x, y=self.y)
        self.canvas = sc
        mainGrid.addWidget(sc, 0, 0, 3, 2)

        mainGrid.setColumnMinimumWidth(0, self.image_dim - 150)

        # Show current cursor position
        self.cursor_label = QLabel('Cursor Position:', self, alignment=QtCore.Qt.AlignRight)
        mainGrid.addWidget(self.cursor_label, 3, 0)

        self.cursor_textbox = QLineEdit(self)
        self.cursor_textbox.setPlaceholderText('Move cursor into axes')
        mainGrid.addWidget(self.cursor_textbox, 3, 1)

        self.canvas.mpl_connect('motion_notify_event', self.move_in_axes)

        # Star selection!
        self.canvas.mpl_connect('button_press_event', self.button_press_callback)

        # Add first column  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

        # Make Group Box to put the label in, because it will misbehave
        instructionsGroupBox = QGroupBox('Instructions', self)
        vBox = QVBoxLayout()

        self.instructions = QLabel('''Click as near to the center of the star as possible. \
The first star that is choosen will be the <span style="background-color: #FFFF00">guide star</span>. \
All additional stars that are clicked on are the \
<span style="background-color: #FFA500">reference stars</span>. Star locations will appear in the list at \
right as they are selected; the guide star can be re-selected \
using the radio buttons. Errors will be shown in the output box below.''', self)
        # self.instructions.setSizePolicy(QSizePolicy.Preferred, QSizePolicy.Maximum)

        self.instructions.setWordWrap(True)

        vBox.addWidget(self.instructions)
        instructionsGroupBox.setLayout(vBox)
        instructionsGroupBox.setSizePolicy(QSizePolicy.Preferred, QSizePolicy.Maximum)
        mainGrid.addWidget(instructionsGroupBox, 0, 2, 1, 2, alignment=QtCore.Qt.AlignTop)

        # Log to update
        self.log_textbox = QTextEdit(self)
        self.log_textbox.setPlaceholderText('No stars selected.')
        self.log_textbox.setSizePolicy(QSizePolicy.Preferred, QSizePolicy.MinimumExpanding)
        self.log_textbox.setMinimumSize(10, 500)
        mainGrid.addWidget(self.log_textbox, 1, 2, alignment=QtCore.Qt.AlignTop)

        mainGrid.setRowMinimumHeight(1, self.image_dim*.5)


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

        mainGrid.addWidget(axGroupBox, 2, 2, 2, 1)


        # Add second column  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

        starsGroupBox = QGroupBox(self)
        col2Form = QFormLayout()

        # Show selected stars
        col2Form.addRow(QLabel('Guide\nStar?', self), QLabel('Star Position', self))
        nStars = 15

        self.star_positions = []
        self.guidestar_buttons = []
        for n in range(nStars):
            # Create radio buttons
            guidestar_button = QRadioButton(starsGroupBox)
            if n == 0: guidestar_button.setChecked(True)
            guidestar_button.toggled.connect(self.make_setGuideStar(guidestar_button))
            if n == 0: guidestar_button.setEnabled(False)
            if n!= 0: guidestar_button.setCheckable(False)

            # Create star position textboxes
            star_position = QLineEdit(starsGroupBox)
            if n == 0: star_position.setPlaceholderText('No stars selected')

            col2Form.addRow(guidestar_button, star_position)

        self.star_positions = [child for child in starsGroupBox.findChildren(QLineEdit)]
        self.guidestar_buttons = [child for child in starsGroupBox.findChildren(QRadioButton)]

        col2Form.setLabelAlignment(QtCore.Qt.AlignHCenter)

        starsGroupBox.setLayout(col2Form)
        mainGrid.addWidget(starsGroupBox, 1, 3, 1, 1)

        # Add "Done" button
        self.done_button = QPushButton('Done', self)
        self.done_button.setToolTip('Close the window')
        self.done_button.clicked.connect(self.fileQuit)
        self.done_button.setMinimumSize(150, 50)
        mainGrid.addWidget(self.done_button, 2, 3, 2, 1, alignment = QtCore.Qt.AlignVCenter)



    @pyqtSlot()
    def update_axes(self):
        '''Changes the axes limits of the matplotlib canvas to the current values of the axis limit textboxes
        '''
        self.canvas.axes.set_xlim(float(self.x1_textbox.text()), float(self.x2_textbox.text()))
        self.canvas.axes.set_ylim(float(self.y2_textbox.text()), float(self.y1_textbox.text()))
        self.canvas.fig.subplots_adjust(top=.95, left = .07)
        self.canvas.draw()

    def update_textboxes(self):
        '''Changes the values of the axis limit textboxes to the current axis limits of the matplotlib canvas'''
        self.x1_textbox.setText(str(self.canvas.axes.get_xlim()[0]))
        self.x2_textbox.setText(str(self.canvas.axes.get_xlim()[1]))
        self.y1_textbox.setText(str(self.canvas.axes.get_ylim()[1]))
        self.y2_textbox.setText(str(self.canvas.axes.get_ylim()[0]))

    def move_in_axes(self, event):
        '''Updates the cursor position textbox when the cursor moves within the matplotlib axis'''
        if event.inaxes:
            # self.cursor_textbox.setEnabled(True)
            self.cursor_textbox.setText('({:.0f}, {:.0f})'.format(event.xdata, event.ydata))


    def button_press_callback(self, event):
        '''WHenever a mouse button is pressed, determines if a button press is within the matplotlib axis;
        if so calls get_ind_under_point'''
        if not event.inaxes: return
        if event.button != 1: return
        self._ind = self.get_ind_under_point(event)

    def get_ind_under_point(self, event):
        '''If user clicks within one epsilon of an identified star, appends that star's position to inds;
        notifies user if no star is nearby or if selected star is already selected'''

        d = np.sqrt((self.x - event.xdata)**2 + (self.y - event.ydata)**2)
        i = np.unravel_index(np.nanargmin(d), d.shape)
        indseq = np.nonzero(np.equal(d, np.amin(d)))[0]
        ind = indseq[0]

        if d[ind] >= self.epsilon:
            if self.printOutput: print('No star within {} pixels. No star selected.'.format(self.epsilon))

            redText = "<br/><span style=\" color:#ff0000;\" >"
            redText += 'No star within {} pixels. No star selected.'.format(self.epsilon)
            redText += "</span>"

            self.log_textbox.setHtml(self.log_textbox.toHtml() + redText)

        elif ind in self.inds:
            if self.printOutput: print('Star already selected, please choose another star')

            redText = "<br/><span style=\" color:#ff0000;\" >"
            redText += 'Star already selected, please choose another star'
            redText += "</span>"

            self.log_textbox.setHtml(self.log_textbox.toHtml() + redText)

        else:
            if self.printOutput: print('Star selected: x={:.1f}, y={:.1f}'.format(self.x[ind],self.y[ind]))
            self.log_textbox.setHtml(self.log_textbox.toHtml() + '<br/>Star selected: x={:.1f}, y={:.1f}'.format(self.x[ind],self.y[ind]))
            self.star_positions[len(self.inds)].setText('x={:.1f}, y={:.1f}'.format(self.x[ind],self.y[ind]))
            self.guidestar_buttons[0].setEnabled(True)
            self.guidestar_buttons[len(self.inds)].setCheckable(True)

            if len(self.inds) == 0:
                color = 'yellow'
            else:
                color = 'orange'
            self.canvas.axes.scatter(self.x[ind], self.y[ind], s=500, facecolors='none', edgecolors=color,  lw=2, marker='o')
            self.canvas.draw()

            self.inds.append(ind)

        return ind

    def make_setGuideStar(self, button):
        # Have to define second function inside to return a function to the toggle.connect method
        def setGuideStar():
            '''
            Change order of self.inds to reflect the new guide star
            Update GUI plotted circles
            '''

            # Plot current guide star as regular star
            self.canvas.axes.scatter(self.x[self.inds[0]], self.y[self.inds[0]], s=500, facecolors='none', edgecolors='orange',  lw=2, marker='o')

            # Update self.inds
            guide_ind = self.guidestar_buttons.index(button)
            self.inds.insert(0, self.inds.pop(guide_ind))

            # Plot new guide star as guide star
            self.canvas.axes.scatter(self.x[self.inds[0]], self.y[self.inds[0]], s=500, facecolors='none', edgecolors='yellow',  lw=2, marker='o')
            self.canvas.draw()

        return setGuideStar


    def fileQuit(self):
        '''Closes the application'''
        self.qApp.exit(0)  # Works only with self.close() after; same as qApp.quit()
        self.close()



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


def run_SelectStars(data, x, y, dist, printOutput=False):
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
    printOutput (boolean)
        Flag enabling output to the terminal

    Returns
    -------
    inds (list)
        List of indices of positions of selected stars
    '''

    # RUN GUI
    qApp = QApplication(sys.argv)
    aw = ApplicationWindow(data=data, x=x, y=y, dist=dist, qApp=qApp, printOutput=printOutput)
    aw.show()
    inds = aw.inds
    qApp.exec_() # Begin interactive session; pauses until qApp.exit() is called

    return inds


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

if __name__ == '__main__':
    data = fits.open('../data/LMCfootprint_80/LMCfootprint_80.31512449419834_-70.45278790693078.fits')[0].data
    data[data == 0] = 0.1 # Adjust so 0 shows up as such with a logNorm colorbar
    gauss_sigma = 5
    dist = 16

    inds = run_SelectStars(data, gauss_sigma, dist)
