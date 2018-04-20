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

# Third Party
import numpy as np
import matplotlib
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure
from matplotlib.colors import LogNorm
import matplotlib.pyplot as plt
from PyQt5 import QtCore
from PyQt5.QtWidgets import (QVBoxLayout, QApplication, QPushButton, QLineEdit,
                             QLabel, QTextEdit, QRadioButton, QGroupBox, QFrame,
                             QSizePolicy, QGridLayout, QDialog, QMessageBox,
                             QListView, QToolButton, QCheckBox, QComboBox)
from PyQt5.QtGui import QStandardItemModel, QStandardItem
from PyQt5.QtCore import pyqtSlot

matplotlib.rcParams['font.family'] = 'serif'
matplotlib.rcParams['font.weight'] = 'light'
matplotlib.rcParams['mathtext.bf'] = 'serif:normal'

from ..star_selector.SelectStarsGUI import StarClickerMatplotlibCanvas

GROUPBOX_TITLE_STYLESHEET = '\
QGroupBox { \
    font-size: 14px; \
    font-weight: bold; \
    margin-top: 30px \
} \
QGroupBox::title { \
    top: -20px \
}'

class SegmentGuidingWindow(QDialog):
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
                 print_output=False):
        '''Defines attributes; calls initUI() method to set up user interface.'''
        self.qApp = qApp
        self.print_output = print_output

        self.image_dim = 800

        self.data = data
        self.x = x
        self.y = y

        self._ind = None
        self.inds = []
        self.inds_of_inds = []
        self.epsilon = dist

        self.in_master_GUI = in_master_GUI
        self.n_orientations = 0

        # Initialize dialog object
        QDialog.__init__(self, modal=True)

        # Initialize user interface
        self.initUI()

    def initUI(self):
        '''Sets up interactive graphical user interface.
        '''

        # Set up star selector dialog window
        self.setWindowTitle("FGS Segment Guide Star Override")
        mainGrid = QGridLayout()  # set grid layout
        self.setLayout(mainGrid)
        self.setFocus()

        # Add plot - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        sc = StarClickerMatplotlibCanvas(self, width=5, height=4, dpi=100,
                         data=self.data, x=self.x, y=self.y, left=0.1, right=0.9)
        self.canvas = sc
        self.canvas.compute_initial_figure(self.canvas.fig, self.data, self.x,
                                           self.y)
        self.canvas.zoom_to_crop()
        self.canvas.setMinimumSize(self.image_dim, 0)
        mainGrid.addWidget(self.canvas, 0, 0, 6, 1)

        # Add cursor-tracking
        cursorFrame = self.create_cursor_section()
        mainGrid.addWidget(cursorFrame, 5, 0, 1, 1, alignment=QtCore.Qt.AlignHCenter)

        # Update cursor position and pixel value under cursor
        self.canvas.mpl_connect('motion_notify_event',
                                self.update_cursor_position)

        # Star selection!
        self.canvas.mpl_connect('button_press_event',
                                self.button_press_callback)

        # Add other widgets  - - - - - - - - - - - - - - - - - - - - - - - - - -

        # Make Group Box to put the label in, because it will misbehave
        instructionsGroupBox = QGroupBox('Instructions', self)
        vBox = QVBoxLayout()

        self.instructions = QLabel('''<b>Guide and Reference Star Selection:</b> Click as near to the center of the star \
        as possible. The first star that is choosen will be the <span style=\
        "background-color: #FFFF00">guide star</span>. All additional stars that are \
        clicked on are the <span style="background-color: #FFA500">reference stars\
        </span>. Star locations will appear in the list below at left as they are selected; \
        the guide star can be re-selected using the radio buttons. Errors will be shown\
        in the output box below.<br><br>
        <b>Segment Guiding Override Commands:</b> Once you are happy with your selection, \
        click "Save Command" to save the stars you have selected as a \
        segment guiding override command. You will see it appear in the list below at right. If you want to select multiple \
        orientations to guide on, repeat the star selection process and again \
        click "Save Command" when you are done. Segment override commands \
        can be deleted, and the order in which they will be commanded can be \
        altered using the buttons in the "Segment Override Commands" box. <br> The \
        center of the segment override pointing (i.e. the location that will \
        match the guide star RA and Dec) can be aligned with either the center \
        of the segment array or with a specific segment by setting \
        the "Center of Override Pointing" parameters.''', self)
        # self.instructions.setSizePolicy(QSizePolicy.Preferred, QSizePolicy.Maximum)

        self.instructions.setWordWrap(True)
        instructionsGroupBox.setStyleSheet(GROUPBOX_TITLE_STYLESHEET)

        vBox.addWidget(self.instructions)
        instructionsGroupBox.setLayout(vBox)
        instructionsGroupBox.setSizePolicy(QSizePolicy.Preferred,
                                           QSizePolicy.Maximum)
        mainGrid.addWidget(instructionsGroupBox, 0, 1, 1, 4,
                           alignment=QtCore.Qt.AlignTop)

        # Log to update
        self.log_textbox = QTextEdit(self)
        self.log_textbox.setPlaceholderText('No stars selected.')
        self.log_textbox.setMinimumSize(350, 50)
        self.log_textbox.setMaximumSize(700, 150)
        mainGrid.addWidget(self.log_textbox, 1, 4, 1, 1,
                           alignment=QtCore.Qt.AlignHCenter)

        # Plot slice of profile under cursor
        prof = StarClickerMatplotlibCanvas(self, width=4, height=3, dpi=100,
                                           data=self.data, left=0.2, right=0.95,
                                           bottom=0.2)
        prof.init_profile()
        self.profile = prof
        self.canvas.mpl_connect('motion_notify_event', self.update_profile)
        self.profile.setMinimumSize(200, 250)
        mainGrid.addWidget(self.profile, 1, 1, 2, 2,
                           alignment=QtCore.Qt.AlignVCenter)

        # Create colorbar-updating section
        cbarGroupBox = self.create_colorbar_section(sc)
        cbarGroupBox.setMaximumSize(200, 150)
        mainGrid.addWidget(cbarGroupBox, 3, 1, 2, 1,
                           alignment=QtCore.Qt.AlignTop|QtCore.Qt.AlignHCenter)

        # Create axes-updating section
        axGroupBox = self.create_axes_section(sc)
        axGroupBox.setMinimumSize(250, 150)
        mainGrid.addWidget(axGroupBox, 3, 2, 2, 1,
                           alignment=QtCore.Qt.AlignTop|QtCore.Qt.AlignHCenter)

        # Create selected stars section
        starsGroupBox = self.create_selectedstars_section()
        starsGroupBox.setMaximumSize(500, 550)
        mainGrid.addWidget(starsGroupBox, 1, 3, 4, 1, alignment=QtCore.Qt.AlignTop)

        # Create segment guiding section
        segmentGuidingGroupBox = self.create_segmentguiding_section()
        segmentGuidingGroupBox.setMaximumSize(500, 225)
        mainGrid.addWidget(segmentGuidingGroupBox, 2, 4, 2, 1,
                           alignment=QtCore.Qt.AlignTop)

        # Create override pointing center section
        overrideCenterGroupBox = self.create_overridecenter_section()
        overrideCenterGroupBox.setMaximumSize(500, 120)
        mainGrid.addWidget(overrideCenterGroupBox, 4, 4, 1, 1,
                           alignment=QtCore.Qt.AlignTop)

        # Buttons to close the GUI
        closeButtonsFrame = self.create_closebuttons_section()
        mainGrid.addWidget(closeButtonsFrame, 5, 1, 1, 4)

        # Show GUI - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        self.show()

    def create_cursor_section(self):
        cursorFrame = QFrame(self)
        cursorGrid = QGridLayout()
        cursorFrame.setLayout(cursorGrid)
        cursorFrame.setMaximumSize(600, 60)
        cursorFrame.setFrameStyle(QFrame.NoFrame)

        # Show current cursor position
        self.cursor_label = QLabel('Cursor Position:', self,
                                   alignment=QtCore.Qt.AlignRight)
        cursorGrid.addWidget(self.cursor_label, 0, 0,
                             alignment=QtCore.Qt.AlignVCenter)

        self.cursor_textbox = QLineEdit(self)
        self.cursor_textbox.setPlaceholderText('Move cursor into axes')
        self.cursor_textbox.setFixedSize(150, 20)
        cursorGrid.addWidget(self.cursor_textbox, 0, 1)

        self.canvas.mpl_connect('motion_notify_event',
                                self.update_cursor_position)

        # Show value under cursor
        self.pixel_label = QLabel('Pixel Value:', self,
                                  alignment=QtCore.Qt.AlignRight)
        cursorGrid.addWidget(self.pixel_label, 0, 2,
                             alignment=QtCore.Qt.AlignVCenter)

        self.pixel_textbox = QLineEdit(self)
        self.pixel_textbox.setPlaceholderText('Move cursor into axes')
        self.pixel_textbox.setFixedSize(150, 20)
        cursorGrid.addWidget(self.pixel_textbox, 0, 3)

        return cursorFrame

    def create_axes_section(self, sc):
        # Create axis-updating section
        axGroupBox = QGroupBox('Axes Limits', self)
        axGrid = QGridLayout()
        axGroupBox.setLayout(axGrid)
        axGroupBox.setStyleSheet(GROUPBOX_TITLE_STYLESHEET)

        axGrid.addWidget(QLabel('X: ( ', self), 0, 0, 1, 1)
        axGrid.addWidget(QLabel(',', self), 0, 2, 1, 1)
        axGrid.addWidget(QLabel(' )', self), 0, 4, 1, 1)

        axGrid.addWidget(QLabel('Y: ( ', self), 1, 0, 1, 1)
        axGrid.addWidget(QLabel(',', self), 1, 2, 1, 1)
        axGrid.addWidget(QLabel(' )', self), 1, 4, 1, 1)

        self.x1_textbox = QLineEdit(str(sc.axes.get_xlim()[0]), self)
        self.x1_textbox.setFixedSize(80, 20)
        axGrid.addWidget(self.x1_textbox, 0, 1, 1, 1)

        self.x2_textbox = QLineEdit(str(sc.axes.get_xlim()[1]), self)
        self.x2_textbox.setFixedSize(80, 20)
        axGrid.addWidget(self.x2_textbox, 0, 3, 1, 1)

        self.y1_textbox = QLineEdit(str(sc.axes.get_ylim()[1]), self)
        self.y1_textbox.setFixedSize(80, 20)
        axGrid.addWidget(self.y1_textbox, 1, 1, 1, 1)

        self.y2_textbox = QLineEdit(str(sc.axes.get_ylim()[0]), self)
        self.y2_textbox.setFixedSize(80, 20)
        axGrid.addWidget(self.y2_textbox, 1, 3, 1, 1)

        self.axlims_button = QPushButton('Update axis limits', self)
        self.axlims_button.clicked.connect(self.update_axes)
        self.axlims_button.setFixedSize(150, 25)
        axGrid.addWidget(self.axlims_button, 2, 0, 1, 5,
                         alignment=QtCore.Qt.AlignHCenter)

        # Add "Zoom Fit" button
        self.zoom_button = QPushButton('Zoom Fit', self)
        self.zoom_button.setToolTip('Zoom to encompass all data')
        self.zoom_button.clicked.connect(sc.zoom_to_fit)
        self.zoom_button.clicked.connect(self.update_textboxes)
        self.zoom_button.setFixedSize(150, 25)
        axGrid.addWidget(self.zoom_button, 3, 0, 1, 5,
                         alignment=QtCore.Qt.AlignHCenter)

        # Add "Crop to Data" button
        self.crop_button = QPushButton('Crop to Data', self)
        self.crop_button.setToolTip('Zoom crop to data')
        self.crop_button.clicked.connect(sc.zoom_to_crop)
        self.crop_button.clicked.connect(self.update_textboxes)
        self.crop_button.setFixedSize(150, 25)
        axGrid.addWidget(self.crop_button, 4, 0, 1, 5,
                         alignment=QtCore.Qt.AlignHCenter)

        return axGroupBox

    def create_colorbar_section(self, sc):
        # Create colorbar-updating section
        cbarGroupBox = QGroupBox('Colorbar Limits', self)
        cbarGrid = QGridLayout()
        cbarGroupBox.setLayout(cbarGrid)
        cbarGroupBox.setStyleSheet(GROUPBOX_TITLE_STYLESHEET)

        cbarGrid.addWidget(QLabel('Min Value:  ', self), 0, 0)
        cbarGrid.addWidget(QLabel('Max Value:  ', self), 1, 0)

        self.vmin_textbox = QLineEdit(str(sc.fitsplot.get_clim()[0]), self)
        self.vmin_textbox.setFixedSize(80, 20)
        cbarGrid.addWidget(self.vmin_textbox, 0, 1)

        self.vmax_textbox = QLineEdit(str(sc.fitsplot.get_clim()[1]), self)
        self.vmax_textbox.setFixedSize(80, 20)
        cbarGrid.addWidget(self.vmax_textbox, 1, 1)

        self.cbarlims_button = QPushButton('Update colorbar limits', self)
        self.cbarlims_button.setStyleSheet('QPushButton {color: black}')
        self.cbarlims_button.clicked.connect(self.update_cbar)
        self.cbarlims_button.setFixedSize(170, 25)
        cbarGrid.addWidget(self.cbarlims_button, 2, 0, 1, 2,
                           alignment=QtCore.Qt.AlignHCenter)

        return cbarGroupBox

    def create_selectedstars_section(self):
        starsGroupBox = QGroupBox('Selected Guide and Reference Stars', self)
        selectStarsGrid = QGridLayout()
        starsGroupBox.setLayout(selectStarsGrid)
        starsGroupBox.setStyleSheet(GROUPBOX_TITLE_STYLESHEET)

        # Show selected stars
        guide_label = QLabel('Guide\nStar?', self)
        guide_label.setFixedSize(40, 30)
        selectStarsGrid.addWidget(guide_label, 0, 0,
                                  alignment=QtCore.Qt.AlignVCenter)

        star_label = QLabel('Star Position', self)
        star_label.setFixedSize(100, 30)
        selectStarsGrid.addWidget(star_label, 0, 1)

        self.nStars = 11
        for n in range(self.nStars):
            # Create radio buttons
            guidestar_button = QRadioButton(starsGroupBox)
            if n == 0: guidestar_button.setChecked(True)
            guidestar_button.released.connect(self.make_setguidestar(guidestar_button))
            if n == 0: guidestar_button.setEnabled(False)
            if n != 0: guidestar_button.setCheckable(False)

            # Create star position textboxes
            star_position = QLineEdit(starsGroupBox)
            star_position.setMinimumSize(170, 20)
            if n == 0: star_position.setPlaceholderText('No stars selected')

            # Create 'remove' buttons
            remove_button = QPushButton('Delete', starsGroupBox)
            remove_button.clicked.connect(self.make_removestar(remove_button))
            remove_button.installEventFilter(self)
            remove_button.setEnabled(False)

            # Add to starsGroupBox grid layout
            selectStarsGrid.addWidget(guidestar_button, n + 1, 0)
            selectStarsGrid.addWidget(star_position, n + 1, 1)
            selectStarsGrid.addWidget(remove_button, n + 1, 2)

        self.star_positions = [child for child in starsGroupBox.findChildren(QLineEdit)]
        self.guidestar_buttons = [child for child in starsGroupBox.findChildren(QRadioButton)]
        self.remove_buttons = [child for child in starsGroupBox.findChildren(QPushButton)]

        save_button = QPushButton('Save Command', starsGroupBox)
        save_button.setFixedSize(200, 25)
        save_button.clicked.connect(self.save_orientation_to_list)
        selectStarsGrid.addWidget(save_button, 12, 0, 1, 3,
                                  alignment=QtCore.Qt.AlignHCenter)

        clear_button = QPushButton('Clear Selection', starsGroupBox)
        clear_button.setFixedSize(200, 25)
        clear_button.clicked.connect(self.clear_selected_stars)
        selectStarsGrid.addWidget(clear_button, 13, 0, 1, 3,
                                  alignment=QtCore.Qt.AlignHCenter)

        return starsGroupBox

    def create_overridecenter_section(self):
        overrideCenterGroupBox = QGroupBox('Center of Override Pointing', self)
        overrideCenterGrid = QGridLayout()
        overrideCenterGroupBox.setLayout(overrideCenterGrid)
        overrideCenterGroupBox.setStyleSheet(GROUPBOX_TITLE_STYLESHEET)

        # Use calculated center of array as pointing
        self.mean_position_checkbox = QCheckBox('Use Center of Segment Array', self)
        self.mean_position_checkbox.toggled.connect(self.update_center_avg)
        overrideCenterGrid.addWidget(self.mean_position_checkbox, 0, 0, 1, 2)

        # Use location of specific array as pointing
        segment_center_label = QLabel('Use Segment as Center:', self)
        overrideCenterGrid.addWidget(segment_center_label, 1, 0)
        self.segment_center_dropdown = QComboBox(self)
        self.segment_center_dropdown.addItem('-Select Segment-')
        for i in range(18):
            self.segment_center_dropdown.addItem(str(i + 1))
        self.segment_center_dropdown.activated.connect(self.update_center_seg)
        overrideCenterGrid.addWidget(self.segment_center_dropdown, 2, 0)

        # Set to default: segNum = 0
        self.mean_position_checkbox.setChecked(True)

        return overrideCenterGroupBox

    def create_segmentguiding_section(self):
        segmentGuidingGroupBox = QGroupBox('Segment Guiding Override Commands', self)
        segmentGuidingGrid = QGridLayout()
        segmentGuidingGroupBox.setLayout(segmentGuidingGrid)
        # segmentGuidingGroupBox.setSizePolicy(QSizePolicy.Preferred, QSizePolicy.Maximum)
        segmentGuidingGroupBox.setStyleSheet(GROUPBOX_TITLE_STYLESHEET)

        # Show segment override commands
        self.segment_override_list = QListView(segmentGuidingGroupBox)
        self.segment_override_model = QStandardItemModel(self.segment_override_list)
        # self.segment_override_list.setMaximumSize(500, 200)
        self.segment_override_list.setModel(self.segment_override_model)
        segmentGuidingGrid.addWidget(self.segment_override_list, 0, 0, 2, 2)

        # Add default text
        item = QStandardItem('No commands saved')
        self.segment_override_model.appendRow(item)
        self.segment_override_list.setEnabled(False)

        # Load button
        load_button = QPushButton('Load', segmentGuidingGroupBox)
        load_button.clicked.connect(self.load_orientation)
        segmentGuidingGrid.addWidget(load_button, 2, 0)

        # Delete button
        remove_button = QPushButton('Delete', segmentGuidingGroupBox)
        remove_button.clicked.connect(self.delete_orientation)
        segmentGuidingGrid.addWidget(remove_button, 2, 1)

        # Up and down arrows
        move_up_button = QToolButton(segmentGuidingGroupBox)
        move_up_button.setArrowType(QtCore.Qt.UpArrow)
        move_up_button.clicked.connect(self.move_orientation_up)
        segmentGuidingGrid.addWidget(move_up_button, 0, 2,
                                     alignment=QtCore.Qt.AlignBottom)

        move_down_button = QToolButton(segmentGuidingGroupBox)
        move_down_button.setArrowType(QtCore.Qt.DownArrow)
        move_down_button.clicked.connect(self.move_orientation_down)
        segmentGuidingGrid.addWidget(move_down_button, 1, 2,
                                     alignment=QtCore.Qt.AlignTop)

        return segmentGuidingGroupBox

    def create_closebuttons_section(self):
        closeButtonsFrame = QFrame(self)
        closeButtonsGrid = QGridLayout()
        closeButtonsFrame.setLayout(closeButtonsGrid)
        closeButtonsFrame.setFrameStyle(QFrame.NoFrame)

        # Add "Done" button
        self.done_button = QPushButton('Done', self)
        self.done_button.setToolTip('Close the window')
        self.done_button.clicked.connect(self.fileQuit)
        self.done_button.setMinimumSize(150, 50)
        closeButtonsGrid.addWidget(self.done_button, 0, 0,
                                   alignment=QtCore.Qt.AlignBottom)

        # Add "Cancel" button
        self.cancel_button = QPushButton('Cancel', self)
        self.cancel_button.setToolTip('Close the window and discard changes')
        self.cancel_button.clicked.connect(self.cancel)
        self.cancel_button.setMinimumSize(150, 50)
        closeButtonsGrid.addWidget(self.cancel_button, 1, 0,
                                   alignment=QtCore.Qt.AlignTop)

        return closeButtonsFrame


    @pyqtSlot()
    def eventFilter(self, remove_button, event):
        '''Parse out the cursor moving in or out of a star deletion button, and
        updating the matplotlib axis with a red highlighted circle accordingly'''
        if event.type() == QtCore.QEvent.Enter:
            if remove_button.isEnabled():
                # Determine index of star corresponding to button
                star_ind = self.remove_buttons.index(remove_button)
                ind_of_star_ind = self.inds_of_inds.index(star_ind)

                # Re-draw selected star as red
                self.canvas.axes.lines[ind_of_star_ind].set_markeredgecolor('red')

                self.canvas.draw()
            return True

        if event.type() == QtCore.QEvent.Leave:
            if remove_button.isEnabled():
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
            self.pixel_textbox.setText('{:.0f}'.format(self.data[int(event.ydata),
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

            countrate = np.sum(self.data[int(event.ydata) - 1:int(event.ydata) + 2,
                                         int(event.xdata) - 1:int(event.xdata) + 2])
            # self.profile.axes.text(0.02, 0.9, '3 x 3 countrate: {:.1f}'.format(countrate),
            #                        transform=self.profile.axes.transAxes)

            # print(dir(self.profile.countrate_label))
            self.profile.countrate_label.set_text('3 x 3 countrate: {:.0f}'.format(countrate))

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
            self.star_positions[len(self.inds)].setText('{} (x={:.1f}, y={:.1f})'.format(ind + 1,
                                                                                         self.x[ind],
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

    def save_orientation_to_list(self):
        if not self.segment_override_list.isEnabled():
            self.segment_override_list.setEnabled(True)
            self.segment_override_model.takeRow(0)

        orientation_summary = ', '.join([str(int(position.text()[:2])) for position in self.star_positions if position.text() != ''])
        item = QStandardItem(orientation_summary)
        self.segment_override_model.appendRow(item)
        # self.segment_override_list.setModel(self.segment_override_model)
        self.clear_selected_stars()
        self.n_orientations += 1

    def clear_selected_stars(self):
        # Reset text boxes and selection buttons
        for i in range(self.nStars):
            self.star_positions[i].setText('')
            if i == 0: self.star_positions[i].setPlaceholderText('No stars selected')
            self.guidestar_buttons[i].setEnabled(False)
            if i == 0: self.guidestar_buttons[i].setChecked(True)
            else: self.guidestar_buttons[i].setCheckable(False)
            self.remove_buttons[i].setEnabled(False)

        # Reset matplotlib canvas
        ind_to_remove = 0
        for i in range(len(self.canvas.axes.lines)):
            # But don't remove the center marker
            if hasattr(self, 'center'):
                if self.canvas.axes.lines[ind_to_remove] == self.center[0]:
                    ind_to_remove = 1
                    continue
            self.canvas.axes.lines.remove(self.canvas.axes.lines[ind_to_remove])
        self.canvas.draw()

        # Reset index list
        self.inds = []

    def load_orientation(self):
        self.clear_selected_stars()

        selected_orientation_index = self.segment_override_list.selectedIndexes()[0].row()
        selected_orientation = self.segment_override_model.item(selected_orientation_index).text()

        selected_indices = [int(s) for s in selected_orientation.split(', ')]

        for i, ind in enumerate(selected_indices):
            # Load text boxes and selection buttons
            self.star_positions[i].setText('{} (x={:.1f}, y={:.1f})'.format(ind,
                                                                            self.x[ind-1],
                                                                            self.y[ind-1]))
            self.guidestar_buttons[i].setEnabled(True)
            self.guidestar_buttons[i].setCheckable(True)
            self.remove_buttons[i].setEnabled(True)

            # Load matplotlib canvas
            if i == 0:
                c = 'yellow'
            else:
                c = 'darkorange'
            self.canvas.axes.plot(self.x[ind-1], self.y[ind-1], 'o', ms=25,
                                  mfc='none', mec=c, mew=2, lw=0)

        self.canvas.draw()

        # Load index list
        self.inds = selected_indices

    def delete_orientation(self):
        selected_orientation_index = self.segment_override_list.selectedIndexes()[0].row()
        self.segment_override_model.takeRow(selected_orientation_index)
        self.n_orientations -= 1

    def move_orientation_up(self):
        # Move the item up
        selected_orientation_index = self.segment_override_list.selectedIndexes()[0].row()
        row = self.segment_override_model.takeRow(selected_orientation_index)
        # Don't move if at the end of the list
        new_index = max(0, selected_orientation_index - 1)
        self.segment_override_model.insertRow(new_index, row)

        # Maintain which item is selected
        ind = self.segment_override_model.index(new_index, 0)
        self.segment_override_list.setCurrentIndex(ind)

    def move_orientation_down(self):
        # Move the item down
        selected_orientation_index = self.segment_override_list.selectedIndexes()[0].row()
        row = self.segment_override_model.takeRow(selected_orientation_index)
        # Don't move if at the end of the list
        new_index = min(self.segment_override_model.rowCount(), selected_orientation_index + 1)
        self.segment_override_model.insertRow(new_index, row)

        # Maintain which item is selected
        ind = self.segment_override_model.index(new_index, 0)
        self.segment_override_list.setCurrentIndex(ind)

    def update_center_seg(self):
        '''Use the location of a specific segment as the pointing center
        '''
        # Uncheck "use segment center" box
        self.mean_position_checkbox.setChecked(False)

        # Remove old center
        try:
            self.canvas.axes.lines.remove(self.center[0])
        except (AttributeError, ValueError) as e:
            pass
        self.segNum = ''

        if self.segment_center_dropdown.currentText() != "-Select Segment-":
            # Replace with new center
            i_seg_center = int(self.segment_center_dropdown.currentText()) - 1
            self.center = self.canvas.axes.plot(self.x[i_seg_center],
                                                self.y[i_seg_center], 'x', ms=20,
                                                alpha=0.8, mfc='white',
                                                mec='white', mew=5, lw=0)

            self.segNum = int(self.segment_center_dropdown.currentText())

        self.canvas.draw()

    def update_center_avg(self):
        '''Use the mean of the array as the pointing center
        '''
        # Remove old center
        try:
            self.canvas.axes.lines.remove(self.center[0])
        except (AttributeError, ValueError) as e:
            pass
        self.segNum = ''

        if self.mean_position_checkbox.isChecked():
            self.segment_center_dropdown.setCurrentIndex(0)

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
        try:
            int(self.segNum)
        except (ValueError, AttributeError) as e:
            no_segNum_selected_dialog = QMessageBox()
            no_segNum_selected_dialog.setText('No center segment number')
            no_segNum_selected_dialog.setInformativeText('The center of override pointing (segNum) has not been defined. Please define before quitting.')
            no_segNum_selected_dialog.setStandardButtons(QMessageBox.Ok)
            no_segNum_selected_dialog.exec()
            return

        self.answer = True
        # If the user didn't choose any stars, ask if they really want to quit.
        if self.inds == [] and self.n_orientations == 0:
            no_stars_selected_dialog = QMessageBox()
            no_stars_selected_dialog.setText('No stars selected' + ' ' * 50)
            no_stars_selected_dialog.setInformativeText('The tool will not be able to continue. Do you want to quit anyway?')
            no_stars_selected_dialog.setStandardButtons(QMessageBox.Yes | QMessageBox.No)
            no_stars_selected_dialog.buttonClicked.connect(self.nostars_dialog)
            no_stars_selected_dialog.exec()

        # If they do, then quit.
        if self.answer:
            # If not being called from the master GUI, exit the whole application
            if not self.in_master_GUI:
                self.qApp.exit(0)  # Works only with self.close() after; same as qApp.quit()

            # Close the star selector dialog window
            self.close()

    def nostars_dialog(self, button):
        if 'No' in button.text():
            self.answer = False

    def cancel(self):
        '''Closes the star selector window and clears indices'''

        # Clear the indices (i.e. don't save user selections)
        self.inds = []

        # If not being called from the master GUI, exit the whole application
        if not self.in_master_GUI:
            self.qApp.exit(0)  # Works only with self.close() after; same as qApp.quit()

        # Close the star selector dialog window
        self.close()

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

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

    # Bring window to front
    try:
        plt.get_current_fig_manager().window.raise_()
    except AttributeError:
        pass

    # Begin interactive session; pauses until window.exit() is called
    if masterGUIapp:
        window.exec_()
    else:
        qApp.exec_()

    # Save indices of selected stars to pass
    inds = []
    for i in range(window.n_orientations):
        orientation = window.segment_override_model.item(i).text()
        selected_indices = [int(s) for s in orientation.split(', ')]
        inds.append(selected_indices)

    # Save index of center segment (pointing)
    segNum = window.segNum

    return inds, segNum
