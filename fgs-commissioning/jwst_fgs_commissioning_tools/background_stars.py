# Add background stars to an FGS image by copying current image
import random
import logging
import requests

import numpy as np
from astropy.stats import sigma_clipped_stats
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.io import ascii as asc
from PyQt5 import QtCore
from PyQt5.QtWidgets import (QVBoxLayout, QApplication, QPushButton, QLineEdit,
                             QLabel, QTextEdit, QRadioButton, QGroupBox, QFrame,
                             QSizePolicy, QGridLayout, QDialog, QMessageBox,
                             QTableWidget, QComboBox)
from PyQt5.QtCore import pyqtSlot, QObject
import matplotlib as mpl
from matplotlib.cm import viridis_r
from astroquery.vizier import Vizier

from . import coordinate_transforms
from .convert_image import counts_to_jmag
from .star_selector.SelectStarsGUI import StarClickerMatplotlibCanvas

GROUPBOX_TITLE_STYLESHEET = 'QGroupBox { font-size: 18px; font-weight: bold; margin-top: 40px } QGroupBox::title { top: -30px }'

# Start logger
LOGGER = logging.getLogger(__name__)


class BackgroundStarsWindow(QDialog):

    def __init__(self, guider, jmag, qApp, in_master_GUI):
        '''Defines attributes; calls initUI() method to set up user interface.'''
        self.qApp = qApp
        self.guider = guider
        self.jmag = jmag
        self.image_dim = 400
        self.in_master_GUI = in_master_GUI

        self.extended = None

        # Initialize dialog object
        QDialog.__init__(self, modal=True)

        # Initialize user interface
        self.initUI()

    def initUI(self):
        '''Sets up interactive graphical user interface.
        '''

        # Set up star selector dialog window
        self.setWindowTitle("Add Background Stars")
        mainGrid = QGridLayout()  # set grid layout
        self.setLayout(mainGrid)
        self.setFocus()

        # Add plot - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        sc = StarClickerMatplotlibCanvas(self, width=5, height=4, dpi=100,
                         data=None, x=None, y=None, bottom=0.23, left=0.15)
        self.canvas = sc

        self.canvas.axes.set_xlim(0, 2048)
        self.canvas.axes.set_ylim(2048, 0)
        self.canvas.axes.set_xlabel('X [pixels]')
        self.canvas.axes.set_ylabel('Y [pixels]')

        # Plot guide star
        self.vmin, self.vmax = (self.jmag + 8, self.jmag - 1)
        self.guide_star = self.canvas.axes.scatter(1024, 1024, marker='*', s=500, c=self.jmag,
                                 cmap=viridis_r,  vmin=self.vmax, vmax=self.vmin)

        # Add colorbar
        self.canvas.cbar_ax = self.canvas.fig.add_axes([0.05, 0.1, 0.9, 0.03])
        norm = mpl.colors.Normalize(vmin=self.vmin, vmax=self.vmax)
        self.canvas.cbar = mpl.colorbar.ColorbarBase(self.canvas.cbar_ax,
                                                     norm=norm, cmap=viridis_r,
                                                     orientation='horizontal')
        self.canvas.cbar.ax.invert_xaxis()
        self.canvas.cbar.set_label('J Magnitude')

        self.canvas.setMinimumSize(self.image_dim, self.image_dim)
        mainGrid.addWidget(sc, 0, 1, 3, 1)

        # Add other widgets  - - - - - - - - - - - - - - - - - - - - - - - - - -
        # Build the random stars section
        randomGroupBox = self.create_random_section()
        mainGrid.addWidget(randomGroupBox, 0, 0, 1, 1)

        # Build the defined stars section
        definedGroupBox = self.create_defined_section()
        mainGrid.addWidget(definedGroupBox, 1, 0, 1, 1)

        # Build the random stars section
        catalogGroupBox = self.create_catalog_section()
        mainGrid.addWidget(catalogGroupBox, 2, 0, 1, 1)

        # Add done button
        button_done = QPushButton("Done", self)
        button_done.clicked.connect(self.quit)
        button_done.setMinimumSize(150, 50)
        mainGrid.addWidget(button_done, 3, 0, 1, 2,
                           alignment=QtCore.Qt.AlignHCenter)

        # Show GUI - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        self.show()

    def create_random_section(self):
        '''
        '''
        # Set up group box
        self.randomGroupBox = QGroupBox('Add Stars Randomly', self)
        self.randomGroupBox.setStyleSheet(GROUPBOX_TITLE_STYLESHEET)
        self.randomGroupBox.setCheckable(True)
        randomGrid = QGridLayout()
        self.randomGroupBox.setLayout(randomGrid)
        self.randomGroupBox.setChecked(True)
        self.randomGroupBox.toggled.connect(self.on_check_section)

        # Add number of stars textbox
        randomGrid.addWidget(QLabel('Number of Stars: ', self), 0, 0, 1, 2)
        self.textbox_n_stars = QLineEdit('5', self)
        self.textbox_n_stars.setMaximumSize(30, 20)
        self.textbox_n_stars.editingFinished.connect(self.draw_random_stars)
        randomGrid.addWidget(self.textbox_n_stars, 0, 2)

        # Add magnitude range texboxes
        # randomGrid.addWidget(QLabel('Magnitudes Dimmer than Guide Star: ', self), 1, 0, 1, 5)
        randomGrid.addWidget(QLabel('From', self), 2, 0)
        self.textbox_jmag_min = QLineEdit('7', self)
        self.textbox_jmag_min.setMaximumSize(30, 20)
        self.textbox_jmag_min.editingFinished.connect(self.draw_random_stars)
        randomGrid.addWidget(self.textbox_jmag_min, 2, 1)
        randomGrid.addWidget(QLabel('mag dimmer to', self), 2, 2)
        self.textbox_jmag_max = QLineEdit('3', self)
        self.textbox_jmag_max.setMaximumSize(30, 20)
        self.textbox_jmag_min.editingFinished.connect(self.draw_random_stars)
        randomGrid.addWidget(self.textbox_jmag_max, 2, 3)
        randomGrid.addWidget(QLabel('mag dimmer than the guide star', self), 2, 4)

        self.draw_random_stars()


        return self.randomGroupBox

    def create_defined_section(self):
        '''
        '''
        # Set up group box
        self.definedGroupBox = QGroupBox('Define Stars to Add', self)
        self.definedGroupBox.setStyleSheet(GROUPBOX_TITLE_STYLESHEET)
        self.definedGroupBox.setCheckable(True)
        definedGrid = QGridLayout()
        self.definedGroupBox.setLayout(definedGrid)
        self.definedGroupBox.setChecked(False)
        self.definedGroupBox.toggled.connect(self.on_check_section)

        # Add table of parameters
        self.tableWidget = QTableWidget(self)
        self.tableWidget.setColumnCount(3)
        self.tableWidget.setRowCount(1)
        self.tableWidget.setHorizontalHeaderLabels(['x [pixels]', 'y [pixels]', 'Jmag'])
        self.tableWidget.setMaximumSize(320, 500)
        definedGrid.addWidget(self.tableWidget, 0, 0)

        self.tableWidget.cellChanged.connect(self.draw_defined_stars)

        # Add button to add another star
        self.button_add_star = QPushButton('Add Another Star', self)
        self.button_add_star.clicked.connect(self.add_star)
        definedGrid.addWidget(self.button_add_star, 1, 0)

        # Add button to deleted star
        self.button_add_star = QPushButton('Delete Star', self)
        self.button_add_star.clicked.connect(self.delete_star)
        definedGrid.addWidget(self.button_add_star, 2, 0)

        return self.definedGroupBox

    def create_catalog_section(self):
        '''
        '''
        # Set up group box
        self.catalogGroupBox = QGroupBox('Add Stars from Guide Star Catalog 2.4.1', self)
        self.catalogGroupBox.setStyleSheet(GROUPBOX_TITLE_STYLESHEET)
        self.catalogGroupBox.setCheckable(True)
        catalogGrid = QGridLayout()
        self.catalogGroupBox.setLayout(catalogGrid)
        self.catalogGroupBox.setChecked(False)
        self.catalogGroupBox.toggled.connect(self.on_check_section)

        # Add RA and Dec text boxes
        catalogGrid.addWidget(QLabel('RA: ', self), 0, 0, alignment=QtCore.Qt.AlignRight)
        self.textbox_RA = QLineEdit(self)
        self.textbox_RA.setMaximumSize(150, 20)
        catalogGrid.addWidget(self.textbox_RA, 0, 1)

        catalogGrid.addWidget(QLabel('Dec: ', self), 1, 0, alignment=QtCore.Qt.AlignRight)
        self.textbox_Dec = QLineEdit(self)
        self.textbox_Dec.setMaximumSize(150, 20)
        catalogGrid.addWidget(self.textbox_Dec, 1, 1)

        # Add units
        self.cb_RAUnits = QComboBox(self)
        self.cb_RAUnits.addItem("-Select-")
        self.cb_RAUnits.addItem("Hours")
        self.cb_RAUnits.addItem("Degrees")
        catalogGrid.addWidget(self.cb_RAUnits, 0, 2)
        catalogGrid.addWidget(QLabel('Degrees', self), 1, 2)

        # Add "query" button
        self.button_query = QPushButton('Query GSC', self)
        self.button_query.clicked.connect(self.draw_catalog_stars)
        catalogGrid.addWidget(self.button_query, 0, 3, 2, 1)

        return self.catalogGroupBox

    def on_check_section(self):
        sections = [self.randomGroupBox, self.catalogGroupBox, self.definedGroupBox]
        if self.sender().isChecked():
            for section in sections:
                if section != self.sender():
                    section.setChecked(False)

    def draw_random_stars(self):
        # Only draw new stars if all the needed parameters exist
        if self.textbox_jmag_min.text() == '' or\
           self.textbox_jmag_max.text() == '' or\
           self.textbox_n_stars.text() == '':
            return

        # Randomly generate x, y, and jmags
        jmag = self.jmag
        size = 2048
        nstars_random = int(self.textbox_n_stars.text())
        vmin = jmag + float(self.textbox_jmag_min.text())
        vmax = jmag + float(self.textbox_jmag_max.text())
        self.x = random.sample(range(size), nstars_random)
        self.y = random.sample(range(size), nstars_random)
        self.jmags = random.sample(
            set(np.linspace(vmin, vmax, 100)),
            nstars_random
            )

        # If the new jmag_min is less than the current colorbar, extend
        # print(self.canvas.cbar.get_clim(), self.vmax, self.vmin)
        if vmin > self.vmin:
            self.vmin = vmin + 1
            norm = mpl.colors.Normalize(vmin=self.vmin, vmax=self.vmax)
            self.canvas.cbar.set_clim(self.vmax, self.vmin)
            self.canvas.cbar.set_norm(norm)


            self.guide_star.set_clim(self.vmax, self.vmin)

        # Remove other stars and lines, if they have already been plotted
        self.clear_plot()

        # Add lines that show the boundaries of random stars
        cbar_vmin = 1 - (vmin - self.vmin) / (self.vmax - self.vmin)
        cbar_vmax = 1 - (vmax - self.vmin) / (self.vmax - self.vmin)
        self.cbar_vmin_line = self.canvas.cbar.ax.axvline(cbar_vmin, c='w')
        self.cbar_vmax_line = self.canvas.cbar.ax.axvline(cbar_vmax, c='w')


        # Plot every star
        self.random_stars = self.canvas.axes.scatter(
            self.x, self.y, c=self.jmags, marker='*', s=500, cmap=viridis_r,
            vmin=self.vmax, vmax=self.vmin
            )

        # Redraw all necessary plot elements
        self.canvas.cbar.draw_all()
        self.canvas.cbar.ax.invert_xaxis()
        self.canvas.draw()

    def draw_defined_stars(self):
        # Only draw stars if the table is full and numeric
        for i_row in range(self.tableWidget.rowCount()):
            for i_col in range(self.tableWidget.columnCount()):
                if not self.tableWidget.item(i_row, i_col):
                    return
                elif self.tableWidget.item(i_row, i_col).text() == '':
                    return
                elif not self.tableWidget.item(i_row, i_col).text().isnumeric():
                    print('There is a cell with non-numeric contents')
                    return

        # Alert user if the coordinates are out of bounds
        for i_row in range(self.tableWidget.rowCount()):
            if 0 > float(self.tableWidget.item(i_row, 0).text()) or\
               float(self.tableWidget.item(i_row, 0).text()) > 2048:
                print('X Location out of bounds.')
                return
            if 0 > float(self.tableWidget.item(i_row, 1).text()) or\
               float(self.tableWidget.item(i_row, 1).text())> 2048:
                print('Y Location out of bounds.')
                return

        # Remove other stars and lines, if they have already been plotted
        self.clear_plot()

        # Read values from table
        self.x = []
        self.y = []
        self.jmags = []
        for i_row in range(self.tableWidget.rowCount()):
            x = float(self.tableWidget.item(i_row, 0).text())
            y = float(self.tableWidget.item(i_row, 1).text())
            jmag = float(self.tableWidget.item(i_row, 2).text())
            self.x.append(x)
            self.y.append(y)
            self.jmags.append(jmag)

        # Plot every star
        self.defined_stars = self.canvas.axes.scatter(
            self.x, self.y, c=self.jmags, marker='*', s=500, cmap=viridis_r,
            vmin=self.vmax, vmax=self.vmin
            )

        # Redraw all necessary plot elements
        self.canvas.cbar.draw_all()
        self.canvas.cbar.ax.invert_xaxis()
        self.canvas.draw()

    def draw_catalog_stars(self):
        # Only draw new stars if all the needed parameters exist
        if self.textbox_RA.text() == '' or\
           self.textbox_Dec.text() == '':
            return

        # Remove other stars and lines, if they have already been plotted
        self.clear_plot()

        # Query guide star catalog (GSC) to find stars around given pointing
        # Convert from RA & Dec to pixel coordinates
        RAunit_index = int(self.cb_RAUnits.currentIndex())
        unit_RA = [None, u.hourangle, u.deg][RAunit_index]
        unit_Dec = u.deg
        coordinates = SkyCoord(self.textbox_RA.text(), self.textbox_Dec.text(), unit=(unit_RA, unit_Dec))
        queried_catalog = self.query_gsc(coordinates, self.guider)

        # Plot every star
        mask = np.array([j is np.ma.masked for j in self.jmags])
        # Plot stars with known jmags
        self.catalog_stars = self.canvas.axes.scatter(
            self.x[~mask], self.y[~mask], c=self.jmags[~mask], marker='*',
            s=500, cmap=viridis_r, vmin=self.vmax, vmax=self.vmin,
            label = None
            )
        # Plot stars with unknown jmags
        if len(self.x[mask]) > 0:
            self.masked_catalog_stars = self.canvas.axes.scatter(
                self.x[mask], self.y[mask], c='white', marker='*', s=500,
                edgecolors='red', label='Unknown J Magnitude'
                )
            self.legend = self.canvas.axes.legend()

        # # Plot extended sources
        # if self.extended:
        #     self.masked_catalog_stars = self.canvas.axes.scatter(
        #         self.extended['x'], self.y[mask], c='white', marker='galaxy_icon.jpg', s=500,
        #         edgecolors='red', label='Extended source'
        #         )
        #     self.legend = self.canvas.axes.legend()

        # Redraw all necessary plot elements
        self.canvas.cbar.draw_all()
        self.canvas.cbar.ax.invert_xaxis()
        self.canvas.draw()


    def clear_plot(self):
        # Remove random stars and lines, if they have already been plotted
        try:
            self.random_stars.remove()
            self.cbar_vmin_line.remove()
            self.cbar_vmax_line.remove()
            self.random_stars = None
        except (AttributeError, ValueError) as e:
            pass

        # Remove defined stars, if they have already been plotted
        try:
            self.defined_stars.remove()
        except (AttributeError, ValueError) as e:
            pass

        # Remove catalog stars, if they have already been plotted
        try:
            self.catalog_stars.remove()
        except (AttributeError, ValueError) as e:
            pass
        # Remove masekd catalog stars, if they have already been plotted
        try:
            self.masked_catalog_stars.remove()
            self.legend.remove()
        except (AttributeError, ValueError) as e:
            pass

    def add_star(self):
        n_rows = self.tableWidget.rowCount()
        self.tableWidget.insertRow(n_rows)

    def delete_star(self):
        # Determine which row is highlighted
        i_row = self.tableWidget.selectedItems()[0].row()

        # Remove that row's x, y, jmag from list of parameters
        self.x.remove(float(self.tableWidget.item(i_row, 0).text()))
        self.y.remove(float(self.tableWidget.item(i_row, 1).text()))
        self.jmags.remove(float(self.tableWidget.item(i_row, 2).text()))

        # Remove that row from the table
        self.tableWidget.removeRow(i_row)

        # Redraw the defined stars
        self.draw_defined_stars()

    def query_gsc(self, coordinates, guider):
        '''Create and parse a web query to GSC 2.4.1 to determine the
        positions and magnitudes of objects around the guide star.

        References
        ----------
            For information about GSC 2:
                https://outerspace.stsci.edu/display/GC
        '''
        # Parse RA and Dec
        RA = coordinates.ra.degree
        Dec = coordinates.dec.degree

        # Query MAST to get GSC 2.4.1 results in CSV form
        radius = 1.6 / 60 # 1.6 arcmin in degrees
        web_query = "http://gsss.stsci.edu/webservices/vo/CatalogSearch.aspx?RA={:f}&DEC={:f}&DSN=+&FORMAT=CSV&CAT=GSC241&SR={:f}&".format(RA, Dec, radius)
        print(web_query)
        page = requests.get(web_query)
        csv_data = page.text
        table = asc.read(csv_data)

        # Only take the necessary columns
        queried_catalog = table['ra', 'dec', 'classification', 'tmassJmag']
        RAs = table['ra']
        Decs = table['dec']
        jmag = table['tmassJmag']

        print('Finshed query; found {} sources.'.format(len(RAs)))

        # Only select sources that are stars
        mask_pointSources = [c == 0 for c in queried_catalog['classification']]
        RAs = RAs[mask_pointSources]
        Decs = Decs[mask_pointSources]
        jmag = jmag[mask_pointSources]

        # self.extended = queried_catalog[~mask_pointSources]

        # Remove the guide star!
        # (Assume that is the star closest to the center, if there is a star
        # within 1" of the pointing)
        distances = [np.sqrt((ra - RA)**2 + (dec - Dec)**2) for (ra, dec) in zip(RAs, Decs)]
        i_mindist = np.where(min(distances))[0][0] # probably just zero
        print(distances)
        print(i_mindist, distances[i_mindist])
        if distances[i_mindist] < 1/60/60:
            print('Removing assumed guide star at {}, {}'.format(RAs[i_mindist], Decs[i_mindist]))
            np.delete(RAs, i_mindist)
            np.delete(Decs, i_mindist)
            np.delete(jmag, i_mindist)
        else:
            print('No guide star found within 1 arcsec of the pointing.')

        # Convert RA/Dec to X/Y pixels
        V2 = (RAs - RA) *60 * 60
        V3 = (Decs - Dec) *60 * 60
        x_dhas, y_dhas = coordinate_transforms.Idl2DHAS(-V2, V3)
        x_raw, y_raw = coordinate_transforms.DHAS2Raw(x_dhas, y_dhas, guider)

        # Only select the sources within the detector frame
        in_detector_frame = []
        for x, y in zip(x_raw, y_raw):
            if (0 < x < 2048) and (0 < y < 2048):
                in_detector_frame.append(True)
            else:
                in_detector_frame.append(False)


        self.x = x_raw[in_detector_frame]
        self.y = y_raw[in_detector_frame]
        self.jmags = jmag[in_detector_frame]




        print('Found {} sources in detector FOV.'.format(len(self.x)))

        return queried_catalog

    def quit(self):
        # If not being called from the master GUI, exit the whole application
        if not self.in_master_GUI:
            self.qApp.exit(0)  # Works only with self.close() after; same as qApp.quit()

        # Close the star selector dialog window
        self.close()


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# EXTERNAL FUNCTIONS
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


def run_background_stars_GUI(guider, jmag, masterGUIapp=None):
    # RUN GUI
    if masterGUIapp:
        qApp = masterGUIapp
        in_master_GUI = True
    else:
        qApp = QtCore.QCoreApplication.instance()
        if qApp is None:
            qApp = QApplication(sys.argv)
        in_master_GUI = False

    window = BackgroundStarsWindow(guider, jmag, qApp=qApp, in_master_GUI=in_master_GUI)

    if masterGUIapp:
        window.exec_()  # Begin interactive session; pauses until window.exit() is called
    else:
        qApp.exec_()

    # Create dictionary to pass to ``add_background_stars``
    stars = {}
    stars['x'] = window.x
    stars['y'] = window.y
    stars['jmag'] = window.jmags

    return stars


def add_background_stars(image, stars, jmag, fgs_counts, guider):
    """Add artificial copies of the input PSF to the image to mimic
    background stars.

    Arguments
    ---------
    image: numpy array
        The input data with the original star
    stars
        Either a boolean or a dictionary that defines how to add the
        additional stars to the image
    jmag: float
        The desired J magnitude of the original star
    fgs_counts: float
        The desired FGS counts in the original star
    guider: int
        The number of the guider used to take the data

    Returns
    -------
    add_data: numpy array
        Image array with original star and additional stars
    """
    # (Try to) only use the data, not the noise
    mean, median, std = sigma_clipped_stats(image, sigma=0, iters=0)
    image[image < mean] = 0

    # Randomly create 5 locations on the image
    size = 2048
    nstars_random = 5

    # Determine jmag and fgs_counts of guide star
    if not fgs_counts:
        if not jmag:
            LOGGER.warning('No counts or J magnitude given, setting to default')
            jmag = 11
        fgs_counts = counts_to_jmag.jmag_to_fgs_counts(jmag, guider)
    else:
        jmag = counts_to_jmag.fgs_counts_to_jmag(fgs_counts, guider)

    # If the flag is simply set to "True", randomly place 5 stars on the image
    if stars == True:
        x_back = random.sample(range(size), nstars_random)
        y_back = random.sample(range(size), nstars_random)
        # Create the new stars 5 mags or more dimmer
        jmags_back = random.sample(set(np.linspace(jmag + 7, jmag + 4, 100)), nstars_random)

    # If users passed a dictionary to the bkgd_stars argument, add stars
    # according the dictionary
    elif type(stars) == dict:
        input_lengths = [len(stars[key]) for key in stars.keys()]
        if len(set(input_lengths)) != 1:
            raise ValueError('Invalid dictionary provided for background star '
                'positions and magnitudes. Ensure the same number of entries is '
                'provided for all fields.')

        x_back = stars['x']
        y_back = stars['y']
        jmags_back = stars['jmag']

    else:
        raise TypeError('Unfamiliar value passed to bkgd_stars: {} Please pass boolean or dictionary of background star x, y, jmag.'.format(stars))

    # Add stars to image
    add_data = np.copy(image)
    for x, y, jmag in zip(x_back, y_back, jmags_back):
        star_fgs_counts = counts_to_jmag.jmag_to_fgs_counts(jmag, guider)
        scale_factor = star_fgs_counts / fgs_counts

        star_data = image * scale_factor
        psfx = psfy = 2048

        x1 = max(0, int(x) - int(psfx / 2))
        x2 = min(2048, int(x) + int(psfx / 2) + 1)
        y1 = max(0, int(y) - int(psfy / 2))
        y2 = min(2048, int(y) + int(psfy / 2) + 1)

        if x > 1024:
            star_data = star_data[:x2 - x1]
        else:
            star_data = star_data[2048 - (x2 - x1):]
        if y > 1024:
            star_data = star_data[:, :y2 - y1]
        else:
            star_data = star_data[:, 2048 - (y2 - y1):]

        # print('Adding background star with magnitude {:.1f} at location ({}, {}).'.format(jmag,
        #                                                                                   x, y))
        add_data[x1:x2, y1:y2] += star_data

    return add_data
