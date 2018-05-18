"""Add background stars to an FGS image by repeating the current image

Authors
-------
    - Lauren Chambers

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

# Standard Library Imports
import random
import logging
import requests
import sys
import os

# Third Party Imports
import numpy as np
import matplotlib as mpl
from astropy import units as u
from astropy.stats import sigma_clipped_stats
from astropy.coordinates import SkyCoord
from astropy.io import ascii as asc
from PyQt5 import QtCore, uic
from PyQt5.QtWidgets import QApplication, QDialog, QFileDialog, QTableWidgetItem, QMessageBox
import pysiaf

# Local Imports
from . import coordinate_transforms
from .convert_image import counts_to_jmag
from .star_selector.SelectStarsGUI import StarClickerMatplotlibCanvas

__location__ = os.path.realpath(os.path.join(os.getcwd(), os.path.dirname(__file__)))

# Start logger
LOGGER = logging.getLogger(__name__)


class BackgroundStarsWindow(QDialog):

    def __init__(self, guider, jmag, qApp, in_master_GUI):
        '''Defines attributes; calls initUI() method to set up user interface.
        '''
        # Initialize general attributes
        self.qApp = qApp
        self.guider = guider
        self.jmag = jmag
        self.image_dim = 400
        self.in_master_GUI = in_master_GUI
        self.method = None
        self.extended = None
        self.x = []
        self.y = []
        self.jmags = []

        # Initialize matplotlib plotting attributes
        self.cbar_vmin_line = None
        self.cbar_vmax_line = None
        self.random_stars = None
        self.defined_stars = None
        self.catalog_stars = None
        self.masked_catalog_stars = None
        self.legend = None

        # Initialize dialog object
        QDialog.__init__(self, modal=True)

        # Import .ui file
        uic.loadUi(os.path.join(__location__, 'background_stars.ui'), self)

        # Create and load background stars dialog GUI session
        self.setWindowTitle("Add Background Stars")
        self.init_matplotlib()
        self.define_GUI_connections()
        self.show()

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # GUI CONSTRUCTION
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    def init_matplotlib(self):
        '''Set up the matplotlib canvas that will preview the locations
        and magnitudes of the background stars relative to the guide
        star.
        '''
        # Connect main matplotlib canvas and add to layout
        self.canvas = StarClickerMatplotlibCanvas(
            parent=self.frame_canvas,  data=None, x=None, dpi=100,
            y=None, bottom=0.23, left=0.15)
        self.frame_canvas.layout().insertWidget(0, self.canvas)
        self.canvas.axes.set_xlim(0, 2048)
        self.canvas.axes.set_ylim(2048, 0)
        self.canvas.axes.set_xlabel('X [pixels]')
        self.canvas.axes.set_ylabel('Y [pixels]')

        # Plot guide star
        self.vmin, self.vmax = (self.jmag + 8, self.jmag - 1)
        self.guide_star = self.canvas.axes.scatter(1024, 1024, marker='*',
                                                   s=500, c=self.jmag,
                                                   cmap='viridis_r',
                                                   vmin=self.vmax,
                                                   vmax=self.vmin)

        # Add colorbar
        self.canvas.cbar_ax = self.canvas.fig.add_axes([0.05, 0.1, 0.9, 0.03])
        norm = mpl.colors.Normalize(vmin=self.vmin, vmax=self.vmax)
        self.canvas.cbar = mpl.colorbar.ColorbarBase(self.canvas.cbar_ax,
                                                     norm=norm, cmap='viridis_r',
                                                     orientation='horizontal')
        self.canvas.cbar.ax.invert_xaxis()
        self.canvas.cbar.set_label('J Magnitude')

    def define_GUI_connections(self):
        # Main dialog widgets
        self.pushButton_done.clicked.connect(self.quit)
        self.pushButton_cancel.clicked.connect(self.quit)

        # Randomly added stars widgets
        self.groupBox_random.toggled.connect(self.on_check_section)
        self.lineEdit_nStars.editingFinished.connect(self.draw_random_stars)
        self.lineEdit_magMin.editingFinished.connect(self.draw_random_stars)
        self.lineEdit_magMax.editingFinished.connect(self.draw_random_stars)

        # User-defined stars widgets
        self.groupBox_defined.toggled.connect(self.on_check_section)
        self.tableWidget.cellChanged.connect(self.draw_defined_stars)
        self.pushButton_addStar.clicked.connect(self.add_star)
        self.pushButton_deleteStar.clicked.connect(self.delete_star)
        self.pushButton_definedFile.clicked.connect(self.load_file)

        # Catalog query stars widgets
        self.groupBox_catalog.toggled.connect(self.on_check_section)
        self.pushButton_queryGSC.clicked.connect(self.draw_catalog_stars)

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # WIDGET CONNECTIONS
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    def on_check_section(self):
        sections = [self.groupBox_random, self.groupBox_catalog,
                    self.groupBox_defined]
        if self.sender().isChecked():
            for section in sections:
                if section != self.sender():
                    section.setChecked(False)

    def load_file(self):
        '''Raise a dialog box to load and parse a file listing
        background star positions and magnitudes; auto-populate the
        table and plot each star.

        Returns
        -------
        filename : str
            Name of the file cataloging background stars
        '''
        filename, _ = QFileDialog.getOpenFileName(self,
                                                 'Open Background Stars Input File',
                                                  "",
                                                  "Input file (*.txt);;All files (*.*)")
        if filename:
            self.lineEdit_definedFile.setText(filename)

            # Parse the file
            tab = asc.read(filename)
            print(tab)
            for i_row, row in enumerate(tab):
                if i_row + 1 > self.tableWidget.rowCount():
                    self.tableWidget.insertRow(i_row)
                for i_col, value in enumerate(row):
                    item = QTableWidgetItem(str(value))
                    print(i_row, i_col, value, type(value), item)
                    self.tableWidget.setItem(i_row, i_col, item)
            print('through the looop')

            self.draw_defined_stars()

            return filename

    def draw_random_stars(self):
        # Only draw new stars if all the needed parameters exist
        if self.lineEdit_magMin.text() == '' or\
           self.lineEdit_magMax.text() == '' or\
           self.lineEdit_nStars.text() == '':
            return

        # Randomly generate x, y, and jmags
        jmag = self.jmag
        size = 2048
        nstars_random = int(self.lineEdit_nStars.text())
        vmin = jmag + float(self.lineEdit_magMin.text())
        vmax = jmag + float(self.lineEdit_magMax.text())
        self.x = random.sample(range(size), nstars_random)
        self.y = random.sample(range(size), nstars_random)
        self.jmags = random.sample(
            set(np.linspace(vmin, vmax, 100)),
            nstars_random
            )

        # Check if the star magnitudes are outside the colorbar limits
        self.check_colorbar_limits(vmin)
        self.check_colorbar_limits(vmax)

        # Remove other stars and lines, if they have already been plotted
        self.clear_plot()

        # Add lines that show the boundaries of random stars
        cbar_vmin = 1 - (vmin - self.vmin) / (self.vmax - self.vmin)
        cbar_vmax = 1 - (vmax - self.vmin) / (self.vmax - self.vmin)
        self.cbar_vmin_line = self.canvas.cbar.ax.axvline(cbar_vmin, c='w')
        self.cbar_vmax_line = self.canvas.cbar.ax.axvline(cbar_vmax, c='w')

        # Plot every star
        self.random_stars = self.canvas.axes.scatter(
            self.x, self.y, c=self.jmags, marker='*', s=500, cmap='viridis_r',
            vmin=self.vmax, vmax=self.vmin
        )

        # Record what method was used
        self.method = "random"

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
                    LOGGER.warbubg('Background Stars: There is a cell with non-numeric contents')
                    return

        # Alert user if the coordinates are out of bounds
        for i_row in range(self.tableWidget.rowCount()):
            if 0 > float(self.tableWidget.item(i_row, 0).text()) or\
               float(self.tableWidget.item(i_row, 0).text()) > 2048:
                LOGGER.warning('Background Stars: X Location out of bounds.')
                return
            if 0 > float(self.tableWidget.item(i_row, 1).text()) or\
               float(self.tableWidget.item(i_row, 1).text()) > 2048:
                LOGGER.warning('Background Stars: Y Location out of bounds.')
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
            self.check_colorbar_limits(jmag)

        # Plot every star
        self.defined_stars = self.canvas.axes.scatter(
            self.x, self.y, c=self.jmags, marker='*', s=500, cmap='viridis_r',
            vmin=self.vmax, vmax=self.vmin
            )

        # Record what method was used
        self.method = "user-defined"

        # Redraw all necessary plot elements
        self.canvas.cbar.draw_all()
        self.canvas.cbar.ax.invert_xaxis()
        self.canvas.draw()

    def draw_catalog_stars(self):
        # Only draw new stars if all the needed parameters exist
        if self.lineEdit_RA.text() == '' or\
           self.lineEdit_Dec.text() == '':
            return

        # Remove other stars and lines, if they have already been plotted
        self.clear_plot()

        # Read position angle from GUI
        position_angle = self.lineEdit_PA.text()
        if position_angle == '':
            no_PA_dialog = QMessageBox()
            no_PA_dialog.setText('No PA entered' + ' ' * 50)
            no_PA_dialog.setInformativeText('It is not possible to place results of a GSC query onto the detector without specifying the position angle (roll angle).')
            no_PA_dialog.setStandardButtons(QMessageBox.Ok)
            no_PA_dialog.exec()
            return
        else:
            position_angle = float(position_angle)

        # Query guide star catalog (GSC) to find stars around given pointing
        # Convert from RA & Dec to pixel coordinates
        RAunit_index = int(self.comboBox_RAUnits.currentIndex())
        unit_RA = [None, u.hourangle, u.deg][RAunit_index]
        unit_Dec = u.deg
        coordinates = SkyCoord(self.lineEdit_RA.text(), self.lineEdit_Dec.text(),
                               unit=(unit_RA, unit_Dec))
        queried_catalog = self.query_gsc(coordinates, self.guider, position_angle)

        # Plot every star
        mask = np.array([j is np.ma.masked for j in self.jmags])
        LOGGER.info('Background Stars: Plotting {} stars onto GUIDER{} FOV.'.format(len(self.x[~mask]), self.guider))

        # Check if the star magnitudes are outside the colorbar limits
        for jmag in self.jmags[~mask]:
            self.check_colorbar_limits(jmag)

        # Plot stars with known jmags
        self.catalog_stars = self.canvas.axes.scatter(
            self.x[~mask], self.y[~mask], c=self.jmags[~mask], marker='*',
            s=500, cmap='viridis_r', vmin=self.vmax, vmax=self.vmin,
            label=None
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
        #         self.extended['x'], self.y[mask], c='white',
        #         marker='galaxy_icon.jpg', s=500,
        #         edgecolors='red', label='Extended source'
        #         )
        #     self.legend = self.canvas.axes.legend()

        # Record what method was used
        self.method = "catalog"

        # Redraw all necessary plot elements
        self.canvas.cbar.draw_all()
        self.canvas.cbar.ax.invert_xaxis()
        self.canvas.draw()

    def add_star(self):
        n_rows = self.tableWidget.rowCount()
        self.tableWidget.insertRow(n_rows)

    def delete_star(self):
        # Determine which row is highlighted
        i_row = self.tableWidget.currentRow()

        # If the row is not empty, remove that row's x, y, jmag from
        # list of parameters
        try:
            self.x.remove(float(self.tableWidget.item(i_row, 0).text()))
            self.y.remove(float(self.tableWidget.item(i_row, 1).text()))
            self.jmags.remove(float(self.tableWidget.item(i_row, 2).text()))
        except (ValueError, AttributeError) as e:
            pass

        # Remove that row from the table
        self.tableWidget.removeRow(i_row)

        # Redraw the defined stars
        self.draw_defined_stars()

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # HELPER FUNCTIONS
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    def check_colorbar_limits(self, j_magnitude):
        '''Check if a certain magnitude is within the current colorbar
        limits. If not, extend the limits.
        '''

        # If the new jmag_min is less than the current colorbar, extend
        if j_magnitude > self.vmin:
            self.vmin = j_magnitude + 1

        # If the new jmag_max is more than the current colorbar, extend
        if j_magnitude < self.vmax:
            self.vmax = j_magnitude - 1

        norm = mpl.colors.Normalize(vmin=self.vmin, vmax=self.vmax)
        self.canvas.cbar.set_clim(self.vmax, self.vmin)
        self.canvas.cbar.set_norm(norm)
        self.guide_star.set_clim(self.vmax, self.vmin)

        # Redraw the colorbar
        self.canvas.cbar.draw_all()
        self.canvas.cbar.ax.invert_xaxis()

    def clear_plot(self):
        # Remove stars and lines, if they have already been plotted
        names = ['random_stars', 'cbar_vmin_line', 'cbar_vmax_line',
                 'defined_stars', 'catalog_stars', 'masked_catalog_stars',
                 'legend']

        for name in names:
            if getattr(self, name) is not None:
                try:
                    getattr(self, name).remove()
                    setattr(self, name, None)
                except ValueError:
                    print('could not remove')

    def query_gsc(self, coordinates, guider, position_angle):
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
        radius = 1.6 / 60  # 1.6 arcmin in degrees
        web_query = "http://gsss.stsci.edu/webservices/vo/CatalogSearch.aspx?RA={:f}&DEC={:f}&DSN=+&FORMAT=CSV&CAT=GSC241&SR={:f}&".format(RA, Dec, radius)
        LOGGER.info('Background Stars: Querying GSC 2.4.1 at ' + web_query)
        page = requests.get(web_query)
        csv_data = page.text
        table = asc.read(csv_data)

        # Only take the necessary columns
        queried_catalog = table['ra', 'dec', 'classification', 'tmassJmag']
        RAs = table['ra']
        Decs = table['dec']
        jmag = table['tmassJmag']

        LOGGER.info('Background Stars: Finshed query; found {} sources.'.format(len(RAs)))

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
        i_mindist = np.where(min(distances))[0][0]  # Probably 0
        if distances[i_mindist] < 1/60/60:
            LOGGER.info('Background Stars: Removing assumed guide star at {}, {}'.format(RAs[i_mindist], Decs[i_mindist]))
            mask_guidestar = [i != i_mindist for i in range(len(RAs))]
            RAs = RAs[mask_guidestar]
            Decs = Decs[mask_guidestar]
            jmag = jmag[mask_guidestar]
        else:
            LOGGER.warning('Background Stars: No guide star found within 1 arcsec of the pointing.')

        # Convert RA/Dec (sky frame) to X/Y pixels (raw frame)
        siaf = pysiaf.Siaf('FGS')
        guider = siaf['FGS{}_FULL'.format(guider)]
        V2ref_arcsec = guider.V2Ref
        V3ref_arcsec = guider.V3Ref

        attitude_ref = pysiaf.utils.rotations.attitude(V2ref_arcsec, V3ref_arcsec, RA, Dec, position_angle)
        V2, V3 = pysiaf.utils.rotations.getv2v3(attitude_ref, RAs, Decs)
        x_det, y_det = guider.tel_to_det(V2, V3)
        x_raw, y_raw = y_det, x_det

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

        LOGGER.info('Background Stars: Found {} sources in GUIDER{} FOV.'.format(len(self.x), self.guider))

        return queried_catalog

    def quit(self):
        # If "cancel" was selected, don't save the data
        if self.sender() == self.pushButton_cancel:
            self.x = []
            self.y = []
            self.jmags = []
            self.method = None

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
    in_master_GUI = masterGUIapp is not None

    if in_master_GUI:
        qApp = masterGUIapp
    else:
        qApp = QtCore.QCoreApplication.instance()
        if qApp is None:
            qApp = QApplication(sys.argv)

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

    return stars, window.method


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
    if stars is True:
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
                             'positions and magnitudes. Ensure the same number '
                             'of entries is provided for all fields.')

        x_back = stars['x']
        y_back = stars['y']
        jmags_back = stars['jmag']

    else:
        raise TypeError('Unfamiliar value passed to bkgd_stars: {} Please pass boolean or dictionary of background star x, y, jmag.'.format(stars))

    # Add stars to image
    # Copy original data array
    add_data = np.copy(image)

    # (Try to) only use the data for added stars, not the noise
    mean, median, std = sigma_clipped_stats(image, sigma=0, iters=0)
    image[image < mean] = 0
    for x, y, jmag_back in zip(x_back, y_back, jmags_back):
        if not isinstance(jmag_back, np.ma.core.MaskedConstant):
            star_fgs_counts = counts_to_jmag.jmag_to_fgs_counts(jmag_back, guider)
            scale_factor = star_fgs_counts / fgs_counts

            star_data = image * scale_factor
            psfx = psfy = 2048

            # Crop the image to fit on the array at the specified location
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

            LOGGER.info('Background Stars: Adding background star with magnitude {:.1f} at location ({}, {}).'.format(jmag_back,
                                                                                              x, y))
            add_data[x1:x2, y1:y2] += star_data

    return add_data
