"""Interactive GUI for adding background stars to an image.

Builds a GUI with PyQt5 that prompts the user to add background stars
to an image, either (1) randomly, (2) by defining the locations and
magnitudes of the stars to add, or (3) by querying the guide star catalog
(GSC).

Authors
-------
    - Lauren Chambers

Use
---
    This module can be used as such:
    ::
        from jwst_magic.convert_image.background_stars_GUI import BackgroundStarsWindow
        window = BackgroundStarsWindow(guider, fgs_mag, qApp=qApp,
                                       in_master_GUI=in_master_GUI)

    Required arguments:
        ``guider`` - guider number (1 or 2)
        ``fgs_mag`` - brightness of the guide star in FGS Magnitude

    Optional arguments:
        ``aApp`` - qApplication instance of parent GUI
        ``in_master_GUI`` - is this module being called as part of the master GUI?

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
import logging
import os
import random
import requests

# Third Party Imports
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.io import ascii as asc
import fgscountrate
import matplotlib as mpl
import numpy as np
from PyQt5 import uic
from PyQt5.QtWidgets import QDialog, QFileDialog, QTableWidgetItem, QMessageBox
import pysiaf

# Local Imports
from jwst_magic.star_selector.SelectStarsGUI import StarClickerMatplotlibCanvas

__location__ = os.path.realpath(os.path.join(os.getcwd(), os.path.dirname(__file__)))

# Start logger
LOGGER = logging.getLogger(__name__)


class BackgroundStarsWindow(QDialog):
    def __init__(self, guider, fgs_mag, ra, dec, qApp, in_master_GUI):
        """Defines attributes; calls initUI() method to set up user interface.
        """
        # Initialize general attributes
        self.qApp = qApp
        self.guider = guider
        self.fgs_mag = fgs_mag
        self.image_dim = 400
        self.in_master_GUI = in_master_GUI
        self.method = None
        self.extended = None
        self.x = []
        self.y = []
        self.fgs_mags = []

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

        # Populate ra and dec if present
        self.lineEdit_RA.setText(str(ra if ra is not None else ''))
        self.lineEdit_Dec.setText(str(dec if dec is not None else ''))

        # Create and load background stars dialog GUI session
        self.setWindowTitle("Add Background Stars")
        self.init_matplotlib()
        self.define_GUI_connections()
        self.show()

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # GUI CONSTRUCTION
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    def init_matplotlib(self):
        """Set up the matplotlib canvas that will preview the locations
        and magnitudes of the background stars relative to the guide
        star.
        """
        # Connect main matplotlib canvas and add to layout
        self.canvas = StarClickerMatplotlibCanvas(
            parent=self.frame_canvas, data=None, x=None, dpi=100,
            y=None, bottom=0.23, left=0.15)
        self.frame_canvas.layout().insertWidget(0, self.canvas)
        self.canvas.axes.set_xlim(0, 2048)
        self.canvas.axes.set_ylim(2048, 0)
        self.canvas.axes.set_xlabel('X [pixels]')
        self.canvas.axes.set_ylabel('Y [pixels]')

        # Plot guide star
        cmap = mpl.colors.Colormap('viridis_r')
        self.vmin, self.vmax = (self.fgs_mag + 8, self.fgs_mag - 1)
        self.guide_star = self.canvas.axes.scatter(1024, 1024, marker='*',
                                                   s=500, cmap=cmap,
                                                   vmin=self.vmax,
                                                   vmax=self.vmin)

        # Add colorbar
        self.canvas.cbar_ax = self.canvas.fig.add_axes([0.05, 0.1, 0.9, 0.03])
        norm = mpl.colors.Normalize(vmin=self.vmin, vmax=self.vmax)
        self.canvas.cbar = mpl.colorbar.ColorbarBase(self.canvas.cbar_ax,
                                                     norm=norm, cmap=cmap,
                                                     orientation='horizontal')
        self.canvas.cbar.ax.invert_xaxis()
        self.canvas.cbar.set_label('FGS Magnitude')

    def define_GUI_connections(self):
        # Main dialog widgets
        self.pushButton_done.clicked.connect(self.quit)
        self.pushButton_cancel.clicked.connect(self.quit)

        # Randomly added stars widgets
        self.groupBox_random.toggled.connect(self.on_check_section)
        self.pushButton_random.clicked.connect(self.draw_random_stars)

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
        """Raise a dialog box to load and parse a file listing
        background star positions and magnitudes; auto-populate the
        table and plot each star.

        Returns
        -------
        filename : str
            Name of the file cataloging background stars
        """
        filename, _ = QFileDialog.getOpenFileName(self,
                                                  'Open Background Stars Input File',
                                                  "",
                                                  "Input file (*.txt);;All files (*.*)")
        if filename:
            self.lineEdit_definedFile.setText(filename)

            # Parse the file
            tab = asc.read(filename)
            for i_row, row in enumerate(tab):
                if i_row + 1 > self.tableWidget.rowCount():
                    self.tableWidget.insertRow(i_row)
                for i_col, value in enumerate(row):
                    item = QTableWidgetItem(str(value))
                    self.tableWidget.setItem(i_row, i_col, item)
            self.draw_defined_stars()

            return filename

    def draw_random_stars(self):
        # Only draw new stars if all the needed parameters exist
        if self.lineEdit_magMin.text() == '' or \
                        self.lineEdit_magMax.text() == '' or \
                        self.lineEdit_nStars.text() == '':
            return

        # Randomly generate x, y, and fgs_mags
        fgs_mag = self.fgs_mag
        size = 2048
        nstars_random = int(self.lineEdit_nStars.text())
        vmin = fgs_mag + float(self.lineEdit_magMin.text())
        vmax = fgs_mag + float(self.lineEdit_magMax.text())
        self.x = random.sample(range(size), nstars_random)
        self.y = random.sample(range(size), nstars_random)
        self.fgs_mags = random.sample(
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
            self.x, self.y, c=self.fgs_mags, marker='*', s=500, cmap='viridis_r',
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
                    LOGGER.warning('Background Stars: There is a cell with non-numeric contents')
                    return

        # Alert user if the coordinates are out of bounds
        for i_row in range(self.tableWidget.rowCount()):
            if 0 > float(self.tableWidget.item(i_row, 0).text()) or \
                            float(self.tableWidget.item(i_row, 0).text()) > 2048:
                LOGGER.warning('Background Stars: X Location out of bounds.')
                return
            if 0 > float(self.tableWidget.item(i_row, 1).text()) or \
                            float(self.tableWidget.item(i_row, 1).text()) > 2048:
                LOGGER.warning('Background Stars: Y Location out of bounds.')
                return

        # Remove other stars and lines, if they have already been plotted
        self.clear_plot()

        # Read values from table
        self.x = []
        self.y = []
        self.fgs_mags = []
        for i_row in range(self.tableWidget.rowCount()):
            x = float(self.tableWidget.item(i_row, 0).text())
            y = float(self.tableWidget.item(i_row, 1).text())
            fgs_mag = float(self.tableWidget.item(i_row, 2).text())
            self.x.append(x)
            self.y.append(y)
            self.fgs_mags.append(fgs_mag)
            self.check_colorbar_limits(fgs_mag)

        # Plot every star
        self.defined_stars = self.canvas.axes.scatter(
            self.x, self.y, c=self.fgs_mags, marker='*', s=500, cmap='viridis_r',
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
        if self.lineEdit_RA.text() == '' or \
                        self.lineEdit_Dec.text() == '':
            return

        # Remove other stars and lines, if they have already been plotted
        self.clear_plot()

        # Read position angle from GUI
        position_angle = self.lineEdit_PA.text()
        if position_angle == '':
            no_PA_dialog = QMessageBox()
            no_PA_dialog.setText('No PA entered' + ' ' * 50)
            no_PA_dialog.setInformativeText(
                'It is not possible to place results of a GSC query onto the '
                'detector without specifying the position angle (roll angle).')
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
        mask = np.array([j is np.ma.masked for j in self.fgs_mags])
        LOGGER.info('Background Stars: Plotting {} stars onto GUIDER{} FOV.'
                    .format(len(self.x[~mask]), self.guider))

        # Check if the star magnitudes are outside the colorbar limits
        for fgs_mag in self.fgs_mags[~mask]:
            self.check_colorbar_limits(fgs_mag)

        # Plot stars with known fgs_mags
        self.catalog_stars = self.canvas.axes.scatter(
            self.x[~mask], self.y[~mask], c=self.fgs_mags[~mask], marker='*',
            s=500, cmap='viridis_r', vmin=self.vmax, vmax=self.vmin,
            label=None
        )
        # Plot stars with unknown fgs_mags
        if len(self.x[mask]) > 0:
            self.masked_catalog_stars = self.canvas.axes.scatter(
                self.x[mask], self.y[mask], c='white', marker='*', s=500,
                edgecolors='red', label='Unknown FGS Magnitude'
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

        # If the row is not empty, remove that row's x, y, fgs_mag from
        # list of parameters
        try:
            self.x.remove(float(self.tableWidget.item(i_row, 0).text()))
            self.y.remove(float(self.tableWidget.item(i_row, 1).text()))
            self.fgs_mags.remove(float(self.tableWidget.item(i_row, 2).text()))
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
        """Check if a certain magnitude is within the current colorbar
        limits. If not, extend the limits.
        """

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
        """Create and parse a web query to GSC 2.4.1 to determine the
        positions and magnitudes of objects around the guide star.

        References
        ----------
            For information about GSC 2:
                https://outerspace.stsci.edu/display/GC
        """
        # Parse RA and Dec
        ra_gs = coordinates.ra.degree
        dec_gs = coordinates.dec.degree
        
        # Set radius around ra and dec and query default (newest) GSC 
        radius = 1.6 / 60  # 1.6 arcmin in degrees
        df = fgscountrate.query_gsc(ra=ra_gs, dec=dec_gs, cone_radius=radius)
        LOGGER.info('Background Stars: Querying Newest Guide Star Catalog')

        # Only take the necessary columns
        queried_catalog = df[['ra', 'dec', 'classification', 'tmassJMag']]
        ra_list = df['ra'].values
        dec_list = df['dec'].values
        fgs_mag_list = []
        for i in range(len(df)):
            row = df.iloc[[i]]
            fgs = fgscountrate.FGSCountrate('', guider)
            try:
                _, _, fgs_magnitude, _ = fgs.query_fgs_countrate_magnitude(data_frame=row)
                fgs_mag_list.append(fgs_magnitude)
            except ValueError:
                fgs_mag_list.append(0)
        fgs_mag_list = np.array(fgs_mag_list)

        LOGGER.info('Background Stars: Finished query; found {} sources.'.format(len(ra_list)))

        # Only select sources that are stars
        mask_point_sources = [c == 0 for c in queried_catalog['classification']]
        ra_list = ra_list[mask_point_sources]
        dec_list = dec_list[mask_point_sources]
        fgs_mag_list = fgs_mag_list[mask_point_sources]

        # self.extended = queried_catalog[~mask_pointSources]

        # Remove the guide star!
        # (Assume that is the star closest to the center, if there is a star
        # within 1" of the pointing)
        distances = [np.sqrt((ra - ra_gs) ** 2 + (dec - dec_gs) ** 2) for (ra, dec) in zip(ra_list, dec_list)]
        i_mindist = np.where(min(distances) == distances)[0][0]  # Probably 0
        if distances[i_mindist] < 1 / 60 / 60:
            LOGGER.info(
                'Background Stars: Removing assumed guide star at {}, {}'.
                format(ra_list[i_mindist], dec_list[i_mindist]))
            mask_guidestar = [i != i_mindist for i in range(len(ra_list))]
            ra_list = ra_list[mask_guidestar]
            dec_list = dec_list[mask_guidestar]
            fgs_mag_list = fgs_mag_list[mask_guidestar]
        else:
            LOGGER.warning('Background Stars: No guide star found within 1 arcsec of the pointing.')

        # Convert RA/Dec (sky frame) to X/Y pixels (raw frame)
        siaf = pysiaf.Siaf('FGS')
        guider = siaf['FGS{}_FULL'.format(guider)]
        v2ref_arcsec = guider.V2Ref
        v3ref_arcsec = guider.V3Ref

        attitude_ref = pysiaf.utils.rotations.attitude(v2ref_arcsec, v3ref_arcsec, ra_gs, dec_gs, position_angle)
        v2, v3 = pysiaf.utils.rotations.getv2v3(attitude_ref, ra_list, dec_list)
        x_raw, y_raw = guider.tel_to_raw(v2, v3)

        # Only select the sources within the detector frame
        in_detector_frame = []
        for x, y in zip(x_raw, y_raw):
            if (0 < x < 2048) and (0 < y < 2048):
                in_detector_frame.append(True)
            else:
                in_detector_frame.append(False)

        self.x = x_raw[in_detector_frame]
        self.y = y_raw[in_detector_frame]
        self.fgs_mags = fgs_mag_list[in_detector_frame]

        LOGGER.info('Background Stars: Found {} sources in GUIDER{} FOV.'
                    .format(len(self.x), self.guider))

        return queried_catalog

    def quit(self):
        # If "cancel" was selected, don't save the data
        if self.sender() == self.pushButton_cancel:
            self.x = []
            self.y = []
            self.fgs_mags = []
            self.method = None

        # If not being called from the master GUI, exit the whole application
        if not self.in_master_GUI:
            self.qApp.exit(0)  # Works only with self.close() after; same as qApp.quit()

        # Close the star selector dialog window
        self.close()
