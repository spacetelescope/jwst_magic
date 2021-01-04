"""Collection of unit tests to verify the correct function of the background
stars module.

Authors
-------
    - Lauren Chambers
    - Shannon Osborne

Use
---
    ::
        pytest test_background_stars.py
"""
import os
import sys

from astropy.io import fits
JENKINS = 'jenkins' in os.getcwd()
if not JENKINS:
    from PyQt5 import QtCore
    from PyQt5.QtWidgets import QTableWidgetItem, QDialogButtonBox, QApplication
import pytest

# Local Imports
from ..utils import utils
from ..convert_image.background_stars import add_background_stars

if not JENKINS:
    from ..convert_image.background_stars_GUI import BackgroundStarsDialog
    from ..masterGUI import MasterGui

SOGS = utils.on_sogs_network()
if not SOGS:
    from pytestqt import qtbot

ROOT = 'test_background_stars'
import pathlib
__location__ = str(pathlib.Path(__file__).parent.absolute())
NIRCAM_IM = os.path.join(__location__, 'data', 'nircam_data_1_ga.fits')


@pytest.fixture()
def master_gui():
    """Set up QApplication object for the Master GUI"""
    #global app
    app = QApplication.instance()
    if app is None:
        app = QApplication(sys.argv)

    master_gui = MasterGui(root=ROOT, in_file=None, out_dir=__location__,
                           segment_guiding=True, app=app, itm=False)

    return master_gui

@pytest.fixture()
def bkgdstars_dialog():
    """Set up QApplication object for the Master GUI"""
    # Set defaults
    guider = 1
    fgs_mag = 12

    #global app
    app = QApplication.instance()
    if app is None:
        app = QApplication(sys.argv)

    bkgdstars_dialog = BackgroundStarsDialog(guider, fgs_mag, ra=None, dec=None,
                                             in_master_GUI=False)
    return bkgdstars_dialog


def test_add_background_stars():
    """Basic test to make sure add_background_stars works
    """
    with fits.open(NIRCAM_IM) as hdulist:
        data = hdulist[1].data

    stars = {'x': [1500], 'y': [200], 'fgs_mag': [10.0]}

    add_background_stars(data, stars, 11.0, 'FGS Magnitude', 1)


@pytest.mark.skipif(JENKINS, reason="Can't import PyQt5 on Jenkins server.")
def test_init_background_stars(bkgdstars_dialog):
    """Make sure the background_stars GUI can launch without errors.
    """
    # Schedule press of "Cancel" button
    QtCore.QTimer.singleShot(0, bkgdstars_dialog.buttonBox.button(QDialogButtonBox.Cancel).clicked)


@pytest.mark.skipif(JENKINS, reason="Can't import PyQt5 on Jenkins server.")
@pytest.mark.skipif(SOGS, reason="Can't import pytest-qt on SOGS machine.")
def test_run_background_stars_gui_populated(qtbot, master_gui):
    """Test that the RA and DEC information from APT get passed through and
    correctly populated in the background stars GUI.
    """
    # Initialize main window
    qtbot.addWidget(master_gui)

    # Set general input
    qtbot.keyClicks(master_gui.lineEdit_inputImage, NIRCAM_IM)
    qtbot.mouseClick(master_gui.buttonGroup_name.buttons()[1], QtCore.Qt.LeftButton)  # set naming method
    qtbot.mouseClick(master_gui.buttonGroup_guider.buttons()[1], QtCore.Qt.LeftButton)  # set to guider 1
    qtbot.keyClicks(master_gui.lineEdit_root, ROOT)  # set root
    qtbot.keyClicks(master_gui.textEdit_out, __location__)  # set out directory
    qtbot.mouseClick(master_gui.buttonGroup_guider.buttons()[0], QtCore.Qt.LeftButton)
    qtbot.mouseClick(master_gui.buttonGroup_guider.buttons()[1], QtCore.Qt.LeftButton)
    qtbot.keyClicks(master_gui.lineEdit_manualid, '1148')
    qtbot.keyClicks(master_gui.lineEdit_manualobs, '01')
    qtbot.mouseClick(master_gui.pushButton_manualid, QtCore.Qt.LeftButton)

    # Check the RA and DEC have been populated, but PA has not
    def handle_dialog():
        try:
            assert master_gui._bkgdstars_dialog.lineEdit_RA.text() != ''
            assert master_gui._bkgdstars_dialog.lineEdit_Dec.text() != ''
            assert master_gui._bkgdstars_dialog.lineEdit_PA.text() == ''
            qtbot.mouseClick(master_gui._bkgdstars_dialog.buttonBox.button(QDialogButtonBox.Ok), QtCore.Qt.LeftButton)
        except AssertionError:
            # If something raising an error above, need to close the pop up gui anyway
            qtbot.mouseClick(master_gui._bkgdstars_dialog.buttonBox.button(QDialogButtonBox.Ok), QtCore.Qt.LeftButton)

    with qtbot.capture_exceptions() as exceptions:
        QtCore.QTimer.singleShot(500, handle_dialog)
        qtbot.mouseClick(master_gui.pushButton_backgroundStars, QtCore.Qt.LeftButton, delay=1)

    # An incomplete entry pressing done will error
    expected_err = 'No background stars selected'
    assert expected_err in str(exceptions[0][1]), "Wrong error captured. Caught: '{}', Expected: '{}'".format(
        str(exceptions[0][1]), expected_err)


@pytest.mark.skipif(JENKINS, reason="Can't import PyQt5 on Jenkins server.")
def test_run_background_stars_gui_query(bkgdstars_dialog):
    """Test running the GUI works when using the GSC query method"""

    # Set parameters (will have ~4 stars in FOV)
    bkgdstars_dialog.lineEdit_RA.setText('1')
    bkgdstars_dialog.lineEdit_Dec.setText('5')
    bkgdstars_dialog.lineEdit_PA.setText('5')

    # Query and then press OK
    bkgdstars_dialog.pushButton_queryGSC.click()
    bkgdstars_dialog.buttonBox.button(QDialogButtonBox.Ok).click()

    # Get parameters
    assert len(bkgdstars_dialog.x) != 0
    assert len(bkgdstars_dialog.y) != 0
    assert len(list(bkgdstars_dialog.fgs_mags)) != 0


@pytest.mark.skipif(JENKINS, reason="Can't import PyQt5 on Jenkins server.")
def test_run_background_stars_gui_random(bkgdstars_dialog):
    """Test running the GUI works when using the random method"""
    fgs_mag = 12 # to match the default in the bkgdstars_dialog call

    # Set parameters
    bkgdstars_dialog.groupBox_catalog.setChecked(False)
    bkgdstars_dialog.groupBox_random.setChecked(True)
    bkgdstars_dialog.lineEdit_nStars.setText('5')
    bkgdstars_dialog.lineEdit_magMax.setText('10')
    bkgdstars_dialog.lineEdit_magMin.setText('9')

    # Query and then press OK
    bkgdstars_dialog.pushButton_random.click()
    bkgdstars_dialog.buttonBox.button(QDialogButtonBox.Ok).click()

    # Check number of stars is correct
    assert len(bkgdstars_dialog.x) == 5
    assert len(bkgdstars_dialog.y) == 5
    assert len(list(bkgdstars_dialog.fgs_mags)) == 5

    # Check limits used for colorbar are correct
    assert bkgdstars_dialog.vmin >= fgs_mag + 10

    # Check range of star magnitude is correct
    assert all(fgs_mag + 9  <= mag <= fgs_mag + 10 for mag in list(bkgdstars_dialog.fgs_mags))


@pytest.mark.skipif(JENKINS, reason="Can't import PyQt5 on Jenkins server.")
def test_run_background_stars_gui_user_input(bkgdstars_dialog):
    """Test running the GUI works when using the user-chosen method"""
    # Create a empty row at bottom of table and add data to row
    bkgdstars_dialog.groupBox_catalog.setChecked(False)
    bkgdstars_dialog.groupBox_defined.setChecked(True)

    # Check there's 1 row at present
    assert bkgdstars_dialog.tableWidget.rowCount() == 1
    assert bkgdstars_dialog.tableWidget.columnCount() == 3

    # Set the table values
    bkgdstars_dialog.tableWidget.setItem(0, 0, QTableWidgetItem(str(600)))
    bkgdstars_dialog.tableWidget.setItem(0, 1, QTableWidgetItem(str(700)))
    bkgdstars_dialog.tableWidget.setItem(0, 2, QTableWidgetItem(str(15)))

    # Press OK
    bkgdstars_dialog.buttonBox.button(QDialogButtonBox.Ok).click()

    # Check number of stars is correct
    assert bkgdstars_dialog.x == [600]
    assert bkgdstars_dialog.y == [700]
    assert len(list(bkgdstars_dialog.fgs_mags)) == 1

    # Check the star magnitude is correct
    assert list(bkgdstars_dialog.fgs_mags)[0] == 15
