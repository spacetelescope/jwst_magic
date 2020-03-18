"""Collection of unit tests to verify the correct function of the Main GUI
of the MAGIC tool.

Authors
-------
    - Shannon Osborne

Use
---
    ::
        pytest test_master_gui.py

Notes
-----
    Add the line `@pytest.mark.skipif(JENKINS, reason="Can't import PyQt5 on Jenkins server.")`
    above tests that call pyqt5
"""
# Standard Library Imports
from unittest.mock import patch
import os
import shutil
import sys

# Third Party Imports
import numpy as np
JENKINS = 'jenkins' in os.getcwd()
if not JENKINS:
    from PyQt5 import QtCore
    from PyQt5.QtWidgets import QDialogButtonBox, QApplication
import pytest

# Local Imports
from ..utils import utils

if not JENKINS:
    from ..masterGUI import MasterGui

SOGS = utils.on_sogs_network()
if not SOGS:
    from pytestqt import qtbot

ROOT = "test_master"
__location__ = os.path.realpath(os.path.join(os.getcwd(), os.path.dirname(__file__)))
TEST_DIRECTORY = os.path.join(__location__, 'out', ROOT)
DATA_DIRECTORY = os.path.join(__location__, 'data', ROOT)
INPUT_IMAGE = os.path.join(__location__, 'data', 'fgs_data_2_cmimf.fits')
COMMAND_FILE = os.path.join(__location__, 'data', 'guiding_selections_for_obs01_G1.txt')

# Only needed if computer is on SOGS (directory in SOGS_PATH = '/data/jwst/wss/guiding/')
COM_PRACTICE_DIR = 'magic_testing_practice'

@pytest.fixture()
def test_directory(test_dir=TEST_DIRECTORY):
    """Create a test directory for permission management.

    Parameters
    ----------
    test_dir : str
        Path to directory used for testing

    Yields
    -------
    test_dir : str
        Path to directory used for testing
    """
    os.makedirs(test_dir)  # creates directory with default mode=511

    yield test_dir
    print("teardown test directory")
    if os.path.isdir(test_dir):
        shutil.rmtree(test_dir)


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


@pytest.mark.skipif(JENKINS, reason="Can't import PyQt5 on Jenkins server.")
@pytest.mark.skipif(not SOGS, reason="SOGS naming not available")
def test_change_practice_commissioning(master_gui):
    """
    Check that re-setting the practice name in the naming method section
    doesn't re-run the fgscountrate call
    """
    # Set main GUI parameters
    master_gui.buttonGroup_guider.buttons()[1].setChecked(True)  # set to guider 1

    # Set naming method
    master_gui.buttonGroup_name.buttons()[0].setChecked(True)  # set commissioning naming method
    master_gui.comboBox_practice.setCurrentIndex(0)
    master_gui.comboBox_car.setCurrentText('OTE-07')
    master_gui.comboBox_obs.setCurrentText('01')

    # Re-set practice
    with patch.object(master_gui, 'update_apt_gs_values') as mock:
        master_gui.comboBox_practice.setCurrentIndex(1)

    assert not mock.called, 'method should not have been called'''


@pytest.mark.skipif(JENKINS, reason="Can't import PyQt5 on Jenkins server.")
@pytest.mark.skipif(not SOGS, reason="SOGS naming not available")
def test_change_car_and_obs_commissioning(master_gui):
    """
    Check that updating the car and obs in the naming method section
    re-runs the fgscountrate call and updates values

    Note: This test is meant to be run on SOGS and assumings the APT file hasn't been
    updated. If this test fails due to a mis-match of APT information, check the APT
    file by hand
    """
    # Set main GUI parameters
    master_gui.buttonGroup_guider.buttons()[1].setChecked(True)  # set to guider 1

    # Set naming method
    master_gui.buttonGroup_name.buttons()[0].setChecked(True)  # set commissioning naming method
    master_gui.comboBox_practice.setCurrentText(COM_PRACTICE_DIR)
    master_gui.comboBox_car.setCurrentText('OTE-07')
    master_gui.comboBox_obs.setCurrentText('01')

    # Check Attributes
    assert master_gui.program_id == 1141
    assert master_gui.observation_num == 1
    assert master_gui.visit_num == 1

    assert master_gui.lineEdit_normalize.text() == 'N13I018276'
    assert master_gui.gs_id == 'N13I018276'
    np.testing.assert_almost_equal(master_gui.gs_ra, 273.14572081885797)
    np.testing.assert_almost_equal(master_gui.gs_dec, 65.5301991668692)

    # Re-set CAR & OBS
    master_gui.comboBox_car.setCurrentText('OTE-13')
    master_gui.comboBox_obs.setCurrentText('01')

    # Re-Check Attributes
    assert master_gui.program_id == 1148
    assert master_gui.observation_num == 1
    assert master_gui.visit_num == 1

    assert master_gui.lineEdit_normalize.text() == 'N13I000018'
    assert master_gui.gs_id == 'N13I000018'
    np.testing.assert_almost_equal(master_gui.gs_ra, 273.206729760604)
    np.testing.assert_almost_equal(master_gui.gs_dec, 65.5335149359777)


@pytest.mark.skipif(JENKINS, reason="Can't import PyQt5 on Jenkins server.")
@pytest.mark.skipif(not SOGS, reason="SOGS naming not available")
def test_update_apt_button_commissioning(master_gui):
    """
    Check that if the "update apt" button is pressed in the naming method section
    re-runs the fgscountrate call and updates values

    Note: This test is meant to be run on SOGS and assumings the APT file hasn't been
    updated. If this test fails due to a mis-match of APT information, check the APT
    file by hand
    """
    # Set main GUI parameters
    master_gui.buttonGroup_guider.buttons()[1].setChecked(True)  # set to guider 1

    # Set naming method
    master_gui.buttonGroup_name.buttons()[0].setChecked(True)  # set commissioning naming method
    master_gui.comboBox_practice.setCurrentText(COM_PRACTICE_DIR)
    master_gui.comboBox_car.setCurrentText('OTE-07')
    master_gui.comboBox_obs.setCurrentText('01')

    # Check Attributes
    assert master_gui.program_id == 1141
    assert master_gui.observation_num == 1
    assert master_gui.visit_num == 1

    assert master_gui.lineEdit_normalize.text() == 'N13I018276'
    assert master_gui.gs_id == 'N13I018276'
    np.testing.assert_almost_equal(master_gui.gs_ra, 273.14572081885797)
    np.testing.assert_almost_equal(master_gui.gs_dec, 65.5301991668692)

    # Re-set APT number and press button
    master_gui.lineEdit_commid.setText('1148')
    master_gui.pushButton_commid.click()
    #QtCore.QTimer.singleShot(0, master_gui.pushButton_commid.clicked)

    # Re-Check Attributes
    assert master_gui.program_id == 1148
    assert master_gui.observation_num == 1
    assert master_gui.visit_num == 1

    assert master_gui.lineEdit_normalize.text() == 'N13I000018'
    assert master_gui.gs_id == 'N13I000018'
    np.testing.assert_almost_equal(master_gui.gs_ra, 273.206729760604)
    np.testing.assert_almost_equal(master_gui.gs_dec, 65.5335149359777)


@pytest.mark.skipif(JENKINS, reason="Can't import PyQt5 on Jenkins server.")
def test_use_apt_button_manual(master_gui, test_directory):
    """
    Test that re-setting apt to blank re-sets program id/obs/visit attributes
    """
    # Set main GUI parameters
    master_gui.buttonGroup_guider.buttons()[1].setChecked(True)  # set to guider 1

    # Set basic info
    master_gui.buttonGroup_name.buttons()[1].setChecked(True)  # set manual naming method
    master_gui.lineEdit_root.setText(ROOT)  # set root
    master_gui.textEdit_out.setText(__location__)  # set out directory
    master_gui.lineEdit_manualid.setText('1141')
    master_gui.lineEdit_manualobs.setText('01')
    master_gui.pushButton_manualid.click()

    # Check Attributes
    assert master_gui.program_id == 1141
    assert master_gui.observation_num == 1
    assert master_gui.visit_num == 1

    assert master_gui.lineEdit_normalize.text() == 'N13I018276'
    assert master_gui.gs_id == 'N13I018276'
    np.testing.assert_almost_equal(master_gui.gs_ra, 273.14572081885797)
    np.testing.assert_almost_equal(master_gui.gs_dec, 65.5301991668692)

    # Reset to blank
    master_gui.lineEdit_manualid.setText('')
    master_gui.lineEdit_manualobs.setText('')
    master_gui.pushButton_manualid.click()

    # Re-check attributes
    assert master_gui.program_id == ''
    assert master_gui.observation_num == ''
    assert master_gui.visit_num == ''

    assert master_gui.lineEdit_normalize.text() == ''
    assert master_gui.gs_id == ''
    assert master_gui.gs_ra == ''
    assert master_gui.gs_dec == ''


apt_parameters = [pytest.param("commissioning", 0, marks=pytest.mark.skipif(not SOGS, reason="SOGS naming not available")),
                  (pytest.param("manual", 1, marks=pytest.mark.skipif(SOGS, reason="SOGS naming not available")))]
@pytest.mark.parametrize('type, button_name', apt_parameters)
@pytest.mark.skipif(JENKINS, reason="Can't import PyQt5 on Jenkins server.")
def test_apt_gs_populated(qtbot, master_gui, test_directory, type, button_name):
    """
    Test APT info + GS Info are populated into segment guiding SOF GUI correctly
    for both manual and commissioning naming methods
    """
    # Initialize main window
    qtbot.addWidget(master_gui)

    # Set general input
    qtbot.keyClicks(master_gui.lineEdit_inputImage, INPUT_IMAGE)
    qtbot.mouseClick(master_gui.buttonGroup_name.buttons()[button_name], QtCore.Qt.LeftButton)  # set naming method
    qtbot.mouseClick(master_gui.buttonGroup_guider.buttons()[1], QtCore.Qt.LeftButton)  # set to guider 1

    # Set practice information
    if type == 'commissioning':
        master_gui.comboBox_practice.setCurrentText(COM_PRACTICE_DIR)
        master_gui.comboBox_car.setCurrentText('OTE-13')
        master_gui.comboBox_obs.setCurrentText('01')

        assert master_gui.buttonGroup_name.buttons()[0].isChecked()
        assert master_gui.comboBox_obs.currentText() == '01'
        assert master_gui.lineEdit_commid.text() == '1148'

    elif type == 'manual':
        qtbot.keyClicks(master_gui.lineEdit_root, ROOT)  # set root
        qtbot.keyClicks(master_gui.textEdit_out, __location__)  # set out directory
        qtbot.keyClicks(master_gui.lineEdit_manualid, '1148')
        qtbot.keyClicks(master_gui.lineEdit_manualobs, '01')
        qtbot.mouseClick(master_gui.pushButton_manualid, QtCore.Qt.LeftButton)

        assert master_gui.buttonGroup_name.buttons()[1].isChecked()
        assert master_gui.lineEdit_manualid.text() == '1148'
        assert master_gui.lineEdit_manualobs.text() == '01'
        assert master_gui.lineEdit_root.text() == ROOT
        assert master_gui.textEdit_out.toPlainText() == __location__

    assert master_gui.buttonGroup_guider.checkedButton().text() == '1'
    assert master_gui.program_id == 1148
    assert master_gui.observation_num == 1
    assert master_gui.visit_num == 1

    # Go to segment guiding
    #qtbot.mouseClick(master_gui.checkBox_useConvertedImage, QtCore.Qt.LeftButton)
    master_gui.groupBox_imageConverter.setChecked(False)
    master_gui.groupBox_starSelector.setChecked(False)
    master_gui.groupBox_fileWriter.setChecked(False)
    master_gui.groupBox_segmentGuiding.setChecked(True)

    # Set up SOF
    qtbot.mouseClick(master_gui.radioButton_regfileSegmentGuiding, QtCore.Qt.LeftButton)
    master_gui.lineEdit_regfileSegmentGuiding.clear()
    qtbot.keyClicks(master_gui.lineEdit_regfileSegmentGuiding, COMMAND_FILE)

    assert master_gui.radioButton_regfileSegmentGuiding.isChecked()
    assert master_gui.lineEdit_regfileSegmentGuiding.text() == COMMAND_FILE
    assert os.path.isfile(COMMAND_FILE)

    def handle_dialog():
        try:
            assert master_gui._test_sg_dialog.lineEdit_programNumber.text() == '1148'
            assert master_gui._test_sg_dialog.lineEdit_observationNumber.text() == '1'
            assert master_gui._test_sg_dialog.lineEdit_visitNumber.text() == '1'
            assert master_gui._test_sg_dialog.lineEdit_RA.text() != ''
            assert master_gui._test_sg_dialog.lineEdit_Dec.text() != ''
            qtbot.mouseClick(master_gui._test_sg_dialog.buttonBox.button(QDialogButtonBox.Ok), QtCore.Qt.LeftButton)
        except AssertionError:
            # If something raising an error above, need to close the pop up gui anyway
            qtbot.mouseClick(master_gui._test_sg_dialog.buttonBox.button(QDialogButtonBox.Ok), QtCore.Qt.LeftButton)

    with qtbot.capture_exceptions() as exceptions:
        QtCore.QTimer.singleShot(500, handle_dialog)
        qtbot.mouseClick(master_gui.pushButton_run, QtCore.Qt.LeftButton, delay=1)

    # check the lack of PA entry causes an error
    expected_err =  'could not convert string to float:'
    assert expected_err in str(exceptions[0][1]), "Wrong error captured. Caught: '{}', Expected: '{}'".format(
        str(exceptions[0][1]), expected_err)


@pytest.mark.skipif(JENKINS, reason="Can't import PyQt5 on Jenkins server.")
@pytest.mark.skipif(not SOGS, reason="SOGS naming not available")
@patch('jwst_magic.masterGUI.MasterGui.mismatched_apt_guider_dialog', autospec=True)
def test_apt_guider_disagree_commissioning(mock_dialog, master_gui):
    """
    Check that if the APT program expects a guider that doesn't match
    the guider chosen in MAGIC, a pop up appears
    """
    # For commissioning method
    # Guider set after choosing CAR/Obs
    master_gui.buttonGroup_name.buttons()[0].setChecked(True)  # set commissioning naming method
    master_gui.comboBox_practice.setCurrentText(COM_PRACTICE_DIR)
    master_gui.comboBox_car.setCurrentText('OTE-07')
    master_gui.comboBox_obs.setCurrentText('01')

    # Check setting to the matching guider is fine
    master_gui.buttonGroup_guider.buttons()[1].click()  # set to guider 1
    assert mock_dialog.call_count == 0  # Check dialog box doesn't pop up

    # Check setting the wrong guider fails
    master_gui.buttonGroup_guider.buttons()[0].click()  # set to guider 2
    assert mock_dialog.called  # Check dialog box pops up

    # Check it works setting the guider before choosing CAR/Obs
    master_gui.comboBox_obs.setCurrentText('02')
    assert mock_dialog.called  # Check dialog box pops up


@pytest.mark.skipif(JENKINS, reason="Can't import PyQt5 on Jenkins server.")
@patch("jwst_magic.masterGUI.MasterGui.mismatched_apt_guider_dialog", autospec=True)
def test_apt_guider_disagree_manual(mock_dialog, master_gui, test_directory):
    """
    Check that if the APT program expects a guider that doesn't match
    the guider chosen in MAGIC, a pop up appears
    """
    # For manual method
    # Guider set after choosing CAR/Obs
    master_gui.buttonGroup_name.buttons()[1].setChecked(True)  # set manual naming method
    master_gui.lineEdit_root.setText(ROOT)  # set root
    master_gui.textEdit_out.setText(__location__)  # set out directory
    master_gui.lineEdit_manualid.setText('1141')
    master_gui.lineEdit_manualobs.setText('01')
    master_gui.pushButton_manualid.click()

    # Check setting to the matching guider is fine
    master_gui.buttonGroup_guider.buttons()[1].click()  # set to guider 1
    assert mock_dialog.call_count == 0  # Check dialog box doesn't pop up

    # Check setting the wrong guider fails
    master_gui.buttonGroup_guider.buttons()[0].click()  # set to guider 2
    assert mock_dialog.called  # Check dialog box pops up

    # Check it works setting the guider before choosing CAR/Obs
    master_gui.lineEdit_manualobs.setText('02')
    master_gui.pushButton_manualid.click()
    assert mock_dialog.called  # Check dialog box pops up
