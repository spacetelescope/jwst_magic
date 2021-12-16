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
import os
import shutil
import sys
from unittest.mock import patch

# Third Party Imports
import numpy as np
JENKINS = '/home/developer/workspace/' in os.getcwd()
if not JENKINS:
    from PyQt5 import QtCore
    from PyQt5.QtWidgets import QDialogButtonBox, QApplication
import pytest

# Local Imports
from jwst_magic.utils import utils

if not JENKINS:
    from jwst_magic.masterGUI import MasterGui

SOGS = utils.on_sogs_network()
if not SOGS:
    from pytestqt import qtbot

ROOT = "test_master"
__location__ = os.path.realpath(os.path.join(os.getcwd(), os.path.dirname(__file__)))
TEST_DIRECTORY = os.path.join(__location__, 'out', ROOT)
DATA_DIRECTORY = os.path.join(__location__, 'data', ROOT)
INPUT_IMAGE = os.path.join(__location__, 'data', 'fgs_data_2_cmimf.fits')
INPUT_IMAGE2 = os.path.join(__location__, 'data', 'nircam_mimf_cal.fits')
COMMAND_FILE = os.path.join(__location__, 'data', 'guiding_selections_test_master_G1.txt')

# Only needed if computer is on SOGS (directory in SOGS_PATH = '***REMOVED***/guiding/')
COM_PRACTICE_DIR = 'magic_pytest_practice'


def delete_contents(directory):
    if os.path.isdir(directory):
        shutil.rmtree(directory)


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
    # global app
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
    master_gui.buttonGroup_guider.buttons()[0].setChecked(True)  # set to guider 2

    # Set naming method
    master_gui.buttonGroup_name.buttons()[0].setChecked(True)  # set commissioning naming method
    master_gui.comboBox_practice.setCurrentIndex(0)
    master_gui.comboBox_car.setCurrentText('OTE-07')
    master_gui.lineEdit_obs.setText('01')

    # Re-set practice
    with patch.object(master_gui, 'update_apt_gs_values') as mock:
        master_gui.comboBox_practice.setCurrentIndex(1)

    assert not mock.called, 'method should not have been called'


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
    master_gui.buttonGroup_guider.buttons()[0].setChecked(True)  # set to guider 2

    # Set naming method
    master_gui.buttonGroup_name.buttons()[0].setChecked(True)  # set commissioning naming method
    master_gui.comboBox_practice.setCurrentText(COM_PRACTICE_DIR)
    master_gui.comboBox_car.setCurrentText('OTE-16')
    master_gui.lineEdit_obs.setText('01')
    master_gui.pushButton_commid.click()  # handle editingfinished not being called in code

    # Check Attributes
    assert master_gui.program_id == 1151
    assert master_gui.observation_num == 1
    assert master_gui.visit_num == 1

    assert master_gui.lineEdit_normalize.text() == 'N4EJ000427'
    assert master_gui.gs_id == 'N4EJ000427'
    np.testing.assert_almost_equal(master_gui.gs_ra, 267.812553, decimal=4)
    np.testing.assert_almost_equal(master_gui.gs_dec, 63.111531, decimal=4)

    # Re-set CAR & OBS
    master_gui.comboBox_car.setCurrentText('OTE-13')
    master_gui.lineEdit_obs.setText('01')
    master_gui.pushButton_commid.click()

    # Re-Check Attributes
    assert master_gui.program_id == 1148
    assert master_gui.observation_num == 1
    assert master_gui.visit_num == 1

    assert master_gui.lineEdit_normalize.text() == 'N4LO000097'
    assert master_gui.gs_id == 'N4LO000097'
    np.testing.assert_almost_equal(master_gui.gs_ra, 254.326768, decimal=4)
    np.testing.assert_almost_equal(master_gui.gs_dec, 68.856895, decimal=4)


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
    master_gui.buttonGroup_guider.buttons()[0].setChecked(True)  # set to guider 2

    # Set naming method
    master_gui.buttonGroup_name.buttons()[0].setChecked(True)  # set commissioning naming method
    master_gui.comboBox_practice.setCurrentText(COM_PRACTICE_DIR)
    master_gui.comboBox_car.setCurrentText('OTE-16')
    master_gui.lineEdit_obs.setText('01')
    master_gui.pushButton_commid.click()

    # Check Attributes
    assert master_gui.program_id == 1151
    assert master_gui.observation_num == 1
    assert master_gui.visit_num == 1

    assert master_gui.lineEdit_normalize.text() == 'N4EJ000427'
    assert master_gui.gs_id == 'N4EJ000427'
    np.testing.assert_almost_equal(master_gui.gs_ra, 267.812553, decimal=4)
    np.testing.assert_almost_equal(master_gui.gs_dec, 63.111531, decimal=4)

    # Re-set APT number and press button
    master_gui.lineEdit_commid.setText('1148')
    master_gui.pushButton_commid.click()
    #QtCore.QTimer.singleShot(0, master_gui.pushButton_commid.clicked)

    # Re-Check Attributes
    assert master_gui.program_id == 1148
    assert master_gui.observation_num == 1
    assert master_gui.visit_num == 1

    assert master_gui.lineEdit_normalize.text() == 'N4LO000097'
    assert master_gui.gs_id == 'N4LO000097'
    np.testing.assert_almost_equal(master_gui.gs_ra, 254.326768, decimal=4)
    np.testing.assert_almost_equal(master_gui.gs_dec, 68.856895, decimal=4)


@pytest.mark.skipif(JENKINS, reason="Can't import PyQt5 on Jenkins server.")
def test_use_apt_button_manual(master_gui, test_directory):
    """
    Test that re-setting apt to blank re-sets program id/obs/visit attributes
    """
    # Set main GUI parameters
    master_gui.buttonGroup_guider.buttons()[0].setChecked(True)  # set to guider 2

    # Set basic info
    master_gui.buttonGroup_name.buttons()[1].setChecked(True)  # set manual naming method
    master_gui.lineEdit_root.setText(ROOT)  # set root
    master_gui.textEdit_out.setText(__location__)  # set out directory
    master_gui.lineEdit_manualid.setText('1151')
    master_gui.lineEdit_manualobs.setText('01')
    master_gui.pushButton_manualid.click()

    # Check Attributes
    assert master_gui.program_id == 1151
    assert master_gui.observation_num == 1
    assert master_gui.visit_num == 1

    assert master_gui.lineEdit_normalize.text() == 'N4EJ000427'
    assert master_gui.gs_id == 'N4EJ000427'
    np.testing.assert_almost_equal(master_gui.gs_ra, 267.812553, decimal=4)
    np.testing.assert_almost_equal(master_gui.gs_dec, 63.111531, decimal=4)

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


apt_parameters = [(pytest.param("commissioning", 0, 'SOF', marks=pytest.mark.skipif(not SOGS, reason="SOGS naming not available"))),
                  (pytest.param("commissioning", 0, 'POF', marks=pytest.mark.skipif(not SOGS, reason="SOGS naming not available"))),
                  (pytest.param("manual", 1, 'SOF', marks=pytest.mark.skipif(SOGS, reason="SOGS naming not available"))),
                  (pytest.param("manual", 1, 'POF', marks=pytest.mark.skipif(SOGS, reason="SOGS naming not available")))]
@pytest.mark.parametrize('type, button_name , filetype', apt_parameters)
@pytest.mark.skipif(JENKINS, reason="Can't import PyQt5 on Jenkins server.")
def test_apt_gs_populated(qtbot, master_gui, test_directory, type, button_name, filetype):
    """
    Test APT info + GS Info are populated into segment guiding SOFand POF GUIs
    correctly for both manual and commissioning naming methods
    """
    # Initialize main window
    qtbot.addWidget(master_gui)

    # Set general input
    qtbot.keyClicks(master_gui.lineEdit_inputImage, INPUT_IMAGE)
    qtbot.mouseClick(master_gui.buttonGroup_name.buttons()[button_name], QtCore.Qt.LeftButton)  # set naming method
    qtbot.mouseClick(master_gui.buttonGroup_guider.buttons()[0], QtCore.Qt.LeftButton)  # set to guider 2

    # Set practice information
    if type == 'commissioning':
        master_gui.comboBox_practice.setCurrentText(COM_PRACTICE_DIR)
        master_gui.comboBox_car.setCurrentText('OTE-13')
        master_gui.lineEdit_obs.setText('01')
        master_gui.pushButton_commid.click()

        assert master_gui.buttonGroup_name.buttons()[0].isChecked()
        assert master_gui.lineEdit_obs.text() == '01'
        assert master_gui.lineEdit_commid.text() == '1148'

        # Delete contents of directory to avoid auto-populating with old data
        delete_contents(master_gui.textEdit_name_preview.toPlainText())

    elif type == 'manual':
        qtbot.keyClicks(master_gui.lineEdit_root, ROOT)  # set root
        qtbot.keyClicks(master_gui.textEdit_out, __location__)  # set out directory
        qtbot.mouseClick(master_gui.buttonGroup_guider.buttons()[1], QtCore.Qt.LeftButton)
        qtbot.mouseClick(master_gui.buttonGroup_guider.buttons()[0], QtCore.Qt.LeftButton)
        qtbot.keyClicks(master_gui.lineEdit_manualid, '1148')
        qtbot.keyClicks(master_gui.lineEdit_manualobs, '01')
        qtbot.mouseClick(master_gui.pushButton_manualid, QtCore.Qt.LeftButton)

        assert master_gui.buttonGroup_name.buttons()[1].isChecked()
        assert master_gui.lineEdit_manualid.text() == '1148'
        assert master_gui.lineEdit_manualobs.text() == '01'
        assert master_gui.lineEdit_root.text() == ROOT
        assert master_gui.textEdit_out.toPlainText() == __location__

    assert master_gui.buttonGroup_guider.checkedButton().text() == '2'
    assert master_gui.program_id == 1148
    assert master_gui.observation_num == 1
    assert master_gui.visit_num == 1

    # Set detector
    master_gui.radioButton_NIRCam.setChecked(True)
    master_gui.comboBox_detector.setCurrentText('B2')
    assert master_gui.comboBox_detector.currentText() == 'B2'

    # Set threshold
    master_gui.lineEdit_threshold.clear()
    qtbot.keyClicks(master_gui.lineEdit_threshold, '0.9')

    # Go to segment guiding
    master_gui.groupBox_imageConverter.setChecked(False)
    master_gui.groupBox_starSelector.setChecked(False)
    master_gui.groupBox_fileWriter.setChecked(False)
    master_gui.groupBox_segmentGuiding.setChecked(True)

    # Set up SOF
    if filetype == 'SOF':
        qtbot.mouseClick(master_gui.radioButton_regfileSegmentGuiding, QtCore.Qt.LeftButton)
        assert master_gui.radioButton_regfileSegmentGuiding.isChecked()

        # Check the pre-populated data (path to root_dir, contents are 0 txt files)
        if type == 'commissioning':
            assert master_gui.lineEdit_regfileSegmentGuiding.text() == master_gui.textEdit_name_preview.toPlainText()
        elif type == 'manual':
            assert master_gui.lineEdit_regfileSegmentGuiding.text() == os.path.join(__location__, 'out', ROOT)
        assert master_gui.comboBox_guidingcommands.count() == 0

        # Change to path that has guiding selections files
        master_gui.lineEdit_regfileSegmentGuiding.clear()
        qtbot.keyClicks(master_gui.lineEdit_regfileSegmentGuiding, os.path.join(__location__, 'data'))
        qtbot.keyClick(master_gui.lineEdit_regfileSegmentGuiding, '\r')  # hit enter
        assert master_gui.lineEdit_regfileSegmentGuiding.text() == os.path.join(__location__, 'data')

        # Check the box that contains the COMMAND_FILE
        i = [i for i in range(master_gui.comboBox_guidingcommands.count()) if COMMAND_FILE.split('/')[-1] in
             master_gui.comboBox_guidingcommands.itemText(i)]
        master_gui.comboBox_guidingcommands.model().item(i[0], 0).setCheckState(QtCore.Qt.Checked)
        assert len(master_gui.comboBox_guidingcommands.checkedItems()) == 1
        assert COMMAND_FILE.split('/')[-1] in master_gui.comboBox_guidingcommands.checkedItems()[0].text()
        assert os.path.isfile(COMMAND_FILE)

    elif filetype == 'POF':
        qtbot.mouseClick(master_gui.radioButton_photometryOverride, QtCore.Qt.LeftButton)

    def handle_dialog():
        try:
            assert master_gui._test_sg_dialog.lineEdit_programNumber.text() == '1148'
            assert master_gui._test_sg_dialog.lineEdit_observationNumber.text() == '1'
            assert master_gui._test_sg_dialog.lineEdit_visitNumber.text() == '1'
            if filetype == 'SOF':
                assert master_gui._test_sg_dialog.lineEdit_RA.text() != ''
                assert master_gui._test_sg_dialog.lineEdit_Dec.text() != ''
                assert master_gui._test_sg_dialog.lineEdit_countrateUncertainty.text() == '0.9'
                assert master_gui._test_sg_dialog.comboBox_detector.currentText() == 'NRCB2'
            qtbot.mouseClick(master_gui._test_sg_dialog.buttonBox.button(QDialogButtonBox.Ok), QtCore.Qt.LeftButton)
        except AssertionError:
            # If something raising an error above, need to close the pop up gui anyway
            qtbot.mouseClick(master_gui._test_sg_dialog.buttonBox.button(QDialogButtonBox.Ok), QtCore.Qt.LeftButton)

    with qtbot.capture_exceptions() as exceptions:
        QtCore.QTimer.singleShot(500, handle_dialog)
        qtbot.mouseClick(master_gui.pushButton_run, QtCore.Qt.LeftButton, delay=1)

    if filetype == 'SOF':
        # check the lack of PA entry causes an error
        expected_err = 'could not convert string to float:'

    elif filetype == 'POF':
        # check the bad countrate value causes an error
        expected_err = 'Countrate factor out of range for count_rate_factor'

    assert expected_err in str(exceptions[0][1]), "Wrong error captured. Caught: '{}', Expected: '{}'".format(
        str(exceptions[0][1]), expected_err)


apt_gs_parameters = [('1451', 1, None, '', '', '', 'gslimits only'),  # Only GS Limits defined
                     ('1141', 1, 1, 'S4OL071716', 159.506895090423, -69.6997532917772, 'gsid and gslimits'),  # Both Limits & ID defined
                     ('1141', 13, None, '', '', '', 'nothing'),  # Nothing defined
                     ]


@pytest.mark.parametrize('apt_number, obs_number, apt_guider, expected_normalize, '
                         'expected_ra, expected_dec, test_type', apt_gs_parameters)
@pytest.mark.skipif(JENKINS, reason="Can't import PyQt5 on Jenkins server.")
def test_apt_gsid_vs_gslimits(qtbot, master_gui, test_directory, apt_number, obs_number, apt_guider,
                              expected_normalize, expected_ra, expected_dec, test_type):
    """Test that the APT files populates the GUI correctly depending on
    if the Guide Star ID or Guide Star Limits are defined
    """
    # Set the guider to be what is in the APT file. If no guider set, set to guider 1
    if apt_guider == 1:
        guider_int = 1
    elif apt_guider == 2:
        guider_int = 0
    else:
        guider_int = 1  # Set guider button to guider 1
    master_gui.buttonGroup_guider.buttons()[guider_int].setChecked(True)  # set to the same guider as APT

    # Set basic info
    master_gui.buttonGroup_name.buttons()[1].setChecked(True)  # set manual naming method
    master_gui.lineEdit_root.setText(ROOT)  # set root
    master_gui.textEdit_out.setText(__location__)  # set out directory
    master_gui.lineEdit_manualid.setText(str(apt_number))
    master_gui.lineEdit_manualobs.setText(str(obs_number))

    if test_type == 'nothing':
        with qtbot.capture_exceptions() as exceptions:
            qtbot.mouseClick(master_gui.pushButton_manualid, QtCore.Qt.LeftButton, delay=1)
        expected_err = "This observation doesn't have a Guide Star Special Requirement"
        assert expected_err in str(exceptions[0][1]), \
            f"Wrong error captured. Caught: '{str(exceptions[0][1])}', Expected: '{expected_err}'"
    else:
        master_gui.pushButton_manualid.click()

    # Check Attributes
    assert master_gui.program_id == int(apt_number)
    assert master_gui.observation_num == int(obs_number)
    assert master_gui.visit_num == 1

    # Check that the guide star ID is populated in the normalize line correctly
    assert master_gui.lineEdit_normalize.text() == expected_normalize
    assert master_gui.gs_id == expected_normalize

    if test_type == 'gslimits only':
        # Check the guider hasn't been changed
        assert master_gui.buttonGroup_name.buttons()[guider_int].isChecked()

        # Check that RA/DEC are empty strings
        assert master_gui.gs_ra == expected_ra
        assert master_gui.gs_dec == expected_dec
    elif test_type == 'gsid and gslimits':
        # Check the guider hasn't been changed
        assert master_gui.buttonGroup_name.buttons()[guider_int].isChecked()

        # Check the RA/DEC come from the GS ID
        np.testing.assert_almost_equal(master_gui.gs_ra, expected_ra, decimal=4)
        np.testing.assert_almost_equal(master_gui.gs_dec, expected_dec, decimal=4)
    elif test_type == 'nothing':
        # Check the guider hasn't been changed
        assert master_gui.buttonGroup_name.buttons()[guider_int].isChecked()

        # Check that RA/DEC are empty strings
        assert master_gui.gs_ra == expected_ra
        assert master_gui.gs_dec == expected_dec


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
    master_gui.lineEdit_obs.setText('01')
    master_gui.pushButton_commid.click()

    # Check setting to the matching guider is fine
    master_gui.buttonGroup_guider.buttons()[1].click()  # set to guider 1
    assert mock_dialog.call_count == 0  # Check dialog box doesn't pop up

    # Check setting the wrong guider fails
    master_gui.buttonGroup_guider.buttons()[0].click()  # set to guider 2
    assert mock_dialog.called  # Check dialog box pops up

    # Check it works setting the guider before choosing CAR/Obs
    master_gui.lineEdit_obs.setText('02')
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
    master_gui.lineEdit_manualid.setText('1151')
    master_gui.lineEdit_manualobs.setText('01')
    master_gui.pushButton_manualid.click()

    # Check setting to the matching guider is fine
    master_gui.buttonGroup_guider.buttons()[0].click()  # set to guider 2
    assert mock_dialog.call_count == 0  # Check dialog box doesn't pop up

    # Check setting the wrong guider fails
    master_gui.buttonGroup_guider.buttons()[1].click()  # set to guider 1
    assert mock_dialog.called  # Check dialog box pops up

    # Check it works setting the guider before choosing CAR/Obs
    master_gui.lineEdit_manualobs.setText('02')
    master_gui.pushButton_manualid.click()
    assert mock_dialog.called  # Check dialog box pops up


@pytest.mark.skipif(JENKINS, reason="Can't import PyQt5 on Jenkins server.")
def test_sg_commands(qtbot, master_gui):
    """
    Test that the segment guiding section of the GUI behaves as expected,
    particularly in terms of populating the checkable combobox
    """

    shifted_file = 'shifted_guiding_selections_test_master_G1_config1.txt'
    unshifted_file = 'unshifted_guiding_selections_test_master_G1_config1.txt'

    # Set General Input
    master_gui.buttonGroup_name.buttons()[1].setChecked(True)  # set manual naming method
    master_gui.lineEdit_root.setText(ROOT)  # set root
    master_gui.textEdit_out.setText(os.path.join(__location__, 'data'))  # set out directory
    master_gui.buttonGroup_guider.buttons()[1].click()  # set to guider 1

    # Go to segment guiding
    master_gui.groupBox_imageConverter.setChecked(False)
    master_gui.groupBox_starSelector.setChecked(False)
    master_gui.groupBox_fileWriter.setChecked(False)
    master_gui.groupBox_segmentGuiding.setChecked(True)

    # Check that the path in the lineedit_regfile_guiding is correct
    assert master_gui.lineEdit_regfileSegmentGuiding.text() == os.path.join(__location__, 'data', 'out', ROOT)

    # Check that the contents of the regfile star selector combobox are correct
    assert master_gui.comboBox_regfileStarSelector.count() == 0

    # Check that the contents of the checkable combobox is shifted only
    cmds = [master_gui.comboBox_guidingcommands.itemText(i) for i in range(master_gui.comboBox_guidingcommands.count())]
    assert len(cmds) == 1
    assert cmds[0] == 'Command 1: ' + shifted_file

    # Switch radio button to unshifted
    master_gui.buttonGroup_segmentGuiding_idAttitude.buttons()[1].click()
    assert master_gui.radioButton_unshifted.isChecked()

    # Check that the path in the lineedit_regfile_guiding hasn't changed
    assert master_gui.lineEdit_regfileSegmentGuiding.text() == os.path.join(__location__, 'data', 'out', ROOT)

    # Check that the contents of the checkable combobox is unshifted only
    cmds = [master_gui.comboBox_guidingcommands.itemText(i) for i in range(master_gui.comboBox_guidingcommands.count())]
    assert len(cmds) == 1
    assert cmds[0] == 'Command 1: ' + unshifted_file

    # Change the path and select a command
    master_gui.lineEdit_regfileSegmentGuiding.clear()
    qtbot.keyClicks(master_gui.lineEdit_regfileSegmentGuiding, os.path.join(__location__, 'data'))
    qtbot.keyClick(master_gui.lineEdit_regfileSegmentGuiding, '\r')  # hit enter
    assert master_gui.lineEdit_regfileSegmentGuiding.text() == os.path.join(__location__, 'data')

    i = [i for i in range(master_gui.comboBox_guidingcommands.count()) if COMMAND_FILE.split('/')[-1] in
         master_gui.comboBox_guidingcommands.itemText(i)]
    master_gui.comboBox_guidingcommands.model().item(i[0], 0).setCheckState(QtCore.Qt.Checked)
    assert len(master_gui.comboBox_guidingcommands.checkedItems()) == 1

    # Flip back to shifted and confirm the path is back to the original
    master_gui.buttonGroup_segmentGuiding_idAttitude.buttons()[0].click()
    assert master_gui.radioButton_shifted.isChecked()
    assert master_gui.lineEdit_regfileSegmentGuiding.text() == os.path.join(__location__, 'data', 'out', ROOT)

    # Flip back to unshifted and confirm the path is also back to the original and nothing is checked
    master_gui.buttonGroup_segmentGuiding_idAttitude.buttons()[1].click()
    assert master_gui.lineEdit_regfileSegmentGuiding.text() == os.path.join(__location__, 'data', 'out', ROOT)
    assert len(master_gui.comboBox_guidingcommands.checkedItems()) == 0


@pytest.mark.skipif(JENKINS, reason="Can't import PyQt5 on Jenkins server.")
def test_shifted_data(master_gui):
    """
    Check that if shifted data is created, the following parts of
    the GUI are updated:
     -  "Use shifted to ID attitude" radio button in SG is selected
     - Contents of combobox are all shifted
     - Shifted file preview is interactive
    """
    shifted_file = 'shifted_guiding_selections_test_master_G1_config1.txt'
    unshifted_file = 'unshifted_guiding_selections_test_master_G1_config1.txt'

    # Set General Input
    master_gui.buttonGroup_name.buttons()[1].setChecked(True)  # set manual naming method
    master_gui.lineEdit_root.setText(ROOT)  # set root
    master_gui.textEdit_out.setText(os.path.join(__location__, 'data'))  # set out directory
    master_gui.buttonGroup_guider.buttons()[1].click()  # set to guider 1

    # Go to segment guiding
    master_gui.groupBox_imageConverter.setChecked(False)
    master_gui.groupBox_starSelector.setChecked(False)
    master_gui.groupBox_fileWriter.setChecked(False)
    master_gui.groupBox_segmentGuiding.setChecked(True)

    # Check that since shifted data is in the out/root directory, that's the default radio button selected
    assert master_gui.radioButton_shifted.isChecked()

    # Check that the path in the lineedit_regfile_guiding is correct
    assert master_gui.lineEdit_regfileSegmentGuiding.text() == os.path.join(__location__, 'data', 'out', ROOT)

    # Check that the contents of the checkable combobox is shifted only
    cmds = [master_gui.comboBox_guidingcommands.itemText(i) for i in range(master_gui.comboBox_guidingcommands.count())]
    assert len(cmds) == 1
    assert cmds[0] == 'Command 1: ' + shifted_file

    # Check file previews
    assert master_gui.comboBox_showcommandsconverted.isEnabled()
    assert master_gui.comboBox_showcommandsshifted.isEnabled()

    # Check 1 option is available in converted image preview
    convert_cmds = [master_gui.comboBox_showcommandsconverted.itemText(i) for i in
                    range(master_gui.comboBox_showcommandsconverted.count())]
    assert len(convert_cmds) == 2
    assert convert_cmds[0] == '- Guiding Command -'
    assert convert_cmds[1] == 'Command 1: ' + unshifted_file

    # Check 1 option is available in shifted image preview
    shift_cmds = [master_gui.comboBox_showcommandsshifted.itemText(i) for i in
                  range(master_gui.comboBox_showcommandsshifted.count())]
    assert len(shift_cmds) == 1
    assert shift_cmds[0] == 'Command 1: ' + shifted_file


@pytest.mark.skipif(JENKINS, reason="Can't import PyQt5 on Jenkins server.")
@pytest.mark.skipif(not SOGS, reason="SOGS naming not available")
def test_commissioning_naming_input_errors(qtbot, master_gui):
    """Test that the correct errors come up when bad inputs are loaded in commissioning naming"""

    # Initialize main window
    qtbot.addWidget(master_gui)

    # Set general input
    qtbot.mouseClick(master_gui.buttonGroup_name.buttons()[0], QtCore.Qt.LeftButton)  # set naming method
    qtbot.mouseClick(master_gui.buttonGroup_guider.buttons()[0], QtCore.Qt.LeftButton)  # set to guider 2
    master_gui.comboBox_practice.setCurrentText(COM_PRACTICE_DIR)
    master_gui.comboBox_car.setCurrentText('OTE-13')
    master_gui.lineEdit_obs.setText('01')

    # Set up file name issue (with a good normalization)
    master_gui.lineEdit_inputImage.setText(INPUT_IMAGE)
    master_gui.comboBox_normalize.setCurrentText('FGS Magnitude')
    master_gui.lineEdit_normalize.setText('14')
    master_gui.pushButton_delbackgroundStars.click()

    with qtbot.capture_exceptions() as exceptions:
        qtbot.mouseClick(master_gui.pushButton_run, QtCore.Qt.LeftButton, delay=1)

    expected_err = "but the user has read in a file that is not a cal.fits file or a padded trk image"
    assert expected_err in str(exceptions[0][1]), \
        f"Wrong error captured. Caught: '{str(exceptions[0][1])}', Expected: '{expected_err}'"

    # Set up count rate issue (with a good image)
    master_gui.lineEdit_inputImage.setText(INPUT_IMAGE2)
    master_gui.comboBox_normalize.setCurrentText('FGS countrate')
    master_gui.lineEdit_normalize.setText('4000000')
    master_gui.pushButton_delbackgroundStars.click()

    with qtbot.capture_exceptions() as exceptions:
        qtbot.mouseClick(master_gui.pushButton_run, QtCore.Qt.LeftButton, delay=1)

    expected_err = "but the user has set the normalization to FGS Countrate"
    assert expected_err in str(exceptions[0][1]), \
        f"Wrong error captured. Caught: '{str(exceptions[0][1])}', Expected: '{expected_err}'"
