"""Collection of unit tests to verify the correct function of the Main GUI
of the MAGIC tool.

Authors
-------
    - Shannon Osborne

Use
---
    ::
        pytest test_main_gui.py

Notes
-----
    Add the line `@pytest.mark.skipif(JENKINS, reason="Can't import PyQt5 on Jenkins server.")`
    above tests that call pyqt5
"""
# Standard Library Imports
from mock import patch
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
import jwst_magic

if not JENKINS:
    from jwst_magic.masterGUI import MasterGui

ROOT = "test_main"
__location__ = os.path.realpath(os.path.join(os.getcwd(), os.path.dirname(__file__)))
TEST_DIRECTORY = os.path.join(__location__, 'out', ROOT)

# Set this to True only if computer is on SOGS
COMMISSIONING = True

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


@pytest.mark.skipif(JENKINS, reason="Can't import PyQt5 on Jenkins server.")
@pytest.mark.skipif(not COMMISSIONING, reason="Can't import PyQt5 on Jenkins server.")
def test_update_practice(test_directory):
    """
    Check that re-setting the practice name in the naming method section
    doesn't re-run the fgscountrate call
    """
    # Initialize main window
    app = QApplication(sys.argv)
    master_gui = MasterGui(root='test', in_file=None, out_dir=test_directory,
                           segment_guiding=True, app=app, itm=False)

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

    assert not mock.called, 'method should not have been called'


@pytest.mark.skipif(JENKINS, reason="Can't import PyQt5 on Jenkins server.")
@pytest.mark.skipif(not COMMISSIONING, reason="Can't import PyQt5 on Jenkins server.")
def test_update_car_and_obs(test_directory):
    """
    Check that updating the car and obs in the naming method section
    re-runs the fgscountrate call and updates values

    Note: This test is meant to be run on SOGS and assumings the APT file hasn't been
    updated. If this test fails due to a mis-match of APT information, check the APT
    file by hand
    """
    # Initialize main window
    app = QApplication(sys.argv)
    master_gui = MasterGui(root='test', in_file=None, out_dir=test_directory,
                           segment_guiding=True, app=app, itm=False)

    # Set main GUI parameters
    master_gui.buttonGroup_guider.buttons()[1].setChecked(True)  # set to guider 1

    # Set naming method
    master_gui.buttonGroup_name.buttons()[0].setChecked(True)  # set commissioning naming method
    master_gui.comboBox_practice.setCurrentText('test_practice')
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
@pytest.mark.skipif(not COMMISSIONING, reason="Can't import PyQt5 on Jenkins server.")
def test_update_apt(test_directory):
    """
    Check that if the "update apt" button is pressed in the naming method section
    re-runs the fgscountrate call and updates values

    Note: This test is meant to be run on SOGS and assumings the APT file hasn't been
    updated. If this test fails due to a mis-match of APT information, check the APT
    file by hand
    """
    # Initialize main window
    app = QApplication(sys.argv)
    master_gui = MasterGui(root='test', in_file=None, out_dir=test_directory,
                           segment_guiding=True, app=app, itm=False)

    # Set main GUI parameters
    master_gui.buttonGroup_guider.buttons()[1].setChecked(True)  # set to guider 1

    # Set naming method
    master_gui.buttonGroup_name.buttons()[0].setChecked(True)  # set commissioning naming method
    master_gui.comboBox_practice.setCurrentText('test_practice')
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
@pytest.mark.skipif(not COMMISSIONING, reason="Can't import PyQt5 on Jenkins server.")
def test_apt_guider_disagree(test_directory):
        """
        Check that if the APT program expects a guider that doesn't match
        the guider choosen in MAGIC, a pop up appears
        """
        # Initialize main window
        app = QApplication(sys.argv)
        master_gui = MasterGui(root='test', in_file=None, out_dir=test_directory,
                               segment_guiding=True, app=app, itm=False)

        # For commissioning method


        # For manual method
        master_gui.buttonGroup_name.buttons()[1].setChecked(True)  # set commissioning naming method
        master_gui.lineEdit_root.setText('test')  # set root
        master_gui.textEdit_out.setText(test_directory)  # set out directory


@pytest.mark.skipif(JENKINS, reason="Can't import PyQt5 on Jenkins server.")
def test_apt_gs_populated(test_directory):
    """Test APT info + GS Info are populated into segment guiding SOF GUI correctly"""
    #TODO: this test will be similar to my POF one - so wait for robel to respond
    #TODO: with insight on getting from 1 gui to another
    # Initialize main window
    app = QApplication(sys.argv)
    master_gui = MasterGui(root='test', in_file=None, out_dir=test_directory,
                           segment_guiding=True, app=app, itm=False)

    # Set input image info

    # Set basic info
    master_gui.buttonGroup_name.buttons()[1].setChecked(True)  # set commissioning naming method
    master_gui.lineEdit_root.setText('test')  # set root
    master_gui.textEdit_out.setText(test_directory)  # set out directory
