"""Collection of unit tests to verify the correct function of the Segment
Guiding tool.

TODO
(In the future, this will) Runs the following tests:
    1. Test all of the different parameters in run_tool individually.
    2. Test with the params dict with ra, dec, and/or pa set to None
    3. Test with the params dict with ra, dec, and/or pa set to floats
    4. Test with the params dict with ra, dec, and/or pa set to ints
    5. Test with the params dict with ra, dec, and/or pa set out of range
    6. Test with parameter_dialog turned on with ra, dec, and/or pa set to None
    7. Test with parameter_dialog turned on with ra, dec, and/or pa set to floats
    8. Test with parameter_dialog turned on with ra, dec, and/or pa set to ints
    9. Test with parameter_dialog turned on with ra, dec, and/or pa set out of range
    10. Test with program_id, observation_num, and/or visit_num set to ints
    11. Test with program_id, observation__num, and/or visit_num set to strings
    12. Test with program_id, observation__num, and/or visit_num set to floats
    13. Test with SEGMENT_INFILE and/or guide set to None

Authors
-------
    - Keira Brooks
    - Lauren Chambers

Use
---
    ::
        pytest test_segment_guiding.py

Notes
-----

The following are ways to run segment_guiding:

    from jwst_magic.segment_guiding import segment_guiding
    segment_guiding.run_tool(segment_infile=SEGMENT_INFILE, guider=guider)

    Or with the segment dialog box:
    ::
    from jwst_magic.segment_guiding import segment_guiding
    segment_guiding.run_tool(program_id=program_id, observation=observation,
                             visit=visit, parameter_dialog=True)

    Or from dictionary of parameters:
    ::
    from jwst_magic.segment_guiding import segment_guiding
    segment_guiding.run_tool(program_id=program_id, observation=observation,
                             visit=visit, guide_star_params_dict=guide_star_params_dict,
                             parameter_dialog=False)

As of 10/23/2019 there are 39 tests for this module. 4 tests cannot be run with
Jenkins because they test the GUI.
"""
# Standard Library Imports
from datetime import datetime
import os
import shutil

# Third Party Imports
import numpy as np
JENKINS = 'jenkins' in os.getcwd()
if not JENKINS:
    from PyQt5 import QtCore
    from PyQt5.QtWidgets import QDialogButtonBox
import pytest

# Local Imports
from utils import parametrized_data
from jwst_magic.segment_guiding.segment_guiding import (generate_segment_override_file,
                                                        SegmentGuidingCalculator,
                                                        generate_photometry_override_file)
if not JENKINS:
    from jwst_magic.segment_guiding.SegmentGuidingGUI import SegmentGuidingDialog


__location__ = os.path.realpath(os.path.join(os.getcwd(), os.path.dirname(__file__)))
ROOT = "test_sgt"
SEGMENT_INFILE = os.path.join(__location__, 'data', '{}_ALLpsfs.txt'.format(ROOT))
SELECTED_SEGS = os.path.join(__location__, 'data', '{}_regfile.txt'.format(ROOT))
PROGRAM_ID = 1141
OBSERVATION_NUM = 7
VISIT_NUM = 1
TEST_DIRECTORY = os.path.join(__location__, 'out', ROOT)

PARAMETRIZED_DATA = parametrized_data()['test_segment_guiding']


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


test_data = PARAMETRIZED_DATA['test_generate_segment_override_file']
sof_parameters = [(0, None, test_data[0]),
                  (4, None, test_data[1]),
                  (0, np.array([[1, 12, 6]]), test_data[2]),
                  (0, np.array([[1, 2, 3], [4, 0, 17, 12, 2]]), test_data[3])]
@pytest.mark.parametrize('seg_num, selected_segs, correct_command', sof_parameters)
def test_generate_segment_override_file(test_directory, seg_num, selected_segs, correct_command):

    # Define the input file locations and parameters
    if selected_segs is None:
        selected_segs = SELECTED_SEGS
    guider = 1


    guide_star_params_dict = {'v2_boff': 0.1,
                              'v3_boff': 0.2,
                              'fgs_num': guider,
                              'ra': 90.9708,
                              'dec': -67.3578,
                              'pa': 157.1234,
                              'seg_num': seg_num}

    generate_segment_override_file(
        SEGMENT_INFILE, guider, PROGRAM_ID, OBSERVATION_NUM, VISIT_NUM, root=ROOT,
        out_dir=__location__, selected_segs=selected_segs, click_to_select_gui=False,
        guide_star_params_dict=guide_star_params_dict, parameter_dialog=False
    )

    # Check to make sure the override file was created, and in the right place
    segment_override_file = os.path.join(test_directory, '{}_gs_override_1141_7_1.txt'.format(datetime.now().strftime('%Y%m%d')))
    assert os.path.isfile(segment_override_file)

    # # Check to make sure the command was written correctly
    with open(segment_override_file) as f:
        segment_override_command = f.read()
    assert segment_override_command == correct_command


sof_valueerror_parameters = [(20, 1, 'out of range'),
                             ('zero', 2, 'invalid literal for int()'),
                             (0, 3, 'Invalid guider number')]
@pytest.mark.parametrize('seg_num, guider, error_text', sof_valueerror_parameters)
def test_segment_guiding_calculator_valueerrors(test_directory, seg_num, guider, error_text):
    # Define the input file locations and parameters
    guide_star_params_dict = {'v2_boff': 0.1,
                              'v3_boff': 0.2,
                              'fgs_num': guider,
                              'ra': 90.9708,
                              'dec': -67.3578,
                              'pa': 157.1234,
                              'seg_num': seg_num}

    with pytest.raises(ValueError) as excinfo:
        generate_segment_override_file(
            SEGMENT_INFILE, guider, PROGRAM_ID, OBSERVATION_NUM, VISIT_NUM,
            root=ROOT, out_dir=__location__, selected_segs=SELECTED_SEGS,
            click_to_select_gui=False, guide_star_params_dict=guide_star_params_dict,
            parameter_dialog=False
        )
    assert error_text in str(excinfo.value)

def test_generate_override_file_valueerrors(test_directory):
    # Define the input file locations and parameters
    guider = 1
    guide_star_params_dict = {'v2_boff': 0.1,
                              'v3_boff': 0.2,
                              'fgs_num': guider,
                              'ra': 90.9708,
                              'dec': -67.3578,
                              'pa': 157.1234,
                              'seg_num': 0}

    # Test error if parameter_dialog=False and guide_star_params_dict=None
    with pytest.raises(ValueError) as excinfo:
        generate_segment_override_file(
            SEGMENT_INFILE, guider, PROGRAM_ID, OBSERVATION_NUM, VISIT_NUM,
            root=ROOT, out_dir=__location__, guide_star_params_dict=None,
            parameter_dialog=False
        )
    assert '`parameter_dialog=False`' in str(excinfo.value)

    # Test error if click_to_select_gui=True and data=None
    with pytest.raises(ValueError) as excinfo:
        generate_segment_override_file(
            SEGMENT_INFILE, guider, PROGRAM_ID, OBSERVATION_NUM, VISIT_NUM,
            root=ROOT, out_dir=__location__,
            click_to_select_gui=True, guide_star_params_dict=guide_star_params_dict,
            parameter_dialog=False
        )
    assert '`click_to_select_GUI=True`' in str(excinfo.value)

    # Test error if click_to_select_gui=False and selected_segs=None
    with pytest.raises(ValueError) as excinfo:
        generate_segment_override_file(
            SEGMENT_INFILE, guider, PROGRAM_ID, OBSERVATION_NUM, VISIT_NUM,
            root=ROOT, out_dir=__location__,
            click_to_select_gui=False, guide_star_params_dict=guide_star_params_dict,
            parameter_dialog=False
        )
    assert '`click_to_select_GUI=False`' in str(excinfo.value)

    # Test error if parameter_dialog=False and countrate_factor=None
    with pytest.raises(ValueError) as excinfo:
        generate_photometry_override_file(
            ROOT, PROGRAM_ID, OBSERVATION_NUM, VISIT_NUM,
            out_dir=__location__, parameter_dialog=False
        )
    assert '`parameter_dialog=False`' in str(excinfo.value)


def test_segment_override_command_out_of_fov(test_directory):
    # Define the input file locations and parameters
    guide_star_params_dict = {'v2_boff': 0.1,
                              'v3_boff': 0.2,
                              'fgs_num': 1,
                              'ra': 90.9708,
                              'dec': -67.3578,
                              'pa': 157.1234,
                              'seg_num': 0}

    sg = SegmentGuidingCalculator(
        "SOF", PROGRAM_ID, OBSERVATION_NUM, VISIT_NUM, ROOT, __location__,
        segment_infile=SEGMENT_INFILE,
        guide_star_params_dict=guide_star_params_dict,
        selected_segs=SELECTED_SEGS
    )

    # Determine the V2/V3 of the pointing center
    sg.get_center_pointing()

    # Calculate the RA/Dec of each segment
    sg.calculate_effective_ra_dec()

    # Make a phony attitude matrix
    attitude = np.random.rand(3, 3)

    # Check with phony attitude matrix
    with pytest.raises(ValueError) as excinfo:
        sg.check_segments_inside_fov(attitude)
    assert 'Incorrect segment guiding calculations' in str(excinfo.value)

    # # Check with phony X/Y locations between -5000 and 5000
    # sg.x_det_segs = np.random.rand(18) * 10000 - 5000
    # sg.y_det_segs = np.random.rand(18) * 10000 - 5000
    # with pytest.raises(ValueError):
    #     sg.check_segments_inside_fov(attitude)
    # # Check that the log is raising a warning
    # assert 'off detector' in caplog.text


def test_generate_photometry_override_file(test_directory):
    generate_photometry_override_file(
        ROOT, PROGRAM_ID, OBSERVATION_NUM, VISIT_NUM, countrate_factor=0.7,
        out_dir=__location__, parameter_dialog=False
    )

    # Check to make sure the override file was created, and in the right place
    segment_override_file = os.path.join(__location__, 'out', ROOT,
                                         '{}_gs_override_1141_7_1.txt'.format(datetime.now().strftime('%Y%m%d')))
    assert os.path.isfile(segment_override_file)

    # Check to make sure the command was written correctly
    with open(segment_override_file) as f:
        segment_override_command = f.read()
    assert segment_override_command == 'sts -gs_select 01141007001 -count_rate_factor=0.700'

    # Try again with an incorrect countrate factor and make sure there's an error
    with pytest.raises(ValueError) as excinfo:
        generate_photometry_override_file(
            ROOT, PROGRAM_ID, OBSERVATION_NUM, VISIT_NUM, countrate_factor=2.0,
            out_dir=__location__, parameter_dialog=False)
    assert 'for count_rate_factor. Expecting between 0.0 and 1.0.' in str(excinfo.value)


@pytest.mark.skipif(JENKINS, reason="Can't import PyQt5 on Jenkins server.")
def test_cancel_file_dialog():
    """Raise a segment override file dialog window and cancel it.
    """
    # Initialize dialog window
    segment_guiding_dialog = SegmentGuidingDialog("SOF", 1, PROGRAM_ID, OBSERVATION_NUM, VISIT_NUM)

    # Schedule press of "Cancel" button
    cancel_button = segment_guiding_dialog.buttonBox.button(QDialogButtonBox.Cancel)
    QtCore.QTimer.singleShot(0, cancel_button.clicked)

    # Run GUI
    accepted = segment_guiding_dialog.exec()

    assert not accepted


@pytest.mark.skipif(JENKINS, reason="Can't import PyQt5 on Jenkins server.")
def test_SOF_parameters_dialog():
    # Initialize dialog window
    segment_guiding_dialog = SegmentGuidingDialog("SOF", 1, PROGRAM_ID, OBSERVATION_NUM, VISIT_NUM)

    # Set parameters
    segment_guiding_dialog.lineEdit_RA.setText('90.9708')
    segment_guiding_dialog.comboBox_RAUnit.setCurrentText('Degrees')
    segment_guiding_dialog.lineEdit_Dec.setText('-67.3578')
    segment_guiding_dialog.lineEdit_PA.setText('157.1234')

    # Schedule press of "Ok" button
    ok_button = segment_guiding_dialog.buttonBox.button(QDialogButtonBox.Ok)
    QtCore.QTimer.singleShot(0, ok_button.clicked)

    # Run GUI
    segment_guiding_dialog.exec()

    # Get parameters
    params = segment_guiding_dialog.get_dialog_parameters()

    assert params == (
        {'v3_boff': 0.0, 'seg_num': 0, 'v2_boff': 0.0, 'ra': 90.9708, 'fgs_num': 1, 'dec': -67.3578, 'pa': 157.1234},
        '1141', '7', '1', 0.6, None
    )

pof_dialog_parameters = [('1142', '8', '2', None, (None, '1142', '8', '2', None, 0.0)),
                         ('1142', '8', '2', 0.0123, (None, '1142', '8', '2', None, 0.0123))]
@pytest.mark.parametrize('program_id, obs_num, visit_num, countrate_factor, out_params', pof_dialog_parameters)
@pytest.mark.skipif(JENKINS, reason="Can't import PyQt5 on Jenkins server.")
def test_POF_parameters_dialog(program_id, obs_num, visit_num, countrate_factor, out_params):
    # Initialize dialog window
    segment_guiding_dialog = SegmentGuidingDialog("POF", 1, PROGRAM_ID, OBSERVATION_NUM, VISIT_NUM)

    # Set parameters
    segment_guiding_dialog.lineEdit_programNumber.setText(program_id)
    segment_guiding_dialog.lineEdit_observationNumber.setText(obs_num)
    segment_guiding_dialog.lineEdit_visitNumber.setText(visit_num)
    if countrate_factor is not None:
        segment_guiding_dialog.doubleSpinBox_countrateFactor.setValue(countrate_factor)

    # Schedule press of "Ok" button
    ok_button = segment_guiding_dialog.buttonBox.button(QDialogButtonBox.Ok)
    QtCore.QTimer.singleShot(0, ok_button.clicked)

    # Run GUI
    segment_guiding_dialog.exec()

    # Get parameters
    params = segment_guiding_dialog.get_dialog_parameters()

    # params = (guide_star_params_dict, program_id, observation_num, visit_num,
    #           threshold_factor, countrate_factor)
    assert params == out_params

def test_write_override_report(test_directory):
    # Define the input file locations and parameters
    guider = 1
    guide_star_params_dict = {'v2_boff': 0.1,
                              'v3_boff': 0.2,
                              'fgs_num': guider,
                              'ra': 90.9708,
                              'dec': -67.3578,
                              'pa': 157.1234,
                              'seg_num': 0}
    selected_segs = np.array([[1, 2, 3], [4, 0, 17, 12, 2]])

    generate_segment_override_file(
        SEGMENT_INFILE, guider, PROGRAM_ID, OBSERVATION_NUM, VISIT_NUM, root=ROOT,
        out_dir=__location__, selected_segs=selected_segs, click_to_select_gui=False,
        guide_star_params_dict=guide_star_params_dict, parameter_dialog=False
    )

    report_file = '{}_gs_override_{}_{}_{}_REPORT.txt'.format(datetime.now().strftime('%Y%m%d'),
                                                              PROGRAM_ID,
                                                              OBSERVATION_NUM,
                                                              VISIT_NUM)
    report_file = os.path.join(test_directory, report_file)
    assert os.path.isfile(report_file)

    correct_file = '''Guide Star Override Report
Generated on 2018/12/07 at 11:59:16 by lchambers

Program ID    : 1141
Observation # : 7
Visit #       : 1
FGS Detector  : 1
Guide Star RA : 90.970800
Guide Star Dec: -67.357800
V3 PA @ GS    : 157.123400

  Star Name  |   MAGIC ID   |      RA      |      Dec     |    Ideal X   |    Ideal Y   |  OSS Ideal X |  OSS Ideal Y |  Detector X  |  Detector Y
-----------------------------------------------------------------------------------------------------------------------------------------------------
star1        | 1            | 90.987850    | -67.355022   | -12.648238   | -22.174651   | 12.648238    | -22.174651   | 840.000187   | 1344.998823
star2        | 4            | 90.986120    | -67.362103   | -0.034250    | 0.106095     | 0.034250     | 0.106095     | 1024.000008  | 1023.000013
ref_only1    | 2            | 90.986990    | -67.358569   | -6.338517    | -11.010128   | 6.338517     | -11.010128   | 932.000192   | 1183.999172
ref_only2    | 3            | 90.980194    | -67.352832   | -6.188804    | -33.708446   | 6.188804     | -33.708446   | 933.999283   | 1504.998775
ref_only3    | 12           | 90.963123    | -67.355541   | 19.393433    | -34.482825   | -19.393433   | -34.482825   | 1304.997893  | 1503.997675
ref_only4    | 17           | 90.954592    | -67.356915   | 32.206292    | -34.798421   | -32.206292   | -34.798421   | 1489.998298  | 1501.997517
ref_only5    | 0            | 90.953775    | -67.360493   | 38.500751    | -23.505004   | -38.500751   | -23.505004   | 1581.997741  | 1339.996887
'''.split('\n')

    with open(report_file) as f:
        lines = f.read().split('\n')

    for l, c in zip(lines[2:], correct_file[2:]):
        assert l.rstrip() == c.rstrip()


test_data = PARAMETRIZED_DATA['test_segment_override_file_wo_obs_visit']
sof_wo_obs_visit_parameters = [(OBSERVATION_NUM, '', 'gs_override_1141_7.txt', test_data[0]),
                               ('', '', 'gs_override_1141.txt', test_data[1]),
                               ('', 1, 'gs_override_1141.txt', test_data[1]), #This test will have the same result as the previous
                               ('1-2', 3, 'gs_override_1141_1-2_1.txt', test_data[2])]
@pytest.mark.parametrize('obsnum, visitnum, out_file, correct_command', sof_wo_obs_visit_parameters)
def test_segment_override_file_wo_obs_visit(test_directory, obsnum, visitnum, out_file, correct_command):

    guider = 1
    guide_star_params_dict = {'v2_boff': 0.1,
                              'v3_boff': 0.2,
                              'fgs_num': guider,
                              'ra': 90.9708,
                              'dec': -67.3578,
                              'pa': 157.1234,
                              'seg_num': 0}

    generate_segment_override_file(
        SEGMENT_INFILE, guider, PROGRAM_ID, obsnum, visitnum, root=ROOT,
        out_dir=__location__, selected_segs=SELECTED_SEGS, click_to_select_gui=False,
        guide_star_params_dict=guide_star_params_dict, parameter_dialog=False
    )

    # Check to make sure the override file was created, and in the right place
    segment_override_file = os.path.join(test_directory, '{}_{}'.format(datetime.now().strftime('%Y%m%d'), out_file))
    assert os.path.isfile(segment_override_file)

    # # Check to make sure the command was written correctly
    with open(segment_override_file) as f:
        segment_override_command = f.read()
    assert segment_override_command == correct_command

test_data = PARAMETRIZED_DATA['test_split_obs_num']
split_obs_num_parameters = [('2', *test_data[0]),
                            ('2, 4, 6', *test_data[1]),
                            ('13-15', *test_data[2]),
                            ('5, 7, 10-12', *test_data[3]),
                            ('5, 3, 1', *test_data[4]),
                            ('10-12, 5, 7', *test_data[5])]
@pytest.mark.parametrize('obs_num, correct_list, correct_string', split_obs_num_parameters)
def test_split_obs_num(obs_num, correct_list, correct_string):
    # Define the input file locations and parameters
    guide_star_params_dict = {'v2_boff': 0.1,
                              'v3_boff': 0.2,
                              'fgs_num': 1,
                              'ra': 90.9708,
                              'dec': -67.3578,
                              'pa': 157.1234,
                              'seg_num': 0}

    sg = SegmentGuidingCalculator(
        "SOF", PROGRAM_ID, obs_num, VISIT_NUM, ROOT, __location__,
        segment_infile=SEGMENT_INFILE,
        guide_star_params_dict=guide_star_params_dict,
        selected_segs=SELECTED_SEGS
    )

    final_num_list, obs_list_string = sg._split_obs_num(obs_num)

    assert final_num_list == correct_list, 'Incorrect observation number list'
    assert obs_list_string == correct_string, 'Incorrect observation number string'

test_data = PARAMETRIZED_DATA['test_segment_override_file_multiple_obs']
sof_multiple_obs_parameters = [('2', 'gs_override_1141_2_1.txt', test_data[0]),
                               ('2, 4, 6', 'gs_override_1141_2,4,6_1.txt', test_data[1]),
                               ('13-15', 'gs_override_1141_13-15_1.txt', test_data[2]),
                               ('5, 7, 10-12', 'gs_override_1141_5,7,10-12_1.txt', test_data[3]),
                               ('5, 3, 1', 'gs_override_1141_1,3,5_1.txt', test_data[4]),
                               ('5, 3, 1, 5', 'gs_override_1141_1,3,5_1.txt', test_data[4]),
                               ('10-12, 5, 7', 'gs_override_1141_5,7,10-12_1.txt', test_data[5])]
@pytest.mark.parametrize('obs_num, correct_file_name, correct_command', sof_multiple_obs_parameters)
def test_segment_override_file_multiple_obs(test_directory, obs_num, correct_file_name, correct_command):
    guider = 1
    guide_star_params_dict = {'v2_boff': 0.1,
                              'v3_boff': 0.2,
                              'fgs_num': guider,
                              'ra': 90.9708,
                              'dec': -67.3578,
                              'pa': 157.1234,
                              'seg_num': 0}

    generate_segment_override_file(
        SEGMENT_INFILE, guider, PROGRAM_ID, obs_num, 1, root=ROOT,
        out_dir=__location__, selected_segs=SELECTED_SEGS, click_to_select_gui=False,
        guide_star_params_dict=guide_star_params_dict, parameter_dialog=False
    )

    # Check to make sure the override file was created, and in the right place
    segment_override_file = os.path.join(
        test_directory, '{}_{}'.format(datetime.now().strftime('%Y%m%d'), correct_file_name)
    )
    assert os.path.isfile(segment_override_file)

    # Also check for the report
    segment_override_report = os.path.join(
        test_directory, '{}_{}'.format(datetime.now().strftime('%Y%m%d'),
                                       correct_file_name.split('.txt')[0] + '_REPORT.txt')
    )
    assert os.path.isfile(segment_override_report)

    # # Check to make sure the command was written correctly
    with open(segment_override_file) as f:
        segment_override_command = f.read()
    assert segment_override_command == correct_command


test_data = PARAMETRIZED_DATA['test_photometry_override_file_multiple_obs']
pof_multiple_obs_parameters = [('2', 'gs_override_1141_2_1.txt', test_data[0]),
                               ('2, 4, 6', 'gs_override_1141_2,4,6_1.txt', test_data[1]),
                               ('13-15', 'gs_override_1141_13-15_1.txt', test_data[2]),
                               ('5, 7, 10-12', 'gs_override_1141_5,7,10-12_1.txt', test_data[3]),
                               ('5, 3, 1', 'gs_override_1141_1,3,5_1.txt', test_data[4]),
                               ('5, 3, 1, 5', 'gs_override_1141_1,3,5_1.txt', test_data[4]),
                               ('10-12, 5, 7', 'gs_override_1141_5,7,10-12_1.txt', test_data[5])]
@pytest.mark.parametrize('obs_num, correct_file_name, correct_command', pof_multiple_obs_parameters)
def test_photometry_override_file_multiple_obs(test_directory, obs_num, correct_file_name, correct_command):

    generate_photometry_override_file(
        ROOT, PROGRAM_ID, obs_num, VISIT_NUM, countrate_factor=0.7,
        out_dir=__location__, parameter_dialog=False
    )

    # Check to make sure the override file was created, and in the right place
    photometry_override_file = os.path.join(
        test_directory, '{}_{}'.format(datetime.now().strftime('%Y%m%d'), correct_file_name)
    )
    assert os.path.isfile(photometry_override_file)

    # # Check to make sure the command was written correctly
    with open(photometry_override_file) as f:
        photometry_override_command = f.read()
    assert photometry_override_command == correct_command
