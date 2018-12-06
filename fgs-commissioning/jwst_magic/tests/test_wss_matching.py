"""Collection of unit tests to verify the correct function of the WSS segment numbering
 matching to PSFs.

TODO


Authors
-------
    - Keira Brooks
    - Lauren Chambers

Use
---
    ::
        pytest test_wss_matching.py

Notes
-----

The following are ways to run segment_guiding:

    from jwst_magic.match_to_wss import MatchToWss
    segment_guiding.run_tool(segment_infile=SEGMENT_INFILE, guider=guider)


"""
import os
import shutil


import numpy as np
from PyQt5 import QtCore
from PyQt5.QtWidgets import QDialogButtonBox

from jwst_magic.match_to_wss import MatchToWss

import pytest

__location__ = os.path.realpath(os.path.join(os.getcwd(), os.path.dirname(__file__)))
FGS_IM_GA = os.path.join(__location__, 'data', 'fgs_data_1_beforelos2ga.fits')
FGS_IM_CMIMF = os.path.join(__location__, 'data', 'fgs_data_2_cmimf.fits')

# https://docs.pytest.org/en/latest/getting-started.html#group-multiple-tests-in-a-class

def test_wss_segs_dict():
    obj = MatchToWss(FGS_IM_GA, global_alignment=True)
    obj.dictionary = {}
    wss_seg_dict = obj.create_wss_seg_dict()






sof_valueerror_parameters = [(20, 1, 'out of range'),
                             ('zero', 2, 'invalid literal for int()'),
                             (0, 3, 'Invalid guider number')]
@pytest.mark.parametrize('seg_num, guider, error_text', sof_valueerror_parameters)
def test_segment_guiding_calculator_valueerrors(seg_num, guider, error_text):
    # Make sure the test output directory does not already exist
    test_output_directory = os.path.join(__location__, 'out', ROOT)
    if os.path.isdir(test_output_directory):
        shutil.rmtree(test_output_directory)

    # Define the input file locations and parameters
    guide_star_params_dict = {'v2_boff': 0.1,
                              'v3_boff': 0.2,
                              'fgs_num': guider,
                              'ra': 90.9708,
                              'dec': -67.3578,
                              'pa': 157.1234,
                              'seg_num': seg_num}

    with pytest.raises(ValueError) as excinfo:
        segment_guiding.generate_segment_override_file(
            SEGMENT_INFILE, guider, PROGRAM_ID, OBSERVATION_NUM, VISIT_NUM,
            root=ROOT, out_dir=__location__, selected_segs=SELECTED_SEGS,
            click_to_select_gui=False, guide_star_params_dict=guide_star_params_dict,
            parameter_dialog=False
        )
    assert error_text in str(excinfo)

    # Remove the test_sgt directory
    shutil.rmtree(test_output_directory)

def test_generate_override_file_valueerrors():
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
        segment_guiding.generate_segment_override_file(
            SEGMENT_INFILE, guider, PROGRAM_ID, OBSERVATION_NUM, VISIT_NUM,
            root=ROOT, out_dir=__location__, guide_star_params_dict=None,
            parameter_dialog=False
        )
    assert '`parameter_dialog=False`' in str(excinfo)

    # Test error if click_to_select_gui=True and data=None
    with pytest.raises(ValueError) as excinfo:
        segment_guiding.generate_segment_override_file(
            SEGMENT_INFILE, guider, PROGRAM_ID, OBSERVATION_NUM, VISIT_NUM,
            root=ROOT, out_dir=__location__,
            click_to_select_gui=True, guide_star_params_dict=guide_star_params_dict,
            parameter_dialog=False
        )
    assert '`click_to_select_GUI=True`' in str(excinfo)

    # Test error if click_to_select_gui=False and selected_segs=None
    with pytest.raises(ValueError) as excinfo:
        segment_guiding.generate_segment_override_file(
            SEGMENT_INFILE, guider, PROGRAM_ID, OBSERVATION_NUM, VISIT_NUM,
            root=ROOT, out_dir=__location__,
            click_to_select_gui=False, guide_star_params_dict=guide_star_params_dict,
            parameter_dialog=False
        )
    assert '`click_to_select_GUI=False`' in str(excinfo)

    # Test error if parameter_dialog=False and countrate_factor=None
    with pytest.raises(ValueError) as excinfo:
        segment_guiding.generate_photometry_override_file(
            ROOT, PROGRAM_ID, OBSERVATION_NUM, VISIT_NUM,
            out_dir=__location__, parameter_dialog=False
        )
    assert '`parameter_dialog=False`' in str(excinfo)


def test_segment_override_command_out_of_fov(caplog, capsys):
    # Define the input file locations and parameters
    guide_star_params_dict = {'v2_boff': 0.1,
                              'v3_boff': 0.2,
                              'fgs_num': 1,
                              'ra': 90.9708,
                              'dec': -67.3578,
                              'pa': 157.1234,
                              'seg_num': 0}

    sg = segment_guiding.SegmentGuidingCalculator(
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
    assert 'Incorrect segment guiding calculations' in str(excinfo)

    # # Check with phony X/Y locations between -5000 and 5000
    # sg.x_det_segs = np.random.rand(18) * 10000 - 5000
    # sg.y_det_segs = np.random.rand(18) * 10000 - 5000
    # with pytest.raises(ValueError):
    #     sg.check_segments_inside_fov(attitude)
    # # Check that the log is raising a warning
    # assert 'off detector' in caplog.text


def test_generate_photometry_override_file():
    # Make sure the test output directory does not already exist
    test_output_directory = os.path.join(__location__, 'out', ROOT)
    if os.path.isdir(test_output_directory):
        shutil.rmtree(test_output_directory)

    segment_guiding.generate_photometry_override_file(
        ROOT, PROGRAM_ID, OBSERVATION_NUM, VISIT_NUM, countrate_factor=0.7,
        out_dir=__location__, parameter_dialog=False
    )

    # Check to make sure the override file was created, and in the right place
    segment_override_file = os.path.join(__location__, 'out', ROOT, 'gs-override-1141_7_1.txt')
    assert os.path.isfile(segment_override_file)

    # Check to make sure the command was written correctly
    with open(segment_override_file) as f:
        segment_override_command = f.read()
    assert segment_override_command == 'sts -gs_select 1141:7:1 -count_rate_factor=0.700'

    # Try again with an incorrect countrate factor and make sure there's an error
    with pytest.raises(ValueError) as excinfo:
        segment_guiding.generate_photometry_override_file(
            ROOT, PROGRAM_ID, OBSERVATION_NUM, VISIT_NUM, countrate_factor=2.0,
            out_dir=__location__, parameter_dialog=False)
    assert 'for count_rate_factor. Expecting between 0.0 and 1.0.' in str(excinfo)

    # Remove the test_sgt directory
    shutil.rmtree(test_output_directory)

def test_convert_boresight_to_v2v3():
    v2v3_offset = segment_guiding._convert_nrca3pixel_offset_to_v2v3_offset(-20.4, -140.53)
    assert v2v3_offset == (-0.638167896, -4.4203823924000005)

def test_cancel_file_dialog():
    """Raise a segment override file dialog window and cancel it.
    """
    # Initialize dialog window
    segment_guiding_dialog = segment_guiding.SegmentGuidingDialog("SOF", 1, PROGRAM_ID, OBSERVATION_NUM, VISIT_NUM)

    # Schedule press of "Cancel" button
    cancel_button = segment_guiding_dialog.buttonBox.button(QDialogButtonBox.Cancel)
    QtCore.QTimer.singleShot(0, cancel_button.clicked)

    # Run GUI
    accepted = segment_guiding_dialog.exec()

    assert not accepted

def test_SOF_parameters_dialog():
    # Initialize dialog window
    segment_guiding_dialog = segment_guiding.SegmentGuidingDialog("SOF", 1, PROGRAM_ID, OBSERVATION_NUM, VISIT_NUM)

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
        '1141', '7', '1', 0.9, None
    )

def test_POF_parameters_dialog():
    # Initialize dialog window
    segment_guiding_dialog = segment_guiding.SegmentGuidingDialog("POF", 1, PROGRAM_ID, OBSERVATION_NUM, VISIT_NUM)

    # Set parameters
    segment_guiding_dialog.lineEdit_programNumber.setText('1142')
    segment_guiding_dialog.lineEdit_observationNumber.setText('8')
    segment_guiding_dialog.lineEdit_visitNumber.setText('2')

    # Schedule press of "Ok" button
    ok_button = segment_guiding_dialog.buttonBox.button(QDialogButtonBox.Ok)
    QtCore.QTimer.singleShot(0, ok_button.clicked)

    # Run GUI
    segment_guiding_dialog.exec()

    # Get parameters
    params = segment_guiding_dialog.get_dialog_parameters()

    assert params == (None, '1142', '8', '2', None, 0.0)

# def test_dialog_window():




# ROOT = "data_beforelos2"
# GUIDER = 1
# SEGMENT_INFILE = "data_beforelos2_G1_ALLpsfs.txt"
# SELECTED_SEGS = "data_beforelos2_G1_regfile.txt"
# GUIDE_STAR_PARAMS_DICT = {'v2_boff': 0.0,
#                           'v3_boff': 0.0,
#                           'fgs_num': 1,
#                           'ra': 55,
#                           'dec': 7,
#                           'pa': 158.,
#                           'seg_num': 0}
#
# PROGRAM_ID = "1326"
# OBSERVATION_NUM = "1"
# VISIT_NUM = "1"
#
#
# def test_guide_star_params_no_ra_dec():
#     segment_guiding.run_tool(SEGMENT_INFILE, GUIDER, program_id=PROGRAM_ID,
#                              observation_num=OBSERVATION_NUM,
#                              visit_num=VISIT_NUM,
#                              click_to_select_GUI=False,
#                              selected_segs=SELECTED_SEGS,
#                              refonly=True,
#                              guide_star_params_dict=GUIDE_STAR_PARAMS_DICT)
#
# def test2():
#     # Parameter dialog off
#     # Params dict with sra, dec, and/or pa set to None
#     guide_star_params_dict = {'v2_boff': 0.0,
#                               'v3_boff': 0.0,
#                               'fgs_num': 1,
#                               'ra': None,
#                               'dec': None,
#                               'pa': None,
#                               'seg_num': 0}
#     segment_guiding.run_tool(SEGMENT_INFILE, GUIDER,
#                              program_id=PROGRAM_ID,
#                              observation_num=OBSERVATION_NUM,
#                              visit_num=VISIT_NUM,
#                              click_to_select_GUI=False,
#                              selected_segs=SELECTED_SEGS,
#                              refonly=True,
#                              guide_star_params_dict=guide_star_params_dict,
#                              parameter_dialog=False)
#
# def test3():
#     # Parameter dialog off
#     # Params dict with sra, dec, and/or pa set to floats
#     guide_star_params_dict = {'v2_boff': 0.0,
#                               'v3_boff': 0.0,
#                               'fgs_num': 1,
#                               'ra': 55.0,
#                               'dec': 7.0,
#                               'pa': 155.0,
#                               'seg_num': 0}
#     segment_guiding.run_tool(SEGMENT_INFILE, GUIDER,
#                              program_id=PROGRAM_ID,
#                              observation_num=OBSERVATION_NUM,
#                              visit_num=VISIT_NUM,
#                              click_to_select_GUI=False,
#                              selected_segs=SELECTED_SEGS,
#                              refonly=True,
#                              guide_star_params_dict=guide_star_params_dict,
#                              parameter_dialog=False)
#
#
# def test4():
#     # Parameter dialog off
#     # Params dict with some ints
#     guide_star_params_dict = {'v2_boff': 0.0,
#                               'v3_boff': 0.0,
#                               'fgs_num': 1,
#                               'ra': 55,
#                               'dec': 7,
#                               'pa': 158.,
#                               'seg_num': 0}
#     segment_guiding.run_tool(SEGMENT_INFILE, GUIDER,
#                              program_id=PROGRAM_ID,
#                              observation_num=OBSERVATION_NUM,
#                              visit_num=VISIT_NUM,
#                              click_to_select_GUI=False,
#                              selected_segs=SELECTED_SEGS,
#                              refonly=True,
#                              guide_star_params_dict=guide_star_params_dict,
#                              parameter_dialog=False)
#
# def test5():
#     # Parameter dialog off
#     # Params dict with sra, dec, and/or pa out of range
#     guide_star_params_dict = {'v2_boff': 0.0,
#                               'v3_boff': 0.0,
#                               'fgs_num': 1,
#                               'ra': 55,
#                               'dec': -91,
#                               'pa': 158.,
#                               'seg_num': 0}
#     segment_guiding.run_tool(SEGMENT_INFILE, GUIDER,
#                              program_id=PROGRAM_ID,
#                              observation_num=OBSERVATION_NUM,
#                              visit_num=VISIT_NUM,
#                              click_to_select_GUI=False,
#                              selected_segs=SELECTED_SEGS,
#                              refonly=True,
#                              guide_star_params_dict=guide_star_params_dict,
#                              parameter_dialog=False)
