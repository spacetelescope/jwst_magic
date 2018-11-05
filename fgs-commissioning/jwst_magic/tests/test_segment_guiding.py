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
"""

import glob
import logging
import os
import shutil
import sys

import numpy as np
import pytest
from astropy.io import ascii as asc
from astropy.io import fits

from jwst_magic import utils, coordinate_transforms
from jwst_magic.segment_guiding import segment_guiding

import pytest

__location__ = os.path.realpath(os.path.join(os.getcwd(), os.path.dirname(__file__)))
ROOT = "test_sgt"
SEGMENT_INFILE = os.path.join(__location__, 'data', '{}_ALLpsfs.txt'.format(ROOT))
SELECTED_SEGS = os.path.join(__location__, 'data', '{}_regfile.txt'.format(ROOT))
PROGRAM_ID = 1141
OBSERVATION_NUM = 7
VISIT_NUM = 1


sof_parameters = [(0, None, 'sts -gs_select 1141:7:1 -star1 = 4, 90.986120, -67.362103, 439821.7, 395839.5, 3, 10 -ref_only1 = 3, 90.980194, -67.352832, 369223.3, 332301.0 -ref_only2 = 10, 90.963969, -67.351986, 397375.0, 357637.5'),
                  (4, None, 'sts -gs_select 1141:7:1 -star1 = 4, 90.970810, -67.357738, 439821.7, 395839.5, 3, 10 -ref_only1 = 3, 90.964892, -67.348467, 369223.3, 332301.0 -ref_only2 = 10, 90.948671, -67.347619, 397375.0, 357637.5'),
                  (0, np.array([[1, 12, 6]]), 'sts -gs_select 1141:7:1 -star1 = 1, 90.987850, -67.355022, 316763.3, 285087.0, 12, 6 -ref_only1 = 12, 90.963123, -67.355541, 532725.0, 479452.5 -ref_only2 = 6, 90.972558, -67.350638, 549830.0, 494847.0'),
                  (0, np.array([[1, 2, 3], [4, 0, 17, 12, 2]]), 'sts -gs_select 1141:7:1 -star1 = 1, 90.987850, -67.355022, 316763.3, 285087.0, 2, 3 -star2 = 4, 90.986120, -67.362103, 439821.7, 395839.5, 0, 17, 12, 2 -ref_only1 = 2, 90.986990, -67.358569, 485390.0, 436851.0 -ref_only2 = 3, 90.980194, -67.352832, 369223.3, 332301.0 -ref_only3 = 12, 90.963123, -67.355541, 532725.0, 479452.5 -ref_only4 = 17, 90.954592, -67.356915, 297443.3, 267699.0 -ref_only5 = 0, 90.953775, -67.360493, 408348.3, 367513.5')]
@pytest.mark.parametrize('seg_num, selected_segs, correct_command', sof_parameters)
def test_generate_segment_override_file(seg_num, selected_segs, correct_command):
    # Make sure the test output directory does not already exist
    test_output_directory = os.path.join(__location__, 'out', ROOT)
    if os.path.isdir(test_output_directory):
        shutil.rmtree(test_output_directory)

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

    segment_guiding.generate_segment_override_file(
        SEGMENT_INFILE, guider, PROGRAM_ID, OBSERVATION_NUM, VISIT_NUM, root=ROOT,
        out_dir=__location__, selected_segs=selected_segs, click_to_select_gui=False,
        guide_star_params_dict=guide_star_params_dict, parameter_dialog=False
    )

    # Check to make sure the override file was created, and in the right place
    segment_override_file = os.path.join(__location__, 'out', ROOT, 'gs-override-1141_7_1.txt')
    assert os.path.isfile(segment_override_file)

    # # Check to make sure the command was written correctly
    with open(segment_override_file) as f:
        segment_override_command = f.read()
    assert segment_override_command == correct_command

    print(segment_override_command)

    # Remove the test_sgt directory
    shutil.rmtree(test_output_directory)


sof_valueerror_parameters = [(20, 1, 'out of range'),
                             ('zero', 2, 'Guide star parameter'),
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

    # Check with phony X/Y locations between -5000 and 5000
    sg.x_det_segs = np.random.rand(18) * 10000 - 5000
    sg.y_det_segs = np.random.rand(18) * 10000 - 5000
    with pytest.raises(ValueError):
        sg.check_segments_inside_fov(attitude)
    # Check that the log is raising a warning
    assert 'off detector' in caplog.text


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
