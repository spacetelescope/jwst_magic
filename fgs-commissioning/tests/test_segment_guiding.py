'''Collection of unit tests to verify the correct function of the Segment
    Guiding tool'''
import glob
import os

import numpy as np
import pytest
from astropy.io import fits

from jwst_magic import utils, coordinate_transforms
from jwst_magic.segment_guiding import segment_guiding





# 1. Test all of the different parameters in run_tool individually.
# 2. Test with the params dict with ra, dec, and/or pa set to None
# 3. Test with the params dict with ra, dec, and/or pa set to floats
# 4. Test with the params dict with ra, dec, and/or pa set to ints
# 5. Test with the params dict with ra, dec, and/or pa set out of range
# 6. Test with parameter_dialog turned on with ra, dec, and/or pa set to None
# 7. Test with parameter_dialog turned on with ra, dec, and/or pa set to floats
# 8. Test with parameter_dialog turned on with ra, dec, and/or pa set to ints
# 9. Test with parameter_dialog turned on with ra, dec, and/or pa set out of range
# 10. Test with program_id, obsercaiotn_num, and/or visit_num set to ints
# 11. Test with program_id, obsercaiotn_num, and/or visit_num set to strings
# 12. Test with program_id, obsercaiotn_num, and/or visit_num set to floats
# 13. Test with segment_infile and/or guide set to None
# 14.
#
#
#
# from jwst_magic.segment_guiding import segment_guiding
# segment_guiding.run_tool(segment_infile=segment_infile, guider=guider)
#
# Or with the segment dialog box:
# ::
# from jwst_magic.segment_guiding import segment_guiding
# segment_guiding.run_tool(program_id=program_id, observation=observation,
#                          visit=visit, parameter_dialog=True)
#
# Or from dictionary of parameters:
# ::
# from jwst_magic.segment_guiding import segment_guiding
# segment_guiding.run_tool(program_id=program_id, observation=observation,
#                          visit=visit, guide_star_params_dict=guide_star_params_dict,
#                          parameter_dialog=False)
#

ROOT = "data_beforelos2"
GUIDER = 1
SEGMENT_INFILE = "data_beforelos2_G1_ALLpsfs.txt"
SELECTED_SEGS = "data_beforelos2_G1_regfile.txt"
GUIDE_STAR_PARAMS_DICT = {'v2_boff': 0.0,
                          'v3_boff': 0.0,
                          'fgs_num': 1,
                          'ra': 55,
                          'dec': 7,
                          'pa': 158.,
                          'seg_num': 0}

PROGRAM_ID = "1326"
OBSERVATION_NUM = "1"
VISIT_NUM = "1"


def test_guide_star_params_no_ra_dec():
    segment_guiding.run_tool(SEGMENT_INFILE, GUIDER, program_id=PROGRAM_ID,
                             observation_num=OBSERVATION_NUM,
                             visit_num=VISIT_NUM,
                             click_to_select_GUI=False,
                             selected_segs=SELECTED_SEGS,
                             refonly=True,
                             guide_star_params_dict=GUIDE_STAR_PARAMS_DICT)

def test2():
    # Parameter dialog off
    # Params dict with sra, dec, and/or pa set to None
    guide_star_params_dict = {'v2_boff': 0.0,
                              'v3_boff': 0.0,
                              'fgs_num': 1,
                              'ra': None,
                              'dec': None,
                              'pa': None,
                              'seg_num': 0}
    segment_guiding.run_tool(SEGMENT_INFILE, GUIDER,
                             program_id=PROGRAM_ID,
                             observation_num=OBSERVATION_NUM,
                             visit_num=VISIT_NUM,
                             click_to_select_GUI=False,
                             selected_segs=SELECTED_SEGS,
                             refonly=True,
                             guide_star_params_dict=guide_star_params_dict,
                             parameter_dialog=False)

def test3():
    # Parameter dialog off
    # Params dict with sra, dec, and/or pa set to floats
    guide_star_params_dict = {'v2_boff': 0.0,
                              'v3_boff': 0.0,
                              'fgs_num': 1,
                              'ra': 55.0,
                              'dec': 7.0,
                              'pa': 155.0,
                              'seg_num': 0}
    segment_guiding.run_tool(SEGMENT_INFILE, GUIDER,
                             program_id=PROGRAM_ID,
                             observation_num=OBSERVATION_NUM,
                             visit_num=VISIT_NUM,
                             click_to_select_GUI=False,
                             selected_segs=SELECTED_SEGS,
                             refonly=True,
                             guide_star_params_dict=guide_star_params_dict,
                             parameter_dialog=False)


def test4():
    # Parameter dialog off
    # Params dict with some ints
    guide_star_params_dict = {'v2_boff': 0.0,
                              'v3_boff': 0.0,
                              'fgs_num': 1,
                              'ra': 55,
                              'dec': 7,
                              'pa': 158.,
                              'seg_num': 0}
    segment_guiding.run_tool(SEGMENT_INFILE, GUIDER,
                             program_id=PROGRAM_ID,
                             observation_num=OBSERVATION_NUM,
                             visit_num=VISIT_NUM,
                             click_to_select_GUI=False,
                             selected_segs=SELECTED_SEGS,
                             refonly=True,
                             guide_star_params_dict=guide_star_params_dict,
                             parameter_dialog=False)

def test5():
    # Parameter dialog off
    # Params dict with sra, dec, and/or pa out of range
    guide_star_params_dict = {'v2_boff': 0.0,
                              'v3_boff': 0.0,
                              'fgs_num': 1,
                              'ra': 55,
                              'dec': -91,
                              'pa': 158.,
                              'seg_num': 0}
    segment_guiding.run_tool(SEGMENT_INFILE, GUIDER,
                             program_id=PROGRAM_ID,
                             observation_num=OBSERVATION_NUM,
                             visit_num=VISIT_NUM,
                             click_to_select_GUI=False,
                             selected_segs=SELECTED_SEGS,
                             refonly=True,
                             guide_star_params_dict=guide_star_params_dict,
                             parameter_dialog=False)
