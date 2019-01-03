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
from scipy import ndimage

from jwst_magic.match_to_wss import MatchToWss

import pytest

__location__ = os.path.realpath(os.path.join(os.getcwd(), os.path.dirname(__file__)))
FGS_IM_GA_RAW = os.path.join(__location__, 'data', 'fgs_data_1_beforelos2ga_raw.fits')
FGS_IM_CMIMF_RAW = os.path.join(__location__, 'data', 'fgs_data_2_cmimf_raw.fits')
FGS_IM_CMIMF_17_RAW = os.path.join(__location__, 'data', 'fgs_data_3_cmimf_raw.fits') # missing one segment

WSS_SEGS_DICT_GA = {}
WSS_SEGS_DICT_CMIMF = {}
WSS_SEGS_DICT_CMIMF_17 = {}

# https://docs.pytest.org/en/latest/getting-started.html#group-multiple-tests-in-a-class
# def test_Class(object):
#     def __init__():
obj_ga = MatchToWss(FGS_IM_GA_RAW, global_alignment=True)
obj_cmimf = MatchToWss(FGS_IM_CMIMF_RAW, global_alignment=True)
ob_cmimf_17 = MatchToWss(FGS_IM_CMIMF_17_RAW, global_alignment=True)



dict_valueerror_parameters = [(obj_ga.dictionary, WSS_SEGS_DICT_GA),
                              (obj_cmimf.dictionary, WSS_SEGS_DICT_CMIMF),
                              (ob_cmimf_17.dictionary, WSS_SEGS_DICT_CMIMF_17)]
@pytest.mark.parametrize('dictionary', dict_valueerror_parameters)
def test_wss_segs_dict(dictionary, truth_dictionary):
    assert dictionary is truth_dictionary



edge_valueerror_parameters = [(obj_ga.top_y_im, obj_ga.bottom_y_im, obj_ga.left_x_im,
                               obj_ga.right_x_im),
                              (obj_ga.top_final, obj_ga.bottom_final, obj_ga.left_final,
                               obj_ga.right_final)]
@pytest.mark.parametrize('top, bottom, left, right', edge_valueerror_parameters)
def test_define_edges(top, bottom, left, right):
    assert top > bottom
    assert right > left


pupil_valueerror_parameters = [(obj_ga.pupil, obj_ga.top_seg, obj_ga.bottom_seg,
                                obj_ga.left_seg, obj_ga.right_seg),
                               (obj_ga.pupil_scaled, obj_ga.top_seg, obj_ga.bottom_seg,
                                obj_ga.left_seg, obj_ga.right_seg),
                               (obj_ga.fill_pupil, obj_ga.top_seg, obj_ga.bottom_seg,
                                obj_ga.left_seg, obj_ga.right_seg),
                               (obj_ga.matched_pupil, obj_ga.top_seg, obj_ga.bottom_seg,
                                obj_ga.left_seg, obj_ga.right_seg)]
@pytest.mark.parametrize('pupil, top_seg, bottom_seg, left_seg, right_seg', pupil_valueerror_parameters)
def test_pupil_edges(pupil, top_seg, bottom_seg, left_seg, right_seg):
    '''
    Make sure that the pupil image hasn't been interpolated or other
    tweaking that would cause the center_of_mass function to give the
    wrong answers
    '''
    #top_final, bottom_final, left_final, right_final = self.define_edges_of_pupil_mask(pupil)
    top_y, top_x = ndimage.measurements.center_of_mass(pupil == top_seg)
    bottom_y, bottom_x = ndimage.measurements.center_of_mass(pupil == bottom_seg)

    left_y, left_x = ndimage.measurements.center_of_mass(pupil == left_seg)
    right_y, right_x = ndimage.measurements.center_of_mass(pupil == right_seg)


    assert top_x == bottom_x or top_y == bottom_y
    assert left_x == right_x or left_y == right_y
