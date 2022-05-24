"""Collection of unit tests to verify the correct function of the
utils.coordinate_transforms module.

Authors
-------
    - Lauren Chambers
    - Shannon Osborne

Use
---
    ::
        pytest test_coordinate_transforms.py
"""
import numpy as np
import pytest

from jwst_magic.tests.utils import parametrized_data
from jwst_magic.utils import coordinate_transforms

PIXEL_COORDS = (1743.3, 241.9)
ANGLE_COORDS = (-55.736935, 8.139518)
PARAMETRIZED_DATA = parametrized_data()['test_segment_guiding']

def test_idl2dhas():
    x_dhas, y_dhas = coordinate_transforms.idl2dhas(*ANGLE_COORDS)
    assert (x_dhas, y_dhas) == (55.736935, 8.139518), \
        "Incorrect conversion from ideal to DHAS coordinates."


raw2dhas_parameters = [(1, (-53.563954655999986, -50.428367427999994)),
                       (2, (-53.149448727999996, 50.25679589999999))]
@pytest.mark.parametrize('guider, correct_conversion', raw2dhas_parameters)
def test_raw2dhas(guider, correct_conversion):
    x_dhas, y_dhas = coordinate_transforms.raw2dhas(*PIXEL_COORDS, guider)
    assert np.isclose(x_dhas, correct_conversion[0]), \
        "Incorrect conversion from FGS raw to DHAS coordinates."
    assert np.isclose(y_dhas, correct_conversion[1]), \
        "Incorrect conversion from FGS raw to DHAS coordinates."


raw2idl_parameters = [(1, (53.563954655999986, -50.428367427999994)),
                      (2, (53.149448727999996, 50.25679589999999))]
@pytest.mark.parametrize('guider, correct_conversion', raw2idl_parameters)
def test_raw2idl(guider, correct_conversion):
    x_idealangle, y_idealangle = coordinate_transforms.raw2idl(*PIXEL_COORDS, guider)
    assert np.isclose(x_idealangle, correct_conversion[0]), \
        "Incorrect conversion from FGS raw to ideal coordinates."
    assert np.isclose(y_idealangle, correct_conversion[1]), \
        "Incorrect conversion from FGS raw to ideal coordinates."


raw2tel_parameters = [(1, (154.01361684103352, -749.5556063870631)),
                      (2, (-30.13725333672268, -649.0559724725167))]
@pytest.mark.parametrize('guider, correct_conversion', raw2tel_parameters)
def test_raw2tel(guider, correct_conversion):
    v2, v3 = coordinate_transforms.raw2tel(*PIXEL_COORDS, guider)
    assert np.isclose(v2, correct_conversion[0]), \
        "Incorrect conversion from FGS raw to V2/V3 coordinates."
    assert np.isclose(v3, correct_conversion[1]), \
        "Incorrect conversion from FGS raw to V2/V3 coordinates."


def test_convert_boresight_to_v2v3():
    v2v3_offset = coordinate_transforms.nrcpixel_offset_to_v2v3_offset(-20.4, -140.53, 'NRCA3')

    assert np.isclose(v2v3_offset[0], -0.638167896), \
        "V2 value does not match expected value for boresight offset."
    assert np.isclose(v2v3_offset[1], -4.4203823924000005), \
        "V3 value does not match expected value for boresight offset."


def test_convert_sky_to_idl_function():
    """Test the conversion of RA and Dec to Idl X and Y
    as done in the segment guiding code for the guide
    star REPORT.txt file
    """
    ra_seg = np.array([71.329352, 71.329352, 71.329352, 71.368611,
                       71.390103, 71.350521, 71.382983, 71.358393])
    dec_seg = np.array([-78.62446, -78.62446, -78.62446, -78.614016,
                        -78.619765, -78.61443, -78.62308, -78.626784])
    pa = 323

    idl_x, idl_y = coordinate_transforms.convert_sky_to_idl(ra_seg[0], dec_seg[0], pa,
                                            ra_list=ra_seg, dec_list=dec_seg,
                                            guider=1, oss=False)

    # Check that the 3 guide stars have the equivalent position to (0,0)
    np.testing.assert_array_almost_equal(idl_x[0:3], np.zeros(3), decimal=6)
    np.testing.assert_array_almost_equal(idl_y[0:3], np.zeros(3), decimal=6)

    truth_x = PARAMETRIZED_DATA['test_convert_sky_to_idl_function'][0]
    truth_y = PARAMETRIZED_DATA['test_convert_sky_to_idl_function'][1]

    np.testing.assert_array_almost_equal(idl_x, truth_x, decimal=6)
    np.testing.assert_array_almost_equal(idl_y, truth_y, decimal=6)
