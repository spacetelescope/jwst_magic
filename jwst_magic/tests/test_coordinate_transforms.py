"""Collection of unit tests to verify the correct function of the
utils.coordinate_transforms module.

Authors
-------
    - Lauren Chambers

Use
---
    ::
        pytest test_coordinate_transforms.py
"""
import numpy as np
import pytest

from jwst_magic.utils import coordinate_transforms

PIXEL_COORDS = (1743.3, 241.9)
ANGLE_COORDS = (-55.736935, 8.139518)


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
