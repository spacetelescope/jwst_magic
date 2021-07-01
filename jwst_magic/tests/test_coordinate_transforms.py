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


def test_Idl2DHAS():
    x_dhas, y_dhas = coordinate_transforms.Idl2DHAS(*ANGLE_COORDS)
    assert (x_dhas, y_dhas) == (55.736935, 8.139518), \
        "Incorrect conversion from ideal to DHAS coordinates."


Raw2DHAS_parameters = [(1, (-53.63248581599999, -50.358308568)),
                       (2, (-53.21744955799999, 50.186975399999994))]
@pytest.mark.parametrize('guider, correct_conversion', Raw2DHAS_parameters)
def test_Raw2DHAS(guider, correct_conversion):
    x_dhas, y_dhas = coordinate_transforms.Raw2DHAS(*PIXEL_COORDS, guider)
    assert np.isclose(x_dhas, correct_conversion[0]), \
        "Incorrect conversion from FGS raw to DHAS coordinates."
    assert np.isclose(y_dhas, correct_conversion[1]), \
        "Incorrect conversion from FGS raw to DHAS coordinates."


Raw2Idl_parameters = [(1, (53.63248581599999, -50.358308568)),
                      (2, (53.21744955799999, 50.186975399999994))]
@pytest.mark.parametrize('guider, correct_conversion', Raw2Idl_parameters)
def test_Raw2Idl(guider, correct_conversion):
    x_idealangle, y_idealangle = coordinate_transforms.Raw2Idl(*PIXEL_COORDS, guider)
    assert np.isclose(x_idealangle, correct_conversion[0]), \
        "Incorrect conversion from FGS raw to ideal coordinates."
    assert np.isclose(y_idealangle, correct_conversion[1]), \
        "Incorrect conversion from FGS raw to ideal coordinates."


Raw2Tel_parameters = [(1, (153.94357268603113, -749.4870601974866)),
                      (2, (-30.20548566488923, -649.1255667574394))]
@pytest.mark.parametrize('guider, correct_conversion', Raw2Tel_parameters)
def test_Raw2Tel(guider, correct_conversion):
    v2, v3 = coordinate_transforms.Raw2Tel(*PIXEL_COORDS, guider)
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
