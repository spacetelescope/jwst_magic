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


DHAS2Raw_parameters = [(1, (928, 206)),
                       (2, (1156, 207))]
@pytest.mark.parametrize('guider, correct_conversion', DHAS2Raw_parameters)
def test_DHAS2Raw(guider, correct_conversion):
    x_raw, y_raw = coordinate_transforms.DHAS2Raw(*ANGLE_COORDS, guider)
    assert (x_raw, y_raw) == correct_conversion, \
        "Incorrect conversion from DHAS to FGS raw coordinates."


def test_Idl2DHAS():
    x_dhas, y_dhas = coordinate_transforms.Idl2DHAS(*ANGLE_COORDS)
    assert (x_dhas, y_dhas) == (55.736935, 8.139518), \
        "Incorrect conversion from ideal to DHAS coordinates."


Idl2Raw_parameters = [(1, (882, 1833)),
                      (2, (1130, 1847))]
@pytest.mark.parametrize('guider, correct_conversion', Idl2Raw_parameters)
def test_Idl2Raw(guider, correct_conversion):
    x_raw, y_raw = coordinate_transforms.Idl2Raw(*ANGLE_COORDS, guider)
    assert (x_raw, y_raw) == correct_conversion, \
        "Incorrect conversion from ideal to FGS raw coordinates."


def test_Raw2Det():
    x_det, y_det = coordinate_transforms.Raw2Det(*PIXEL_COORDS)
    assert (x_det, y_det) == (241.9, 1743.3), \
        "Incorrect conversion from FGS raw to FGS detector coordinates."


Raw2DHAS_parameters = [(1, (-53.70147089878638, -49.332391547578155)),
                       (2, (-53.14152903705431, 48.846349592719555))]
@pytest.mark.parametrize('guider, correct_conversion', Raw2DHAS_parameters)
def test_Raw2DHAS(guider, correct_conversion):
    x_dhas, y_dhas = coordinate_transforms.Raw2DHAS(*PIXEL_COORDS, guider)
    assert np.isclose(x_dhas, correct_conversion[0]), \
        "Incorrect conversion from FGS raw to DHAS coordinates."
    assert np.isclose(y_dhas, correct_conversion[1]), \
        "Incorrect conversion from FGS raw to DHAS coordinates."


Raw2Idl_parameters = [(1, (53.70147089878638, -49.332391547578155)),
                      (2, (53.14152903705431, 48.846349592719555))]
@pytest.mark.parametrize('guider, correct_conversion', Raw2Idl_parameters)
def test_Raw2Idl(guider, correct_conversion):
    x_idealangle, y_idealangle = coordinate_transforms.Raw2Idl(*PIXEL_COORDS, guider)
    assert np.isclose(x_idealangle, correct_conversion[0]), \
        "Incorrect conversion from FGS raw to ideal coordinates."
    assert np.isclose(y_idealangle, correct_conversion[1]), \
        "Incorrect conversion from FGS raw to ideal coordinates."


Raw2Tel_parameters = [(1, (-52.6118059444085, -50.49287779246356)),
                      (2, (-52.97899983353878, 49.02258207773389))]
@pytest.mark.parametrize('guider, correct_conversion', Raw2Tel_parameters)
def test_Raw2Tel(guider, correct_conversion):
    v2, v3 = coordinate_transforms.Raw2Tel(*PIXEL_COORDS, guider)
    assert np.isclose(v2, correct_conversion[0]), \
        "Incorrect conversion from FGS raw to V2/V3 coordinates."
    assert np.isclose(v3, correct_conversion[1]), \
        "Incorrect conversion from FGS raw to V2/V3 coordinates."


def test_convert_boresight_to_v2v3():
    v2v3_offset = coordinate_transforms.nrca3pixel_offset_to_v2v3_offset(-20.4, -140.53)

    assert np.isclose(v2v3_offset[0], -0.638167896), \
        "V2 value does not match expected value for boresight offset."
    assert np.isclose(v2v3_offset[1], -4.4203823924000005), \
        "V3 value does not match expected value for boresight offset."
