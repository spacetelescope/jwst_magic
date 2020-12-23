""" THE MODULE THIS TEST SUITE TESTS HAS BEEN DEPRECATED AS OF 03 JULY, 2019

Collection of unit tests to verify the correct function of the
convert_image.renormalize module.

Authors
-------
    - Lauren Chambers

Use
---
    ::
        pytest test_renormalize.py
"""

import numpy as np
import pytest

from ..convert_image.renormalize import check_norm_value_unit

gsid_parameters = [
    (12, 'FGS Magnitude', 'N13A000158'),
    ('N13I000018', 'Guide Star ID', 'N13I000018')
]
@pytest.mark.parametrize('norm_value, norm_unit, correct_gsid', gsid_parameters)
def test_check_gsid(norm_value, norm_unit, correct_gsid):
    """
    Test how the normalization is done with  in terms of the interface with the
    fgscountrate module
    """
    gsid = check_norm_value_unit(norm_value, norm_unit)

    assert gsid == correct_gsid

#
# from ..convert_image.renormalize import NormalizeToCountrate
#
# norm_parameters = [
#     (12.5, 'FGS Magnitude', 1, 7575858.324193505, 12.5),
#     (2000000, 'FGS countrate', 1, 2000000, 13.94600462171104),
#     (12.5, 'FGS Magnitude', 2, 8396174.193692164, 12.5),
#     (2000000, 'FGS countrate', 2, 2000000, 14.057628611394456)
# ]
# @pytest.mark.parametrize('value, unit, guider, correct_countrate, correct_mag', norm_parameters)
# def test_NormalizeToCountrate(value, unit, guider, correct_countrate, correct_mag):
#     """Test that the NormalizeToCountrate class can be instantiated,
#     and that the conversion methods between count rate and magnitude
#     is working as expected.
#     """
#     ntc = NormalizeToCountrate(value, unit, guider)
#
#     assert np.isclose(ntc.to_countrate(), correct_countrate), 'Incorrect FGS countrate.'
#     assert np.isclose(ntc.to_fgs_mag(), correct_mag), 'Incorrect FGS magnitude.'
#
#
# def test_NormalizeToCountrate_unit_error():
#     """Test that the NormalizeToCountrate class raises the desired
#     errors when an invalid unit is specified.
#     """
#     ntc = NormalizeToCountrate(12.0, "banana", 1)
#
#     with pytest.raises(ValueError) as excinfo:
#         ntc.to_countrate()
#     assert 'Unknown unit' in str(excinfo.value)
#
#     with pytest.raises(ValueError) as excinfo:
#         ntc.to_fgs_mag()
#     assert 'Unknown unit' in str(excinfo.value)
