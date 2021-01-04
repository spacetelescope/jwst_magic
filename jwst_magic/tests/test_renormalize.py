"""

Collection of unit tests to verify the correct function of the
convert_image.renormalize module.

Authors
-------
    - Shannon Osborne
    - Lauren Chambers

Use
---
    ::
        pytest test_renormalize.py
"""

import numpy as np
import pytest

from ..convert_image.renormalize import check_norm_value_unit, convert_to_countrate_fgsmag

gsid_parameters = [
    (11, 'FGS Magnitude', 1, 14292991.482979316, 11.047518323503994),
    (1000000, 'FGS Countrate', 2, 1000000, 13.98190894951814),
    ('N1HL000080', 'Guide Star ID', 2, 13010864.64607799, 11.19577024436347)
]
@pytest.mark.parametrize('norm_value, norm_unit, guider, correct_countrate, correct_magnitude', gsid_parameters)
def test_convert_to_count_rate_mag(norm_value, norm_unit, guider, correct_countrate, correct_magnitude):
    """
    Test how the FGS count rate and magnitude are calculated based on input norm unit and value
    """

    fgs_countrate, fgs_mag = convert_to_countrate_fgsmag(norm_value, norm_unit, guider)

    assert np.isclose(fgs_countrate, correct_countrate, 1e-5)
    assert np.isclose(fgs_mag, correct_magnitude, 1e-5)

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

