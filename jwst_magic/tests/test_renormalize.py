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

from jwst_magic.convert_image.renormalize import convert_to_countrate_fgsmag

gsid_parameters = [
    (11, 'FGS Magnitude', 1, 14932429.201492183, 11),
    (1000000, 'FGS Countrate', 2, 1000000, 13.98190894951814),
    ('N1HL000080', 'Guide Star ID', 2, 12834521.151443688, 11.209282584607347)
]
@pytest.mark.parametrize('norm_value, norm_unit, guider, correct_countrate, correct_magnitude', gsid_parameters)
def test_convert_to_count_rate_mag(norm_value, norm_unit, guider, correct_countrate, correct_magnitude):
    """
    Test how the FGS count rate and magnitude are calculated based on input norm unit and value
    """

    fgs_countrate, fgs_mag = convert_to_countrate_fgsmag(norm_value, norm_unit, guider, gs_catalog='GSC242')

    assert np.isclose(fgs_countrate, correct_countrate, 1e-5)
    assert np.isclose(fgs_mag, correct_magnitude, 1e-5)

