"""Collection of unit tests to verify the correct function of the Image
Conversion tool.

TODO
(In the future, this will) Runs the following tests:
    1. Test all of convert_image

Authors
-------
    - Keira Brooks
    - Lauren Chambers

Use
---
    ::
        pytest test_convert_image.py

Notes
-----

At this moment, convert image can take in a NIRCam or FGS image either created
with WebbPSF postage stamps or created with the ITM simulator.

The following are ways to run convert_image:

    from jwst_magic.convert_image import convert_image_to_raw
    convert_image_to_raw.convert_im(input_im, guider, root)

"""

import os
import yaml

from astropy.io import fits
import numpy as np
from photutils import find_peaks
import pytest

from .utils import parametrized_data
from ..convert_image import convert_image_to_raw_fgs


__location__ = os.path.realpath(os.path.join(os.getcwd(), os.path.dirname(__file__)))
ROOT = "test_convertim"
NIRCAM_IM = os.path.join(__location__, 'data', 'nircam_data_1_ga.fits')
NIRCAM_PED_IM = os.path.join(__location__, 'data', 'nircam_w_ped.fits')
FGS_PED_IM = os.path.join(__location__, 'data', 'fgs_w_ped.fits')
FGS_GA_IM = os.path.join(__location__, 'data', 'fgs_data_1_beforelos2ga.fits')
FGS_CMIMF_IM = os.path.join(__location__, 'data', 'fgs_data_2_cmimf.fits')

PARAMETRIZED_DATA = parametrized_data()['test_convert_image']



norm_parameters = [
    (NIRCAM_IM, 1, True, 2000000, 'FGS countrate', True, 1992.7434250045358),
    (NIRCAM_IM, 2, True, 12, 'FGS Magnitude', True, 5563.853153791558),
    (NIRCAM_IM, 1, True, 'N13I000018', 'Guide Star ID', True, 1770.7861951935322),
    (FGS_GA_IM, 2, False, 12, 'FGS Magnitude', True, 171.37823193763356),
    (NIRCAM_IM, 2, True, '', 'Guide Star ID', True, 5563.853153791558),
    (NIRCAM_PED_IM, 1, True, 12, 'FGS Magnitude', False, 2994.680894443115), # NRC contains TEST keyword - tests ped
    (FGS_PED_IM, 1, False, 12, 'FGS Magnitude', False, 136580.37527457761) # non-ITM FGS image - tests ped
]
@pytest.mark.parametrize('image, guider, nircam, norm_value, norm_unit, itm, data_max', norm_parameters)
def test_convert_im_normalization(image, guider, nircam, norm_value, norm_unit, itm, data_max):
    """
    Test how the normalization is done with  in terms of the interface with the
    fgscountrate module
    """
    data = convert_image_to_raw_fgs.convert_im(image, guider, ROOT, nircam=nircam,
                                               nircam_det=None, normalize=True,
                                               norm_value=norm_value,
                                               norm_unit=norm_unit,
                                               coarse_pointing=False,
                                               jitter_rate_arcsec=None,
                                               logger_passed=False, itm=itm)

    assert np.isclose(data.max(), data_max)


norm_parameters = [
    (12.5, 'FGS Magnitude', ValueError, 'Unacceptable FGS Magnitude value'),
    ('N13I000018', 'FGS Magnitude', TypeError, 'Mismatched normalization value'),
    ('N13I000018', 'FGS countrate', TypeError, 'Mismatched normalization value')
]
@pytest.mark.parametrize('value, unit, error, error_text', norm_parameters)
def test_convert_im_normalization_error(value, unit, error, error_text):
    with pytest.raises(error) as excinfo:
        data = convert_image_to_raw_fgs.convert_im(NIRCAM_IM, 1, ROOT, nircam=True,
                                                   nircam_det=None, normalize=True,
                                                   norm_value=value, norm_unit=unit,
                                                   coarse_pointing=False,
                                                   jitter_rate_arcsec=None,
                                                   logger_passed=False, itm=True)
        assert error_text in str(excinfo)


def test_transform_sci_to_fgs_raw():
    test_data = PARAMETRIZED_DATA['test_transform_sci_to_fgs_raw']

    with fits.open(FGS_CMIMF_IM) as hdulist:
        data = hdulist[1].data

    # Transform sci to FGS1 raw
    image = convert_image_to_raw_fgs.transform_sci_to_fgs_raw(data, 1)
    sources = find_peaks(image, np.max(image) * .05, box_size=5)
    coords = [(x, y) for (x, y) in sources['x_peak', 'y_peak']]
    assert np.array_equal(np.array(coords), test_data[0]), \
        'Incorrect transformation from DMS/science frame to raw FGS1 frame'

    # Transform sci to FGS2 raw
    image = convert_image_to_raw_fgs.transform_sci_to_fgs_raw(data, 2)
    sources = find_peaks(image, np.max(image) * .05, box_size=5)
    coords = [(x, y) for (x, y) in sources['x_peak', 'y_peak']]
    assert np.array_equal(np.array(coords), test_data[1]), \
        'Incorrect transformation from DMS/science frame to raw FGS2 frame'


def test_transform_nircam_raw_to_fgs_raw():
    test_data = PARAMETRIZED_DATA['test_transform_nircam_raw_to_fgs_raw']

    with fits.open(FGS_CMIMF_IM) as hdulist:
        data = hdulist[1].data

    # Transform NRCA3 raw to FGS1 raw
    image = convert_image_to_raw_fgs.transform_nircam_raw_to_fgs_raw(data, 'A3', 1)
    sources = find_peaks(image, np.max(image) * .05, box_size=5)
    coords = [(x, y) for (x, y) in sources['x_peak', 'y_peak']]
    assert np.array_equal(np.array(coords), test_data[0]),\
        'Incorrect transformation from raw NRCA3 (thus also A1, A5, B2, B4) frame to raw FGS1 frame'

    # Transform NRCB1 raw to FGS1 raw
    image = convert_image_to_raw_fgs.transform_nircam_raw_to_fgs_raw(data, 'B1', 1)
    sources = find_peaks(image, np.max(image) * .05, box_size=5)
    coords = [(x, y) for (x, y) in sources['x_peak', 'y_peak']]
    assert np.array_equal(np.array(coords), test_data[1]), \
        'Incorrect transformation from raw NRCB1 (thus also A2, A4, B3, B5) frame to raw FGS1 frame'

    # Transform NRCA5 raw to FGS2 raw
    image = convert_image_to_raw_fgs.transform_nircam_raw_to_fgs_raw(data, 'A5', 2)
    sources = find_peaks(image, np.max(image) * .05, box_size=5)
    coords = [(x, y) for (x, y) in sources['x_peak', 'y_peak']]
    assert np.array_equal(np.array(coords), test_data[2]), \
        'Incorrect transformation from raw NRCA5 (thus also A1, A3, B2, B4) frame to raw FGS2 frame'

    # Transform NRCB1 raw to FGS2 raw
    image = convert_image_to_raw_fgs.transform_nircam_raw_to_fgs_raw(data, 'B1', 2)
    sources = find_peaks(image, np.max(image) * .05, box_size=5)
    coords = [(x, y) for (x, y) in sources['x_peak', 'y_peak']]
    assert np.array_equal(np.array(coords), test_data[3]), \
        'Incorrect transformation from raw NRCB1 (thus also A2, A4, B3, B5) frame to raw FGS2 frame'


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
    gsid = convert_image_to_raw_fgs.check_norm_value_unit(norm_value, norm_unit)

    assert gsid == correct_gsid
