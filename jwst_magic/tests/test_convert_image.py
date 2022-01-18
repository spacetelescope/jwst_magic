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
import shutil
import yaml

from astropy.io import ascii as asc
from astropy.io import fits
import numpy as np
import pytest

from jwst_magic.tests.utils import parametrized_data
from jwst_magic.convert_image import convert_image_to_raw_fgs
from jwst_magic.utils import utils


__location__ = os.path.realpath(os.path.join(os.getcwd(), os.path.dirname(__file__)))
ROOT = "test_convertim"
NIRCAM_IM = os.path.join(__location__, 'data', 'nircam_data_1_ga.fits')
NIRCAM_PED_IM = os.path.join(__location__, 'data', 'nircam_w_ped.fits')
NIRCAM_MIMF_IM = os.path.join(__location__, 'data', 'nircam_mimf_cal.fits')
FGS_PED_IM = os.path.join(__location__, 'data', 'fgs_w_ped.fits')
FGS_GA_IM = os.path.join(__location__, 'data', 'fgs_data_1_beforelos2ga.fits')
FGS_CMIMF_IM = os.path.join(__location__, 'data', 'fgs_data_2_cmimf.fits')
SEGMENT_INFILE = os.path.join(__location__, 'data', 'all_found_psfs_test_select_psfs.txt')

PARAMETRIZED_DATA = parametrized_data()['test_convert_image']
TEST_DIRECTORY = os.path.join(__location__, 'out', ROOT)


@pytest.fixture()
def test_directory(test_dir=TEST_DIRECTORY):
    """Create a test directory for permission management.

    Parameters
    ----------
    test_dir : str
        Path to directory used for testing

    Yields
    -------
    test_dir : str
        Path to directory used for testing
    """
    utils.ensure_dir_exists(test_dir)  # creates directory with default mode=511

    yield test_dir
    print("teardown test directory")
    if os.path.isdir(test_dir):
        shutil.rmtree(test_dir)


norm_parameters = [
    (NIRCAM_IM, 1, True, 2000000, 'FGS countrate', True, 2052.),
    (NIRCAM_IM, 2, True, 12, 'FGS Magnitude', True, 6445.),
    (NIRCAM_IM, 1, True, 'N13I000018', 'Guide Star ID', True, 1812.),
    (FGS_GA_IM, 2, False, 12, 'FGS Magnitude', True, 5770.),
    (NIRCAM_IM, 2, True, '', 'Guide Star ID', True, 6445.),  # uses fgs_mag = 12 by default
    (NIRCAM_PED_IM, 1, True, 12, 'FGS Magnitude', False, 141800.), # NRC contains TEST keyword - tests ped
    (FGS_PED_IM, 1, False, 12, 'FGS Magnitude', False, 163398.), # non-ITM FGS image - tests ped
]
@pytest.mark.parametrize('image, guider, nircam, norm_value, norm_unit, itm, data_max', norm_parameters)
def test_convert_im_normalization(test_directory, image, guider, nircam, norm_value, norm_unit, itm, data_max):
    """
    Test how the normalization is done with  in terms of the interface with the
    fgscountrate module
    """
    data, all_found_psfs_file, \
    psf_center_file, _ = convert_image_to_raw_fgs.convert_im(image, guider, ROOT, nircam=nircam,
                                                             out_dir=__location__, smoothing='default',
                                                             nircam_det=None, normalize=True,
                                                             norm_value=norm_value,
                                                             norm_unit=norm_unit,
                                                             gs_catalog='GSC242',
                                                             coarse_pointing=False,
                                                             jitter_rate_arcsec=None,
                                                             logger_passed=False, itm=itm)

    assert os.path.exists(all_found_psfs_file)
    assert psf_center_file is None
    assert np.isclose(data.max(), data_max, atol=1.)


norm_parameters = [
    (12.5, 'Guide Star ID', TypeError, 'Mismatch: Normalization value for unit of'),
    ('N13I000018', 'FGS Magnitude', TypeError, 'Mismatch: Normalization value for unit of'),
    ('N13I000018', 'FGS countrate', TypeError, 'Mismatch: Normalization value for unit of')
]
@pytest.mark.parametrize('value, unit, error, error_text', norm_parameters)
def test_convert_im_error(test_directory, value, unit, error, error_text):
    with pytest.raises(error) as excinfo:
        data, _, _, _ = convert_image_to_raw_fgs.convert_im(NIRCAM_IM, 1, ROOT, nircam=True,
                                                            out_dir=__location__,
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
    sources = utils.find_peaks(image, box_size=5, threshold=np.max(image) * .05)

    coords = [(x, y) for (x, y) in sources['x_peak', 'y_peak']]
    assert np.array_equal(np.array(coords), test_data[0]), \
        'Incorrect transformation from DMS/science frame to raw FGS1 frame'

    # Transform sci to FGS2 raw
    image = convert_image_to_raw_fgs.transform_sci_to_fgs_raw(data, 2)
    sources = utils.find_peaks(image, box_size=5, threshold=np.max(image) * .05)
    coords = [(x, y) for (x, y) in sources['x_peak', 'y_peak']]
    assert np.array_equal(np.array(coords), test_data[1]), \
        'Incorrect transformation from DMS/science frame to raw FGS2 frame'


def test_transform_nircam_raw_to_fgs_raw():
    test_data = PARAMETRIZED_DATA['test_transform_nircam_raw_to_fgs_raw']

    with fits.open(FGS_CMIMF_IM) as hdulist:
        data = hdulist[1].data

    # Transform NRCA3 raw to FGS1 raw
    image = convert_image_to_raw_fgs.transform_nircam_raw_to_fgs_raw(data, 'A3', 1)
    sources = utils.find_peaks(image, box_size=5, threshold=np.max(image) * .05)
    coords = [(x, y) for (x, y) in sources['x_peak', 'y_peak']]
    assert np.array_equal(np.array(coords), test_data[0]),\
        'Incorrect transformation from raw NRCA3 (thus also A1, A5, B2, B4) frame to raw FGS1 frame'

    # Transform NRCB1 raw to FGS1 raw
    image = convert_image_to_raw_fgs.transform_nircam_raw_to_fgs_raw(data, 'B1', 1)
    sources = utils.find_peaks(image, box_size=5, threshold=np.max(image) * .05)
    coords = [(x, y) for (x, y) in sources['x_peak', 'y_peak']]
    assert np.array_equal(np.array(coords), test_data[1]), \
        'Incorrect transformation from raw NRCB1 (thus also A2, A4, B3, B5) frame to raw FGS1 frame'

    # Transform NRCA5 raw to FGS2 raw
    image = convert_image_to_raw_fgs.transform_nircam_raw_to_fgs_raw(data, 'A5', 2)
    sources = utils.find_peaks(image, box_size=5, threshold=np.max(image) * .05)
    coords = [(x, y) for (x, y) in sources['x_peak', 'y_peak']]
    assert np.array_equal(np.array(coords), test_data[2]), \
        'Incorrect transformation from raw NRCA5 (thus also A1, A3, B2, B4) frame to raw FGS2 frame'

    # Transform NRCB1 raw to FGS2 raw
    image = convert_image_to_raw_fgs.transform_nircam_raw_to_fgs_raw(data, 'B1', 2)
    sources = utils.find_peaks(image, box_size=5, threshold=np.max(image) * .05)
    coords = [(x, y) for (x, y) in sources['x_peak', 'y_peak']]
    assert np.array_equal(np.array(coords), test_data[3]), \
        'Incorrect transformation from raw NRCB1 (thus also A2, A4, B3, B5) frame to raw FGS2 frame'


def test_create_psf_center_file(test_directory):
    """ Test that the right files for the MIMF case are written out and contain the correct data
    """
    image = NIRCAM_MIMF_IM
    guider = 1

    data, all_found_psfs_file, \
    psf_center_file, _ = convert_image_to_raw_fgs.convert_im(image, guider, ROOT, nircam=True,
                                                             out_dir=__location__, smoothing='low',
                                                             nircam_det=None, normalize=True,
                                                             norm_value=12,
                                                             norm_unit='FGS Magnitude',
                                                             gs_catalog='GSC242',
                                                             coarse_pointing=False,
                                                             jitter_rate_arcsec=None,
                                                             logger_passed=False, itm=False)

    assert os.path.exists(all_found_psfs_file)
    assert os.path.exists(psf_center_file)

    # Check PSF location in guiding selections doesn't match psf center
    # guiding selections should have found a knot in the PSF not at the center
    with open(all_found_psfs_file) as f:
        guiding_selections_contents = f.read()
    with open(psf_center_file) as f:
        no_smooth_contents = f.read()

    assert guiding_selections_contents != no_smooth_contents


def test_read_in_all_found_psfs_file(test_directory):
    """Test reading in an all found psfs file into the convert_im function
    """
    image = NIRCAM_IM
    guider = 1
    input_all_found_psfs = SEGMENT_INFILE

    data, all_found_psfs_file, \
    psf_center_file, _ = convert_image_to_raw_fgs.convert_im(image, guider, ROOT, nircam=True,
                                                             out_dir=__location__, smoothing='high',
                                                             nircam_det=None, normalize=True,
                                                             norm_value=12,
                                                             norm_unit='FGS Magnitude',
                                                             all_found_psfs_file=input_all_found_psfs,
                                                             gs_catalog='GSC242',
                                                             coarse_pointing=False,
                                                             jitter_rate_arcsec=None,
                                                             logger_passed=False, itm=False)

    assert os.path.exists(all_found_psfs_file)
    assert psf_center_file is None

    input_table = asc.read(input_all_found_psfs)
    output_table = asc.read(all_found_psfs_file)
    assert all(input_table == output_table)


detection_parameters = [
    ('high', 'standard-deviation', 18),  # high smoothing checked
    ('low', 'standard-deviation', 1),  # no smoothing checked
    ('default', 'pixel-wise', 636),  # pixel-wise threshold checked
    ('high', 'pixel-wise', 28),  # high smoothing and pixel-wise threshold checked
    (26, 'standard-deviation', 18),  # gauss sigma of 26 (same as high smoothing)
]
@pytest.mark.parametrize('smoothing, detection_threshold, correct_number_of_psfs', detection_parameters)
def test_psf_detection_methods(test_directory, smoothing, detection_threshold, correct_number_of_psfs):
    """Test the different options for detecting PSFs"""
    image = NIRCAM_IM
    guider = 1
    data, all_found_psfs_file, \
    psf_center_file, _ = convert_image_to_raw_fgs.convert_im(image, guider, ROOT, nircam=True,
                                                             out_dir=__location__, smoothing=smoothing,
                                                             detection_threshold=detection_threshold,
                                                             nircam_det=None, normalize=True,
                                                             norm_value=12,
                                                             norm_unit='FGS Magnitude',
                                                             all_found_psfs_file=None,
                                                             gs_catalog='GSC242',
                                                             coarse_pointing=False,
                                                             jitter_rate_arcsec=None,
                                                             logger_passed=False, itm=False)

    # Check the saved out all_found_psfs_file
    in_table = asc.read(all_found_psfs_file)
    assert len(in_table) == correct_number_of_psfs
