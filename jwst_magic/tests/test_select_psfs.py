"""Collection of unit tests to verify the correct function of the
select_psfs module.

Authors
-------
    - Lauren Chambers

Use
---
    ::
        pytest test_select_psfs.py
"""
import os
import shutil

from astropy.io import ascii as asc
import pytest

from .utils import parametrized_data
from ..star_selector.select_psfs import select_psfs
from ..utils import utils

__location__ = os.path.realpath(os.path.join(os.getcwd(), os.path.dirname(__file__)))
FGS_CMIMF_IM = os.path.join(__location__, 'data', 'fgs_data_2_cmimf.fits')
NIRCAM_IM = os.path.join(__location__, 'data', 'nircam_data_1_ga.fits')
CONVERTED_NIRCAM_IM_GA = os.path.join(__location__, 'data', 'converted_nircam_data_1_ga.fits')
CONVERTED_NIRCAM_IM_MIMF = os.path.join(__location__, 'data', 'converted_nircam_data_1_mimf.fits')
ROOT = "test_select_psfs"
SEGMENT_INFILE = os.path.join(__location__, 'data', 'all_found_psfs_{}.txt'.format(ROOT))
SELECTED_SEGS = os.path.join(__location__, 'data', 'guiding_selections_{}.txt'.format(ROOT))
VISIT_NUM = 1

TEST_DIRECTORY = os.path.join(__location__, 'out', ROOT)

PARAMETRIZED_DATA = parametrized_data()['test_select_psfs']


@pytest.fixture(scope="module")
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


def test_select_psfs_with_file(test_directory):
    guiding_selections_path, all_found_psfs_path = select_psfs(
        CONVERTED_NIRCAM_IM_GA, ROOT, 2, guiding_selections_file=SELECTED_SEGS,
        out_dir=test_directory
    )

    # Ensure the correct files were generated
    assert os.path.exists(guiding_selections_path), 'guiding_selections_test_select_psfs_G2.txt not generated.'
    assert os.path.exists(all_found_psfs_path), 'all_found_psfs_test_select_psfs_G2.txt not generated.'


test_data = PARAMETRIZED_DATA['test_select_psfs_without_file']
select_psfs_without_file_parameters = [(CONVERTED_NIRCAM_IM_GA, False, False, 21, test_data['non-ga']),
                                       (CONVERTED_NIRCAM_IM_GA, True, False, 18, test_data['ga']),
                                       (CONVERTED_NIRCAM_IM_MIMF, False, True, 1, test_data['mimf'])
                                       ]
@pytest.mark.parametrize('in_data, ga, no_smooth, n_psfs, correct_all_found_psfs_txt', select_psfs_without_file_parameters)
def test_select_psfs_without_file(test_directory, in_data, ga, no_smooth, n_psfs, correct_all_found_psfs_txt):
    guiding_selections_path, all_found_psfs_path = select_psfs(
        in_data, ROOT, 2, global_alignment=ga, no_smoothing=no_smooth, testing=True,
        out_dir=test_directory
    )

    # Ensure the correct files were generated
    assert os.path.exists(guiding_selections_path), 'guiding_selections_test_select_psfs_G2.txt not generated.'
    assert os.path.exists(all_found_psfs_path), 'all_found_psfs_test_select_psfs_G2.txt not generated.'
    if no_smooth is True:
        main_path = guiding_selections_path.split('/guiding_selections')[0]
        no_smooth_path = os.path.join(main_path, 'psf_center_{}_G2.txt'.format(ROOT))
        assert os.path.exists(no_smooth_path), 'psf_center_test_select_psfs_G2.txt not generated.'

    # Test that the contents of all_found_psfs_test_select_psfs_G2.txt is correct
    # Right number of PSFs found?
    all_segment_locations = asc.read(all_found_psfs_path)
    assert len(all_segment_locations) == n_psfs

    # Right file contents
    with open(all_found_psfs_path) as f:
        all_found_psfs_contents = f.read()
    assert all_found_psfs_contents == correct_all_found_psfs_txt

    # Check PSF location in guiding selections doesn't match psf center
    # guiding selections should have found a knot in the PSF not at the center
    if no_smooth is True:
        with open(guiding_selections_path) as f:
            guiding_selections_contents = f.read()
        with open(no_smooth_path) as f:
            no_smooth_contents = f.read()
        assert guiding_selections_contents != no_smooth_contents
