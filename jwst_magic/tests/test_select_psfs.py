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

select_psfs_with_file_parameters = [([SELECTED_SEGS], [1]), ([SELECTED_SEGS, SELECTED_SEGS], [2,3])]
@pytest.mark.parametrize('guiding_selection_file, expected_configs', select_psfs_with_file_parameters)
def test_select_psfs_with_file(test_directory, guiding_selection_file, expected_configs):
    # Run code
    guiding_selections_path_list, all_found_psfs_path, _ = select_psfs(
        CONVERTED_NIRCAM_IM_GA, ROOT, 2, guiding_selections_file=guiding_selection_file,
        out_dir=__location__
    )

    # Check for multiple commands
    assert len(guiding_selection_file) == len(guiding_selections_path_list)

    # Ensure the correct files were generated based on output path information nad have the right file names
    # Check guiding selections file(s)
    for i in range(len(guiding_selection_file)):
        assert os.path.exists(guiding_selections_path_list[i]), \
               'unshifted_guiding_selections_test_select_psfs_G2_config{}.txt not generated.'.format(expected_configs[i])
        assert sorted(guiding_selections_path_list)[i] == \
               os.path.join(TEST_DIRECTORY, 'guiding_config_{}'.format(expected_configs[i]),
              'unshifted_guiding_selections_test_select_psfs_G2_config{}.txt'.format(expected_configs[i]))

    # Check all found PSF file
    assert os.path.exists(all_found_psfs_path), \
        'unshifted_all_found_psfs_test_select_psfs_G2.txt not generated.'

    assert all_found_psfs_path == os.path.join(TEST_DIRECTORY, 'unshifted_all_found_psfs_test_select_psfs_G2.txt')



test_data = PARAMETRIZED_DATA['test_select_psfs_without_file']
select_psfs_without_file_parameters = [(CONVERTED_NIRCAM_IM_GA, 'default', 21, test_data['non-ga']),
                                       (CONVERTED_NIRCAM_IM_GA, 'high', 18, test_data['ga']),
                                       (CONVERTED_NIRCAM_IM_MIMF, 'low', 1, test_data['mimf'])
                                       ]
@pytest.mark.parametrize('in_data, smooth, n_psfs, correct_all_found_psfs_txt', select_psfs_without_file_parameters)
def test_select_psfs_without_file(test_directory, in_data, smooth, n_psfs, correct_all_found_psfs_txt):
    guiding_selections_path_list, all_found_psfs_path, psf_center_path = select_psfs(
        in_data, ROOT, 2, smoothing=smooth, testing=True,
        out_dir=__location__
    )

    # Ensure the correct files were generated (1 selection by default in code for test cases)
    assert len(guiding_selections_path_list) == 1
    assert os.path.exists(guiding_selections_path_list[0]), 'unshifted_guiding_selections_test_select_psfs_G2.txt not generated.'
    assert os.path.exists(all_found_psfs_path), 'unshifted_all_found_psfs_test_select_psfs_G2.txt not generated.'
    if smooth is 'low':
        assert os.path.exists(psf_center_path), 'unshifted_psf_center_test_select_psfs_G2.txt not generated.'

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
    if smooth is 'low':
        with open(guiding_selections_path_list[0]) as f:
            guiding_selections_contents = f.read()
        with open(psf_center_path) as f:
            no_smooth_contents = f.read()
        assert guiding_selections_contents != no_smooth_contents
