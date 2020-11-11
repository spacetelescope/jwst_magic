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
from collections import OrderedDict
import io
import os
import shutil
import yaml

from astropy.io import ascii as asc
import numpy as np
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
        CONVERTED_NIRCAM_IM_GA, ROOT, 2, guiding_selections_file_list=guiding_selection_file,
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

test_data = PARAMETRIZED_DATA['test_yaml_file']
test_yaml_file_parameters = [(test_data['cols_order'], test_data['cols_list'], test_data['all_cols'])]
@pytest.mark.parametrize('cols_order, cols_list, all_cols', test_yaml_file_parameters)
def test_yaml_file(cols_order, cols_list, all_cols):
    guider = 1
    testing_num = 11
    out_dir = os.path.join(__location__, 'out', ROOT)
    new_config_numbers = np.arange(1, testing_num+1)

    # Write guiding selections, all found_psfs, and all_selections_yaml files to fake a previous observation
    guiding_selections_path_list = []
    for (i, cols) in zip(new_config_numbers, cols_list):
        guiding_selections_path = os.path.join(out_dir, 'guiding_config_{}'.format(i),
                                            'unshifted_guiding_selections_{}_G{}_config{}.txt'.format(ROOT, guider, i))
        utils.write_cols_to_file(guiding_selections_path, labels=['y', 'x', 'countrate'], cols=cols)
        guiding_selections_path_list.append(guiding_selections_path)

    all_found_psfs_path = os.path.join(out_dir, 'unshifted_all_found_psfs_{}_G{}.txt'.format(ROOT, guider))
    utils.write_cols_to_file(all_found_psfs_path, labels=['label', 'y', 'x', 'countrate'], cols=all_cols)

    old_yaml_path = os.path.join(out_dir, 'all_guiding_selections.yaml')
    utils.setup_yaml()
    final_data = OrderedDict()
    for i, file in enumerate(guiding_selections_path_list):
        config_data = cols_order[i]
        final_data['guiding_config_{}'.format(i+1)] = config_data
    with io.open(old_yaml_path, 'w', encoding="utf-8") as f:
        yaml.dump(final_data, f, default_flow_style=False, allow_unicode=True)

    # Simulate loading these 11 configs via a file
    guiding_selections_path_list, all_found_psfs_path, _ = select_psfs(
        CONVERTED_NIRCAM_IM_GA, ROOT, guider, guiding_selections_file_list=guiding_selections_path_list,
        out_dir=__location__
    )

    # Check file exists
    yaml_file = os.path.join(out_dir, 'all_guiding_selections.yaml')
    assert os.path.exists(yaml_file), 'all_guiding_selections.yaml not generated.'

    # Check contents of file -> 11 Configs in correct number order
    with open(yaml_file, 'r') as stream:
        data_loaded = yaml.safe_load(stream)
    assert list(data_loaded.keys()) ==  ['guiding_config_{}'.format(num+1) for num in range(testing_num)]

    # Simulate re-loading the same config again + 1 more new config
    guiding_selections_path_list, all_found_psfs_path, _ = select_psfs(
        CONVERTED_NIRCAM_IM_GA, ROOT, guider, guiding_selections_file_list=[guiding_selections_path_list[0], SELECTED_SEGS],
        out_dir=__location__
    )

    # Re-load contents of yaml file
    with open(yaml_file, 'r') as stream:
        data_loaded = yaml.safe_load(stream)

    # Another config was added, but check they are still in number order
    assert len(list(data_loaded.keys())) == testing_num+1
    assert list(data_loaded.keys()) == ['guiding_config_{}'.format(num+1) for num in range(testing_num+1)]

    # Check the contents of the yaml file matches what it should
    for i, config in enumerate(cols_order):
        assert data_loaded['guiding_config_{}'.format(i+1)] == config

    # The new selections file added didn't have a yaml file, so an empty list was appended to this yaml
    assert len(data_loaded['guiding_config_{}'.format(testing_num+1)]) == 0