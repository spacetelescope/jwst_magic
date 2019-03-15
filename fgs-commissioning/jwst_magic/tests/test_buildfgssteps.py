"""Collection of unit tests to verify the correct function of the
buildfgssteps module.

Authors
-------
    - Lauren Chambers

Use
---
    ::
        pytest test_buildfgssteps.py
"""
import os
import shutil
import yaml

from astropy.io import ascii as asc
from astropy.io import fits
import numpy as np
from photutils import find_peaks
import pytest

from jwst_magic.fsw_file_writer import write_files
from jwst_magic.fsw_file_writer.buildfgssteps import BuildFGSSteps

__location__ = os.path.realpath(os.path.join(os.getcwd(), os.path.dirname(__file__)))
FGS_CMIMF_IM = os.path.join(__location__, 'data', 'fgs_data_2_cmimf.fits')
NIRCAM_IM = os.path.join(__location__, 'data', 'nircam_data_1_ga.fits')
ROOT = "test_buildfgssteps"
SEGMENT_INFILE = os.path.join(__location__, 'data', '{}_ALLpsfs.txt'.format(ROOT))
SELECTED_SEGS = os.path.join(__location__, 'data', '{}_regfile.txt'.format(ROOT))
# PROGRAM_ID = 1141
# OBSERVATION_NUM = 7
VISIT_NUM = 1

TEST_DIRECTORY = os.path.join(__location__, 'out', ROOT)

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
    os.makedirs(test_dir)  # creates directory with default mode=511

    yield test_dir
    print("teardown test directory")
    if os.path.isdir(test_dir):
        shutil.rmtree(test_dir)


def parametrized_data():
    """Load parametrized data from file.

    Returns
    -------
    test_data : dict
        Dictionary containing parametrized test data
    """
    parametrized_data_file = os.path.join(__location__, 'data', 'parametrized_test_data.yml')
    with open(parametrized_data_file) as f:
        test_data = yaml.load(f.read())

    return test_data['test_buildfgssteps']


test_data = parametrized_data()['test_shift_to_id_attitude']
shift_to_id_attitude_parameters = [(False, 1, test_data['guiding_selections_coords'][0], test_data['all_found_psfs_coords'][0]),
                                   (True, 1, test_data['guiding_selections_coords'][1], test_data['all_found_psfs_coords'][1]),
                                   (True, 2, test_data['guiding_selections_coords'][2], test_data['all_found_psfs_coords'][2])]
@pytest.mark.parametrize('crowded_field, guider, guiding_selections_coords, all_found_psfs_coords', shift_to_id_attitude_parameters)
def test_shift_to_id_attitude(test_directory, crowded_field, guider, guiding_selections_coords, all_found_psfs_coords):
    with fits.open(FGS_CMIMF_IM) as hdulist:
        fgs_data = hdulist[1].data

    # Run the code
    BFS = BuildFGSSteps(fgs_data, guider, ROOT, 'ID', guiding_selections_file=SELECTED_SEGS,
                        out_dir=__location__, catalog=SEGMENT_INFILE, crowded_field=crowded_field)
    write_files.write_all(BFS)

    # Define filenames
    file_root = '{}_G{}'.format(ROOT, guider)
    guiding_selections_file = os.path.join(test_directory, 'shifted', 'guiding_selections_{}.txt'.format(file_root))
    all_found_psfs_file = os.path.join(test_directory, 'shifted', 'all_found_psfs_{}.txt'.format(file_root))
    FGS_img = os.path.join(test_directory, 'shifted', file_root + '.fits')

    # Check that the right files were put in the right place
    assert os.path.exists(FGS_img)
    assert os.path.exists(guiding_selections_file)
    assert os.path.exists(all_found_psfs_file)

    # Make sure the shifted guiding_selections*.txt is correct
    guiding_selections_cat = asc.read(guiding_selections_file)
    coords = np.array([(x, y) for (x, y) in guiding_selections_cat['x', 'y']])
    assert np.array_equal(coords, guiding_selections_coords)

    # Make sure the shifted all_found_psfs*.txt is correct
    all_found_psfs_cat = asc.read(all_found_psfs_file)
    coords = np.array([(x, y) for (x, y) in all_found_psfs_cat['x', 'y']])
    assert np.array_equal(coords, all_found_psfs_coords)

    # Make sure the location of the PSFs in the image matches the all_found_psfs*.txt
    with fits.open(FGS_img) as hdulist:
        image = hdulist[0].data
        sources = find_peaks(image, np.max(image) * .8, box_size=5)
        img_coords = [(x, y) for (x, y) in sources['x_peak', 'y_peak']]
        assert set(img_coords) == set([(x, y) for (x, y) in all_found_psfs_cat['x', 'y']])
