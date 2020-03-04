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

from astropy.io import ascii as asc
from astropy.io import fits
import numpy as np
from photutils import find_peaks
import pytest

from .utils import parametrized_data
from ..fsw_file_writer import write_files
from ..fsw_file_writer.buildfgssteps import BuildFGSSteps
from ..utils import utils

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

PARAMETRIZED_DATA = parametrized_data()['test_buildfgssteps']


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


@pytest.fixture(scope="module")
def open_image(image=FGS_CMIMF_IM):
    """Open a FITS image into a numpy array.

    Parameters
    ----------
    image : str
        Path to image to open

    Yields
    -------
    fgs_data : 2-D numpy array
        Array of image data
    """
    with fits.open(image) as hdulist:
        fgs_data = hdulist[1].data

    yield fgs_data


test_data = PARAMETRIZED_DATA['test_shift_to_id_attitude']
shift_to_id_attitude_parameters = [
    (False, 1, test_data['guiding_selections_coords'][0], test_data['all_found_psfs_coords'][0]),
    (True, 1, test_data['guiding_selections_coords'][1], test_data['all_found_psfs_coords'][1]),
    (True, 2, test_data['guiding_selections_coords'][2], test_data['all_found_psfs_coords'][2])
]
@pytest.mark.parametrize('crowded_field, guider, guiding_selections_coords, all_found_psfs_coords',
                         shift_to_id_attitude_parameters)
def test_shift_to_id_attitude(open_image, test_directory, crowded_field, guider,
                              guiding_selections_coords, all_found_psfs_coords):
    # Run the code
    BFS = BuildFGSSteps(open_image, guider, ROOT, 'ID', guiding_selections_file=SELECTED_SEGS,
                        out_dir=__location__, catalog=SEGMENT_INFILE, crowded_field=crowded_field)
    write_files.write_all(BFS)

    # Define filenames
    file_root = '{}_G{}'.format(ROOT, guider)
    guiding_selections_file = os.path.join(test_directory, 'shifted',
                                           'guiding_selections_{}.txt'.format(file_root))
    all_found_psfs_file = os.path.join(test_directory, 'shifted',
                                       'all_found_psfs_{}.txt'.format(file_root))
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


correct_count_rate_parameters = []
test_data = PARAMETRIZED_DATA['test_correct_count_rate']
for guider in [1, 2]:
    for step in ['CAL', 'ID', 'ACQ1', 'ACQ2', 'TRK', 'LOSTRK']:
        g = 'guider{}'.format(guider)
        correct_count_rate_parameters.append((guider, step,
                                              test_data[g][step]))
@pytest.mark.parametrize('guider, step, correct_data_dict', correct_count_rate_parameters)
def test_correct_count_rate(open_image, guider, step, correct_data_dict):
    """Check that image data is being generated with counts and count
    rates as expected. Test for all guider steps and both guiders.
    """

    # Input data
    assert np.isclose(np.min(open_image), 0.0), 'Incorrect input data min - changed your file?'
    assert np.isclose(np.max(open_image), 15693.34), 'Incorrect input data max - changed your file?'

    # Run the code
    BFS = BuildFGSSteps(open_image, guider, ROOT, step, guiding_selections_file=SELECTED_SEGS,
                        out_dir=__location__, catalog=SEGMENT_INFILE)

    # Assert ~exactly for time-normalized data (before detector effects are added)
    assert np.isclose(correct_data_dict['time_normed_im'][0], np.min(BFS.time_normed_im)), \
        'Incorrect {} time-normalized counts.'.format(step)
    assert np.isclose(correct_data_dict['time_normed_im'][1], np.max(BFS.time_normed_im)), \
        'Incorrect {} time-normalized counts.'.format(step)

    # Assert ~exactly for LOSTRK (where there is no bias and thus no noise)
    if step == 'LOSTRK':
        assert np.isclose(correct_data_dict[step][0], np.min(BFS.image)), \
            'Incorrect {} counts.'.format(step)
        assert np.isclose(correct_data_dict[step][1], np.max(BFS.image)), \
            'Incorrect {} counts.'.format(step)

    # For other steps/data types, assert within a range, to account for
    # randomness introduced by detector effects
    else:
        assertion_range = 20  # Counts; Note full range is two times this

        # Bias
        assert abs(correct_data_dict['bias'][0] - np.min(BFS.bias)) < assertion_range, \
            '{} bias image counts out of expected range.'.format(step)
        assert abs(correct_data_dict['bias'][1] - np.max(BFS.bias)) < assertion_range, \
            '{} bias image counts out of expected range.'.format(step)

        # CDS
        if 'TRK' not in step:
            assert abs(correct_data_dict['cds'][0] - np.min(BFS.cds)) < assertion_range, \
                '{} CDS counts out of expected range.'.format(step)
            assert abs(correct_data_dict['cds'][1] - np.max(BFS.cds)) < assertion_range, \
                '{} CDS counts out of expected range.'.format(step)

        # Strips
        if step is 'ID':
            assert abs(correct_data_dict['strips'][0] - np.min(BFS.strips)) < assertion_range, \
                'ID strips counts out of expected range.'
            assert abs(correct_data_dict['strips'][1] - np.max(BFS.strips)) < assertion_range, \
                'ID strips counts out of expected range.'

        # Final step product
        assert abs(correct_data_dict[step][0] - np.min(BFS.image)) < assertion_range, \
            '{} counts out of expected range.'.format(step)
        assert abs(correct_data_dict[step][1] - np.max(BFS.image)) < assertion_range, \
            '{} counts out of expected range.'.format(step)

    # Assert exact ount rates
    assert (BFS.countrate == correct_data_dict['countrates']).all(), \
        'Incorrect {} count rate.'.format(step)
