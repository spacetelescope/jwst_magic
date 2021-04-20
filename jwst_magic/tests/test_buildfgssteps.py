"""Collection of unit tests to verify the correct function of the
buildfgssteps module.

Authors
-------
    - Lauren Chambers
    - Shannon Osborne

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

from jwst_magic.tests.utils import parametrized_data
from jwst_magic.fsw_file_writer import write_files
from jwst_magic.fsw_file_writer.buildfgssteps import BuildFGSSteps, shift_to_id_attitude
from jwst_magic.fsw_file_writer.buildfgssteps import OSS_TRIGGER, COUNTRATE_CONVERSION, DIM_STAR_THRESHOLD_FACTOR, \
    BRIGHT_STAR_THRESHOLD_ADDEND
from jwst_magic.fsw_file_writer.rewrite_prc import rewrite_prc
from jwst_magic.utils import utils

__location__ = os.path.realpath(os.path.join(os.getcwd(), os.path.dirname(__file__)))
FGS_CMIMF_IM = os.path.join(__location__, 'data', 'fgs_data_2_cmimf.fits')
NIRCAM_IM = os.path.join(__location__, 'data', 'nircam_data_1_ga.fits')
CONVERTED_NIRCAM_IM_MIMF = os.path.join(__location__, 'data', 'converted_nircam_data_1_mimf.fits')
ROOT = "test_buildfgssteps"
SEGMENT_INFILE_CMIMF = os.path.join(__location__, 'data', 'unshifted_all_found_psfs_{}_G1.txt'.format(ROOT))
SELECTED_SEGS_CMIMF_OLD = os.path.join(__location__, 'data', '{}_regfile.txt'.format(ROOT))
SELECTED_SEGS_CMIMF = os.path.join(__location__, 'data', 'unshifted_guiding_selections_{}_G1_config1.txt'.format(ROOT))
SELECTED_SEGS_MIMF = os.path.join(__location__, 'data', 'guiding_selections_nircam_data_1_mimf.txt')
PSF_CENTER_MIMF = os.path.join(__location__, 'data', 'psf_center_test_buildfgssteps_G1.txt')
ALL_PSFS_MIMF = os.path.join(__location__, 'data', 'all_found_psfs_buildfgssteps_G1.txt')
CENTER_POINTING_1 = os.path.join(__location__, 'data', 'center_pointing_{}_G1.txt'.format(ROOT))
CENTER_POINTING_2 = os.path.join(__location__, 'data', 'center_pointing_{}_2_G1.txt'.format(ROOT))

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
    (SELECTED_SEGS_CMIMF_OLD, False, 1, test_data['guiding_selections_coords'][0], test_data['all_found_psfs_coords'][0], CENTER_POINTING_1),
    (SELECTED_SEGS_CMIMF, False, 1, test_data['guiding_selections_coords'][0], test_data['all_found_psfs_coords'][0], CENTER_POINTING_1),
    (SELECTED_SEGS_CMIMF, True, 1, test_data['guiding_selections_coords'][1], test_data['all_found_psfs_coords'][1], CENTER_POINTING_2),
    (SELECTED_SEGS_CMIMF, True, 2, test_data['guiding_selections_coords'][2], test_data['all_found_psfs_coords'][2], CENTER_POINTING_2)
]
@pytest.mark.parametrize('guiding_selections, crowded_field, guider, guiding_selections_coords, all_found_psfs_coords,'
                         'center_pointing', shift_to_id_attitude_parameters)
def test_shift_to_id_attitude(open_image, test_directory, guiding_selections, crowded_field, guider,
                              guiding_selections_coords, all_found_psfs_coords, center_pointing):

    #Run the prep code that's in run_magic.py
    if 'guiding_config' in guiding_selections:
        out_dir_fsw = os.path.join(TEST_DIRECTORY, 'guiding_config_{}'.format(
            guiding_selections.split('guiding_config_')[1].split('/')[0]))
    else:
        out_dir_fsw = TEST_DIRECTORY

    # Run main function to test
    fgs_im, guiding_selections_file, psf_center_file = shift_to_id_attitude(
        open_image, ROOT, guider, out_dir_fsw, guiding_selections_file=guiding_selections,
        all_found_psfs_file=SEGMENT_INFILE_CMIMF, center_pointing_file=center_pointing,
        psf_center_file=None, crowded_field=crowded_field, logger_passed=True)

    # Define filenames
    file_root = guiding_selections.split('/')[-1]
    if 'regfile' in file_root:
        file_root = file_root.split('.txt')[0].split('_regfile')[0]
    elif 'unshifted_guiding_selections' in file_root:
        file_root = file_root.split('guiding_selections_')[-1].split('.txt')[0]
    if 'G{}'.format(guider) not in file_root:
        file_root += '_G{}'.format(guider)

    guiding_selections_file = os.path.join(out_dir_fsw, 'shifted_guiding_selections_{}.txt'.format(file_root))
    all_found_psfs_file = os.path.join(out_dir_fsw, 'shifted_all_found_psfs_{}.txt'.format(file_root))
    FGS_img = os.path.join(out_dir_fsw,  'FGS_imgs', 'shifted_' + file_root + '.fits')
    center_pointing_file = os.path.join(out_dir_fsw, 'shifted_center_pointing_{}.txt'.format(file_root))

    # Check that the right files were put in the right place
    assert os.path.exists(FGS_img)
    assert os.path.exists(guiding_selections_file)
    assert os.path.exists(all_found_psfs_file)
    assert os.path.exists(center_pointing_file)

    # Make sure the shifted guiding_selections*.txt is correct
    guiding_selections_cat = asc.read(guiding_selections_file)
    coords = np.array([(x, y) for (x, y) in guiding_selections_cat['x', 'y']])
    assert np.array_equal(coords, guiding_selections_coords)

    # Make sure the shifted all_found_psfs*.txt is correct
    all_found_psfs_cat = asc.read(all_found_psfs_file)
    coords = np.array([(x, y) for (x, y) in all_found_psfs_cat['x', 'y']])
    assert np.array_equal(coords, all_found_psfs_coords)

    # Make sure the shifted center_pointing file is correct (if its a list, it's shifted, if it's not, its unchanged)
    center_pointing_cat = asc.read(center_pointing_file, format='commented_header', delimiter=',')
    label = center_pointing_cat.colnames[0]
    center_pointing_cat_original = asc.read(center_pointing, format='commented_header', delimiter=',')
    if isinstance(center_pointing_cat[label][0], str):
        assert center_pointing_cat[label][0] != center_pointing_cat_original[label][0]
    else:
        assert center_pointing_cat[label][0] == center_pointing_cat_original[label][0]

    # Make sure the location of the PSFs in the image matches the all_found_psfs*.txt
    with fits.open(FGS_img) as hdulist:
        image = hdulist[0].data
        sources = find_peaks(image, np.median(image+ (3 * np.std(image))), box_size=5)
        img_coords = [(x, y) for (x, y) in sources['x_peak', 'y_peak']]
        assert set(img_coords) == set([(int(x), int(y)) for (x, y) in all_found_psfs_cat['x', 'y']])


correct_count_rate_parameters = []
test_data = PARAMETRIZED_DATA['test_correct_count_rate']
for guider in [1, 2]:
    for step in ['CAL', 'ID', 'ACQ1', 'ACQ2', 'TRK', 'LOSTRK']:
        g = 'guider{}'.format(guider)
        if step in ['CAL', 'ID', 'LOSTRK']:
            correct_count_rate_parameters.append((guider, step,
                                                  test_data[g][step]))
        else:
            correct_count_rate_parameters.append(pytest.param(guider, step, test_data[g][step],
                                                 marks=pytest.mark.xfail(reason="bias issue: see JWSTFGS-213")))
@pytest.mark.parametrize('guider, step, correct_data_dict', correct_count_rate_parameters)
def test_correct_count_rate(open_image, test_directory, guider, step, correct_data_dict):
    """Check that image data is being generated with counts and count
    rates as expected. Test for all guider steps and both guiders.
    """

    # Input data
    assert np.isclose(np.min(open_image), 0.0), 'Incorrect input data min - changed your file?'
    assert np.isclose(np.max(open_image), 15693.34), 'Incorrect input data max - changed your file?'

    # Run the code
    fgs_im, guiding_selections_file, psf_center_file = shift_to_id_attitude(
        open_image, ROOT, guider, TEST_DIRECTORY, guiding_selections_file=SELECTED_SEGS_CMIMF_OLD,
        all_found_psfs_file=SEGMENT_INFILE_CMIMF, center_pointing_file=CENTER_POINTING_1,
        psf_center_file=None, crowded_field=False, logger_passed=True)
    BFS = BuildFGSSteps(fgs_im, guider, ROOT, step, guiding_selections_file=guiding_selections_file,
                        out_dir=TEST_DIRECTORY, shift_id_attitude=True)

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
        assertion_range = 250  # Counts;

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
        if step == 'ID':
            assert abs(correct_data_dict['strips'][0] - np.min(BFS.strips)) < assertion_range, \
                'ID strips counts out of expected range.'
            assert abs(correct_data_dict['strips'][1] - np.max(BFS.strips)) < assertion_range, \
                'ID strips counts out of expected range.'

        # Final step product
        assert abs(correct_data_dict[step][0] - np.min(BFS.image)) < assertion_range, \
            '{} counts out of expected range.'.format(step)
        assert abs(correct_data_dict[step][1] - np.max(BFS.image)) < assertion_range, \
            '{} counts out of expected range.'.format(step)

    # Assert exact count rates
    assert (BFS.countrate == correct_data_dict['countrates']).all(), \
        'Incorrect {} count rate.'.format(step)


def test_psf_center_file(test_directory):
    """Test that when psf_center_file is set, the array position for TRK
    is pulled from the psf_center file rather than the guiding selections file
    which is used for all other steps.
    """
    image = fits.getdata(CONVERTED_NIRCAM_IM_MIMF, 0)
    guider = 1
    shift_id_attitude = False

    # Run the code
    fileobj_id = BuildFGSSteps(
        image, guider, ROOT, step='ID', guiding_selections_file=SELECTED_SEGS_MIMF,
        out_dir=TEST_DIRECTORY, psf_center_file=PSF_CENTER_MIMF, shift_id_attitude=shift_id_attitude
    )
    fileobj_acq1 = BuildFGSSteps(
        image, guider, ROOT, step='ACQ1', guiding_selections_file=SELECTED_SEGS_MIMF,
        out_dir=TEST_DIRECTORY, psf_center_file=PSF_CENTER_MIMF, shift_id_attitude=shift_id_attitude
    )

    fileobj_trk = BuildFGSSteps(
        image, guider, ROOT, step='TRK', guiding_selections_file=SELECTED_SEGS_MIMF,
        out_dir=TEST_DIRECTORY, psf_center_file=PSF_CENTER_MIMF, shift_id_attitude=shift_id_attitude
    )

    assert (fileobj_id.xarr, fileobj_id.yarr) == (fileobj_acq1.xarr, fileobj_acq1.yarr)
    assert (fileobj_acq1.xarr, fileobj_acq1.yarr) != (fileobj_trk.xarr, fileobj_trk.yarr)
    assert fileobj_acq1.countrate == fileobj_trk.countrate


oss_defaults_parameters = [100000, 1000000]
@pytest.mark.parametrize('catalog_countrate', oss_defaults_parameters)
def test_oss_defaults(test_directory, catalog_countrate):
    """Test that when use_oss_defaults is set to True, the attributes
    of the file object (which will be used when writing out all FSW files)
    are set appropriatly.
    """
    image = fits.getdata(CONVERTED_NIRCAM_IM_MIMF, 0)
    guider = 1
    use_oss_defaults = True

    fileobj = BuildFGSSteps(
        image, guider, ROOT, step='ID', guiding_selections_file=SELECTED_SEGS_MIMF,
        out_dir=TEST_DIRECTORY, psf_center_file=PSF_CENTER_MIMF, shift_id_attitude=False,
        use_oss_defaults=use_oss_defaults, catalog_countrate=catalog_countrate)

    # Compare the countrate and threshold to what's expected
    assert fileobj.countrate == catalog_countrate * COUNTRATE_CONVERSION

    if catalog_countrate < OSS_TRIGGER:
        assert fileobj.threshold == catalog_countrate * COUNTRATE_CONVERSION * DIM_STAR_THRESHOLD_FACTOR
    else:
        assert fileobj.threshold == (catalog_countrate * COUNTRATE_CONVERSION) - BRIGHT_STAR_THRESHOLD_ADDEND


def test_rewrite_prc(open_image, test_directory):
    """Compare the results from reqrite_prc and buildfgsteps -
    check shifted guiding selections and ID prc file"""

    # Delete path if it exists
    if os.path.isdir(TEST_DIRECTORY):
        shutil.rmtree(TEST_DIRECTORY)

    # Define inputs
    inds_list = [[0, 2, 3, 11, 14]] # inds that match the guiding_selections_file
    guiding_selections_file = os.path.join(TEST_DIRECTORY, 'guiding_config_1', SELECTED_SEGS_CMIMF.split('/')[-1])
    all_found_psfs_file = os.path.join(TEST_DIRECTORY, SEGMENT_INFILE_CMIMF.split('/')[-1])
    center_of_pointing = 0
    center_of_pointing_file = CENTER_POINTING_1
    guider = 1
    thresh_factor = 0.5
    shifted = True
    crowded_field = False
    step = 'ACQ1'

    # Copy file for testing
    os.makedirs(os.path.join(TEST_DIRECTORY, 'FGS_imgs'))
    os.makedirs(os.path.join(TEST_DIRECTORY, 'guiding_config_1'))
    shutil.copyfile(FGS_CMIMF_IM, os.path.join(TEST_DIRECTORY, 'FGS_imgs', 'unshifted_test_buildfgssteps_G1.fits'))
    shutil.copyfile(SEGMENT_INFILE_CMIMF, all_found_psfs_file)
    shutil.copyfile(SELECTED_SEGS_CMIMF, guiding_selections_file)

    # Run buildfgssteps
    out_fsw = os.path.join(TEST_DIRECTORY, 'guiding_config_1')
    fgs_im, guiding_selections_file, psf_center_file = shift_to_id_attitude(
        open_image, ROOT, guider, out_fsw, guiding_selections_file=guiding_selections_file,
        all_found_psfs_file=all_found_psfs_file, center_pointing_file=center_of_pointing_file,
        psf_center_file=None, crowded_field=crowded_field, logger_passed=True)
    BFS = BuildFGSSteps(fgs_im, guider, ROOT, step, guiding_selections_file=guiding_selections_file,
                        out_dir=out_fsw, thresh_factor=thresh_factor, shift_id_attitude=shifted)
    write_files.write_prc(BFS)

    # Check for output files
    shifted_guiding_selections = os.path.join(out_fsw,  'shifted_guiding_selections_{}_config1.txt'.format(ROOT+'_G1'))
    shifted_acq_prc = os.path.join(out_fsw, 'dhas_shifted', ROOT+'_G1_ACQ.prc')
    assert os.path.exists(shifted_guiding_selections)
    assert os.path.exists(shifted_acq_prc)
    buildsteps_selections = asc.read(shifted_guiding_selections)
    with open(shifted_acq_prc, 'r') as file:
        buildsteps_prc = file.read()

    # Run rewrite_prc
    rewrite_prc(inds_list, center_of_pointing, guider, ROOT, __location__, thresh_factor, shifted, crowded_field)

    # Check for output files
    shifted_guiding_selections2 = os.path.join(TEST_DIRECTORY, 'guiding_config_2',
                                               'shifted_guiding_selections_{}_config2.txt'.format(ROOT+'_G1'))
    shifted_acq_prc2 = os.path.join(TEST_DIRECTORY, 'guiding_config_2', 'dhas_shifted', ROOT+'_G1_ACQ.prc')
    assert os.path.exists(shifted_guiding_selections2)
    assert os.path.exists(shifted_acq_prc2)
    rewrite_prc_selections = asc.read(shifted_guiding_selections2)
    with open(shifted_acq_prc2, 'r') as file:
        rewrite_prc_prc = file.read()

    # Confirm outputs match
    assert str(rewrite_prc_selections) == str(buildsteps_selections)
    assert rewrite_prc_prc == buildsteps_prc


prc_list = [('ID', 'ID', False, 237576.0000, None),
            ('ACQ1', 'ACQ', False, 237576.0000, None),
            ('ACQ1', 'ACQ', True, None, 100000)]
@pytest.mark.parametrize('step, step_name, use_oss_defaults, guide_star_countrate, catalog_countrate', prc_list)
def test_prc_thresholds(test_directory, step, step_name, use_oss_defaults, guide_star_countrate, catalog_countrate):
    """Check the right thresholds make it into the ACQ prc files"""
    # Delete path if it exists
    if os.path.isdir(TEST_DIRECTORY):
        shutil.rmtree(TEST_DIRECTORY)

    # Define basic inputs
    image = fits.getdata(CONVERTED_NIRCAM_IM_MIMF)
    guider = 1
    shifted = False
    guiding_selections_file = SELECTED_SEGS_MIMF

    # Set the threshold factor variable
    thresh_factor = 0.5

    # Run buildfgssteps
    out_fsw = os.path.join(TEST_DIRECTORY, 'guiding_config_1')
    BFS_factor = BuildFGSSteps(image, guider, ROOT, step, guiding_selections_file=guiding_selections_file,
                               out_dir=out_fsw, thresh_factor=thresh_factor, shift_id_attitude=shifted,
                               use_oss_defaults=use_oss_defaults, catalog_countrate=catalog_countrate)
    write_files.write_prc(BFS_factor)

    # Check threshold in the output files
    if step_name == 'ID':
        directory = 'ground_system'
    elif step_name == 'ACQ':
        directory = 'dhas'

    thresh_factor_prc = os.path.join(out_fsw, directory, ROOT+f'_G1_{step_name}.prc')
    assert os.path.exists(thresh_factor_prc)
    with open(thresh_factor_prc, 'r') as file:
        buildsteps_factor_prc = file.read()
    if use_oss_defaults:
        assert str(int(catalog_countrate * COUNTRATE_CONVERSION * DIM_STAR_THRESHOLD_FACTOR)) in buildsteps_factor_prc
    else:
        assert str(int(guide_star_countrate * thresh_factor)) in buildsteps_factor_prc


star_list = [(False, 237576.0000, None),  # gs countrate is from the above guiding selections file
             (True, None, 100000)]
@pytest.mark.parametrize('use_oss_defaults, guide_star_countrate, catalog_countrate', star_list)
def test_star_thresholds(test_directory, use_oss_defaults, guide_star_countrate, catalog_countrate):
    """Check the right thresholds make it into the ID star files"""
    step = 'ID'

    # Delete path if it exists
    if os.path.isdir(TEST_DIRECTORY):
        shutil.rmtree(TEST_DIRECTORY)

    # Define basic inputs
    image = fits.getdata(CONVERTED_NIRCAM_IM_MIMF)
    guider = 1
    shifted = False
    guiding_selections_file = SELECTED_SEGS_MIMF

    # Set the threshold factor variable
    thresh_factor = 0.5

    # Run buildfgssteps
    out_fsw = os.path.join(TEST_DIRECTORY, 'guiding_config_1')
    BFS_factor = BuildFGSSteps(image, guider, ROOT, step, guiding_selections_file=guiding_selections_file,
                               out_dir=out_fsw, thresh_factor=thresh_factor, shift_id_attitude=shifted,
                               use_oss_defaults=use_oss_defaults, catalog_countrate=catalog_countrate)
    write_files.write_star(BFS_factor)

    # Check threshold in the output files (0th index for GS, 4th index for threshold)
    thresh_factor_star = os.path.join(out_fsw, 'dhas', ROOT+f'_G1_{step}.star')
    assert os.path.exists(thresh_factor_star)
    star_data = asc.read(thresh_factor_star, data_start=1)
    if use_oss_defaults:
        assert np.isclose(int(catalog_countrate * COUNTRATE_CONVERSION * DIM_STAR_THRESHOLD_FACTOR), star_data[0][4])
    else:
        assert np.isclose(int(guide_star_countrate * thresh_factor), star_data[0][4])
