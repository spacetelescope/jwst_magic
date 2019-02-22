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
    os.mkdir(test_dir)  # creates directory with default mode=511

    yield test_dir
    print("teardown test directory")
    if os.path.isdir(test_dir):
        shutil.rmtree(test_dir)



shift_to_id_attitude_parameters = [(False, 1, [(1024.0, 1024.0), (997.0, 1178.0), (1048.0, 1179.0), (1072.0, 1247.0), (1023.0, 1275.0)], [(1024.0, 1024.0), (1021.0, 1162.0), (997.0, 1178.0), (1048.0, 1179.0), (1072.0, 1190.0), (996.0, 1205.0), (1047.0, 1205.0), (1072.0, 1219.0), (997.0, 1233.0), (1048.0, 1233.0), (1022.0, 1247.0), (1072.0, 1247.0), (998.0, 1258.0), (1045.0, 1259.0), (1023.0, 1275.0)]),
                                   (True, 1, [(986.0, 1688.0), (959.0, 1842.0), (1010.0, 1843.0), (1034.0, 1911.0), (985.0, 1939.0)], [(986.0, 1688.0), (983.0, 1826.0), (959.0, 1842.0), (1010.0, 1843.0), (1034.0, 1854.0), (958.0, 1869.0), (1009.0, 1869.0), (1034.0, 1883.0), (959.0, 1897.0), (1010.0, 1897.0), (984.0, 1911.0), (1034.0, 1911.0), (960.0, 1922.0), (1007.0, 1923.0), (985.0, 1939.0)]),
                                   (True, 2, [(1003.0, 1697.0), (976.0, 1851.0), (1027.0, 1852.0), (1051.0, 1920.0), (1002.0, 1948.0)], [(1003.0, 1697.0), (1000.0, 1835.0), (976.0, 1851.0), (1027.0, 1852.0), (1051.0, 1863.0), (975.0, 1878.0), (1026.0, 1878.0), (1051.0, 1892.0), (976.0, 1906.0), (1027.0, 1906.0), (1001.0, 1920.0), (1051.0, 1920.0), (977.0, 1931.0), (1024.0, 1932.0), (1002.0, 1948.0)])]
@pytest.mark.parametrize('crowded_field, guider, guiding_selections_coords, all_found_psfs_coords', shift_to_id_attitude_parameters)
def test_shift_to_id_attitude(test_directory, crowded_field, guider, guiding_selections_coords, all_found_psfs_coords):
    with fits.open(FGS_CMIMF_IM) as hdulist:
        fgs_data = hdulist[1].data

    # Run the code
    BFS = BuildFGSSteps(fgs_data, guider, ROOT, 'ID', guiding_selections_file=SELECTED_SEGS,
                        out_dir=__location__, catalog=SEGMENT_INFILE, crowded_field=crowded_field)

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
    coords = [(x, y) for (x, y) in guiding_selections_cat['x', 'y']]
    assert coords == guiding_selections_coords

    # Make sure the shifted all_found_psfs*.txt is correct
    all_found_psfs_cat = asc.read(all_found_psfs_file)
    coords = [(x, y) for (x, y) in all_found_psfs_cat['x', 'y']]
    assert coords == all_found_psfs_coords

    # Make sure the location of the PSFs in the image matches the all_found_psfs*.txt
    with fits.open(FGS_img) as hdulist:
        image = hdulist[0].data
        sources = find_peaks(image, np.max(image) * .8, box_size=5)
        img_coords = [(x, y) for (x, y) in sources['x_peak', 'y_peak']]
        assert set(img_coords) == set(coords)



