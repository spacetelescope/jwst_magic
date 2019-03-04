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
from astropy.io import fits
import numpy as np
from photutils import find_peaks
import pytest

from jwst_magic.star_selector.select_psfs import select_psfs

__location__ = os.path.realpath(os.path.join(os.getcwd(), os.path.dirname(__file__)))
FGS_CMIMF_IM = os.path.join(__location__, 'data', 'fgs_data_2_cmimf.fits')
NIRCAM_IM = os.path.join(__location__, 'data', 'nircam_data_1_ga.fits')
CONVERTED_NIRCAM_IM = os.path.join(__location__, 'data', 'converted_nircam_data_1_ga.fits')
ROOT = "test_select_psfs"
SEGMENT_INFILE = os.path.join(__location__, 'data', 'all_found_psfs_{}.txt'.format(ROOT))
SELECTED_SEGS = os.path.join(__location__, 'data', 'guiding_selections_{}.txt'.format(ROOT))
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

def test_select_psfs_with_file(test_directory):
    guiding_selections_path, all_found_psfs_path = select_psfs(
        CONVERTED_NIRCAM_IM, ROOT, 2, guiding_selections_file=SELECTED_SEGS,
        out_dir=test_directory
    )

    # Ensure the correct files were generated
    assert os.path.exists(guiding_selections_path), 'guiding_selections_test_select_psfs_G2.txt not generated.'
    assert os.path.exists(all_found_psfs_path), 'all_found_psfs_test_select_psfs_G2.txt not generated.'


all_found_psfs_without_file = '''# label y x countrate
F 707.0000 848.0000 8053.0000
O 716.0000 1258.0000 27255.0000
J 719.0000 1058.0000 15701.0000
F 722.0000 865.0000 5286.0000
F 739.0000 874.0000 4320.0000
H 861.0000 960.0000 47211.0000
Q 879.0000 1306.0000 18428.0000
C 880.0000 771.0000 15725.0000
L 882.0000 1166.0000 5537.0000
L 906.0000 1162.0000 5880.0000
A 1024.0000 692.0000 17285.0000
R 1032.0000 1412.0000 45795.0000
N 1042.0000 1245.0000 106899.0000
E 1048.0000 881.0000 29745.0000
P 1181.0000 1321.0000 30061.0000
K 1187.0000 1152.0000 32794.0000
B 1207.0000 801.0000 25434.0000
G 1223.0000 977.0000 82360.0000
I 1340.0000 1077.0000 41864.0000
M 1342.0000 1258.0000 31596.0000
D 1377.0000 866.0000 34275.0000
'''

all_found_psfs_without_file_ga ='''# label y x countrate
F 723.0000 860.0000 1183.0000
O 723.0000 1259.0000 17949.0000
J 724.0000 1060.0000 15276.0000
H 861.0000 961.0000 46459.0000
C 880.0000 771.0000 15725.0000
Q 880.0000 1306.0000 16855.0000
L 896.0000 1167.0000 3370.0000
A 1026.0000 697.0000 15731.0000
R 1033.0000 1414.0000 40876.0000
N 1042.0000 1246.0000 104626.0000
E 1049.0000 885.0000 18002.0000
P 1182.0000 1322.0000 31191.0000
K 1188.0000 1152.0000 31623.0000
B 1208.0000 801.0000 24279.0000
G 1223.0000 978.0000 77198.0000
I 1340.0000 1077.0000 41864.0000
M 1342.0000 1258.0000 31596.0000
D 1376.0000 868.0000 31524.0000
'''

select_psfs_without_file_parameters = [(False, 21, all_found_psfs_without_file),
                                       (True, 18, all_found_psfs_without_file_ga)]
@pytest.mark.parametrize('ga, n_psfs, correct_all_found_psfs_txt', select_psfs_without_file_parameters)
def test_select_psfs_without_file(test_directory, ga, n_psfs, correct_all_found_psfs_txt):
    guiding_selections_path, all_found_psfs_path = select_psfs(
        CONVERTED_NIRCAM_IM, ROOT, 2, global_alignment=ga, testing=True,
        out_dir=test_directory
    )

    # Ensure the correct files were generated
    assert os.path.exists(guiding_selections_path), 'guiding_selections_test_select_psfs_G2.txt not generated.'
    assert os.path.exists(all_found_psfs_path), 'all_found_psfs_test_select_psfs_G2.txt not generated.'

    # Test that the contents of all_found_psfs_test_select_psfs_G2.txt is correct
    # Right number of PSFs found?
    all_segment_locations = asc.read(all_found_psfs_path)
    assert len(all_segment_locations) == n_psfs

    # Right file contents
    with open(all_found_psfs_path) as f:
        all_found_psfs_contents = f.read()
    assert all_found_psfs_contents == correct_all_found_psfs_txt
