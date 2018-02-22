# Make custom prc files
from astropy.io import ascii as asc
import numpy as np
import time

from jwst_fgs_commissioning_tools.fsw_file_writer.mkproc import Mkproc
from jwst_fgs_commissioning_tools import utils

# If you wanna make a tutorial....
# If you are testing, say, a range of geometries for guide and reference stars on a single global alignment ID image in the DHAS, it would be obnoxious and cumbersome to re-run the entire tool every time you just want to command a different set of guide/reference stars. The only thing that needs to be re-written in this case is the `.prc` file; the output `.fits` image that is fed into DHAS will be the same.

# There is a stand-alone script, `DHAStesting_prc.py` that can be used in cases like this, which just re-write the `.prc` file according to user input.

# The user input is provided as a series of letters, corresponding to segments of JWST as shown below. The first letter corresponds to the guide star segment.
# <img src="JWSTgrid.png" style="width: 400px;">
# guider = 1
# root = 'tool_tutorial_ncga' # MUST match the root of the file you are trying to alter
# out_dir = os.path.join(current_dir, 'out', root) # Where to find the prc

# thresh_factor = 0.9 # default
# # The thresh_factor is the percentage of the PSF countrate that is used as DHAS threshold.
# # (e.g. if thresh_factor is 0.9 and the countrate is 100, the threshold provided to DHAS will be 90)

# order = 'OPHIA'

# DHAStesting_prc.rewrite_prc(order, guider, root, out_dir)


def rewrite_prc(order, guider, root, out_dir, thresh_factor=0.9):
    all_psfs = out_dir + '/{}_G{}_ALLpsfs.txt'.format(root, guider)

    # Open log of all identified PSFs
    print("Reading from (and writing to) {}".format(out_dir))
    all_rows = asc.read(all_psfs)
    labels = all_rows['label'].data

    if not order.isalpha():
        raise TypeError('Must enter only letters as PSF labels.')

    # Match each segment to its label
    inds = [np.argwhere(labels == letter)[0] for letter in order]
    inds = np.concatenate(inds)
    rows = all_rows[inds]

    # Add threshold
    thresh = rows['countrate'] * thresh_factor
    rows['threshold'] = thresh

    # Make CECIL proc file
    print('Threshold: {}'.format(thresh_factor))
    Mkproc(guider, root, rows['x'], rows['y'], rows['countrate'], step='ID',
           out_dir=out_dir, thresh_factor=thresh_factor)

    # # Save log file, too
    # filename='/{}_{}test_{}_{:.2f}thresh.txt'.format(time.strftime("%Y%m%d_%H%M%S"), step, order, thresh_factor)
    # asc.write(rows, output=out_dir + filename)


if __name__ == '__main__':
    guider = 1
    root = 'ga02_NIRCam'
    out_dir = '/Users/lchambers/TEL/FGS/Commissioning-tools/out/GA_DHAStesting/ga02_00001'

    # Prompt user to input PSFs
    order = input('Enter the order of PSFs to test with DHAS (e.g. "AQGMJ"): ')
    order = order.upper()

    rewrite_prc(order, guider, root, out_dir)