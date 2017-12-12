# Make custom prc files
from mkproc import Mkproc
from astropy.io import ascii as asc
import numpy as np
import time
import utils


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
    root = 'ga01_NIRCam'
    out_dir = '/Users/lchambers/TEL/FGS/Commissioning-tools/out/GA_DHAStesting/ga01_00002'

    # Prompt user to input PSFs
    order = input('Enter the order of PSFs to test with DHAS (e.g. "AQGMJ"): ')
    order = order.upper()

    rewrite_prc(order, guider, root, out_dir)