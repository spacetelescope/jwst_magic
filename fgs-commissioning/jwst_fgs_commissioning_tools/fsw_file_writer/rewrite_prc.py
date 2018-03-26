'''Make custom prc files'''
from astropy.io import ascii as asc
import numpy as np

from jwst_fgs_commissioning_tools.star_selector import select_psfs
from jwst_fgs_commissioning_tools.fsw_file_writer.mkproc import Mkproc
from jwst_fgs_commissioning_tools import utils


def rewrite_prc(order, guider, root, out_dir, thresh_factor=0.9,
                prc=True, regfile=True):
    '''For a given dataset, rewrite the PRC and regfile to select a new commanded
    guide star and reference stars'''

    # Open log of all identified PSFs
    print("Reading from (and writing to) {}".format(out_dir))
    all_psfs = out_dir + '/{}_G{}_ALLpsfs.txt'.format(root, guider)
    all_rows = asc.read(all_psfs)
    labels = all_rows['label'].data

    if len(set(labels)) != len(labels):
        raise ValueError('Could not accurately map labels to segments. This is '
                         'most likely due to incorrect labelling in the ALLpsfs.txt file.')

    if not order.isalpha():
        raise TypeError('Must enter only letters as PSF labels.')

    # Match each segment to its label
    inds = [np.argwhere(labels == letter)[0] for letter in order]
    inds = np.concatenate(inds)
    rows = all_rows[inds]

    # Rewrite CECIL proc file
    if prc:
        # Add threshold
        thresh = rows['countrate'] * thresh_factor
        rows['threshold'] = thresh
        print('Threshold: {}'.format(thresh_factor))

        Mkproc(guider, root, rows['x'], rows['y'], rows['countrate'], step='ID',
               out_dir=out_dir, thresh_factor=thresh_factor)

    # Rewrite regfile.txt
    if regfile:
        cols = select_psfs.create_cols_for_coords_counts(all_rows['x'], all_rows['y'],
                                                         all_rows['countrate'], 0,
                                                         inds=inds)
        utils.write_cols_to_file(out_dir,
                                 filename='{0}_G{1}_regfile.txt'.format(root, guider),
                                 labels=['y', 'x', 'countrate'],
                                 cols=cols)
