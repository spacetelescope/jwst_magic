"""Rewrite .prc and regfile.txt files with new guide and reference stars

Allow the user to rewrite the ID.prc, ACQ.prc, and regfile.txt files
without rerunning image conversion and star selection.

Authors
-------
    - Lauren Chambers

Use
---
    This module can be used as such:
    ::
        from jwst_magic.fsw_file_writer import rewrite_prc
        rewrite_prc.rewrite_prc(inds, guider, root, out_dir)

    Required arguments:
        ```inds``` - string listing segment labels or list of segment indices
            indicating which segments to re-write as the guide and
            reference stars
        ```guider``` - guider number (1 or 2)
        ```root``` - name used to create the output directory,
            {out_dir}/out/{root}
        ```out_dir``` - where output files will be saved. If not provided, the
            image(s) will be saved within the repository at
            tools/fgs-commissioning/

    Optional arguments:
        ```thresh_factor``` - factor by which to multiply the countrates
            to determine the threshold count rate
        ```prc``` - denotes whether to rewrite {root}_G{guider}_ID.prc
        ```regfile``` - denotes whether to rewrite {root}_G{guider}_regfile.txt

"""

# Standard Library Imports
import os
import logging

# Third Party Imports
from astropy.io import ascii as asc
import numpy as np

# Local Imports
from .. import utils
from ..star_selector import select_psfs
from ..fsw_file_writer.mkproc import Mkproc

# Start logger
LOGGER = logging.getLogger(__name__)


def rewrite_prc(inds, guider, root, out_dir, thresh_factor=0.9,
                prc=True, regfile=True):
    """For a given dataset, rewrite the PRC and regfile to select a
    new commanded guide star and reference stars

    Parameters
    ----------
    inds : str or list
        String listing segment labels or list of segment indices
        indicating which segments to re-write as the guide and
        reference stars
    guider : int
        Guider number (1 or 2)
    root : str
        Name used to create the output directory, {out_dir}/out/{root}
    out_dir : str
        Where output files will be saved. If not provided, the
        image(s) will be saved within the repository at
        tools/fgs-commissioning/
    thresh_factor : float, optional
        Factor by which to multiply the countrates to determine
        the threshold count rate
    prc : bool, optional
        Denotes whether to rewrite {root}_G{guider}_ID.prc
    regfile : bool, optional
        Denotes whether to rewrite {root}_G{guider}_regfile.txt

    Raises
    ------
    TypeError
        A string of labels was provided that is not only letters.
    ValueError
        Could not match provided segment labels to ALLpsfs.txt.
        Value passed to inds that is neither an alphabetic
            string or a list of integers.
    """

    out_dir = os.path.join(out_dir, 'out', root)

    # Open log of all identified PSFs
    LOGGER.info("Rewrite PRC: Reading from (and writing to) {}".format(out_dir))
    all_psfs = os.path.join(out_dir, '{}_G{}_ALLpsfs.txt'.format(root, guider))
    all_rows = asc.read(all_psfs)

    # If the segments are specified as a string of alphabetic labels
    if isinstance(inds, str):
        labels = all_rows['label'].data
        if len(set(labels)) != len(labels):
            raise ValueError('Could not accurately map labels to segments. This is '
                             'most likely due to incorrect labelling in the ALLpsfs.txt file.')
        if not inds.isalpha():
            raise TypeError('Must enter only letters as PSF labels.')

        # Match each segment to its label
        inds = [np.argwhere(labels == letter)[0] for letter in inds]
        inds = np.concatenate(inds)
    # If the segments are specified as a list of inds, directly from the GUI
    elif isinstance(inds, list):
        pass
    else:
        raise ValueError('Invalid indices provided (must be string of alphabetic labels or list of integer indices): ', inds)

    rows = all_rows[inds]

    # Rewrite CECIL proc file
    if prc:
        # Add threshold
        thresh = rows['countrate'] * thresh_factor
        rows['threshold'] = thresh
        LOGGER.info('Rewrite PRC: Threshold: {}'.format(thresh_factor))

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
