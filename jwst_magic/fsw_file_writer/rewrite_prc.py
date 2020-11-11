"""Rewrite .prc and guiding_selections*.txt files with new guide and reference stars

Allow the user to rewrite the ID.prc, ACQ.prc, and guiding_selections*.txt files
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
        ``inds`` - string listing segment labels or list of segment indices
            indicating which segments to re-write as the guide and
            reference stars
        ``guider`` - guider number (1 or 2)
        ``root`` - name used to create the output directory,
            {out_dir}/out/{root}
        ``out_dir`` - where output files will be saved. If not provided, the
            image(s) will be saved within the repository at
            jwst_magic/

    Optional arguments:
        ``thresh_factor`` - factor by which to multiply the countrates
            to determine the threshold count rate
        ``prc`` - denotes whether to rewrite {root}_G{guider}_ID.prc
        ``guiding_selections_file`` - denotes whether to rewrite
            guiding_selections_{root}_G{guider}.txt

"""

# Standard Library Imports
from collections import OrderedDict
import io
import logging
import os
import yaml

# Third Party Imports
from astropy.io import ascii as asc
import numpy as np

# Local Imports
from jwst_magic.fsw_file_writer import buildfgssteps, write_files
from jwst_magic.star_selector import select_psfs
from jwst_magic.utils import utils

# Start logger
LOGGER = logging.getLogger(__name__)


def rewrite_prc(inds_list, segnum, guider, root, out_dir, threshold, shifted, crowded_field):
    """For a given dataset, rewrite the PRC and guiding_selections*.txt to select a
    new commanded guide star and reference stars

    Parameters
    ----------
    inds_list : list
        List of configurations of segment indices indicating which
        segments to re-write as the guide and reference stars
    segnum : int
        Segment number for the chosen center of pointing
    guider : int
        Guider number (1 or 2)
    root : str
        Name used to create the output directory, {out_dir}/out/{root}
    out_dir : str
        Where output files will be saved. If not provided, the
        image(s) will be saved within the repository at
        jwst_magic/
    threshold : float
        Threshold to use in prc files
    shifted : bool
        If the image has been chosen to be shifted to the ID attitude
    crowded_field : bool
        If the image is a crowded field image

    Raises
    ------
    TypeError
        A string of labels was provided that is not only letters.
    ValueError
        Could not match provided segment labels to all_found_psfs*.txt.
        Value passed to inds that is neither an alphabetic
            string or a list of integers.
    """

    out_path = os.path.join(out_dir, 'out', root)
    file_root = '{}_G{}'.format(root, guider)
    LOGGER.info("Rewrite PRC: Reading from (and writing to) {}".format(out_path))

    # Find converted FGS image
    if os.path.exists(os.path.join(out_path, 'FGS_imgs/unshifted_{}.fits'.format(file_root))):
        fgs_im = os.path.join(out_path, 'FGS_imgs/unshifted_{}.fits'.format(file_root))
    elif os.path.exists(os.path.join(out_path, 'FGS_imgs/{}.fits'.format(file_root))):
        fgs_im = os.path.join(out_path, 'FGS_imgs/{}.fits'.format(file_root))
    else:
        LOGGER.error("Cannot find expected raw FGS image input. Looked in {} for a file named "
                     "unshifted_{}.fits.".format(out_path, file_root))

    # Find psf_center file
    psf_center_file = os.path.join(out_path, 'unshifted_psf_center_{}.txt'.format(file_root))
    if not os.path.exists(psf_center_file):
        LOGGER.warning("Cannot find psf_center file: {}. Assuming it is not relevant to this image type".format(
            os.path.join(out_path, psf_center_file)))
        psf_center_file = None

    # Read in unshifted all psfs file(s)
    all_psfs_unshifted = os.path.join(out_path, 'unshifted_all_found_psfs_{}.txt'.format(file_root))
    all_psfs_old = os.path.join(out_path, 'all_found_psfs_{}.txt'.format(file_root))
    if os.path.exists(all_psfs_unshifted):
        unshifted_all_psfs = all_psfs_unshifted
    elif os.path.exists(all_psfs_old):
        unshifted_all_psfs = all_psfs_old
    else:
        LOGGER.warning('Rewrite PRC: Cannot find all found PSFs file in this directory: {}. '
                       'Looked for {} and {}.'.format(out_path, all_psfs_unshifted, all_psfs_old))
    LOGGER.info("Rewrite PRC: Reading unshifted all psfs file: {}".format(all_psfs_unshifted.split('/')[-1]))
    unshifted_all_rows = asc.read(unshifted_all_psfs)

   # Rewrite new unshifted guiding selections file, with new selections
    cols_list = [select_psfs.create_cols_for_coords_counts(unshifted_all_rows['x'], unshifted_all_rows['y'],
                                                     unshifted_all_rows['countrate'], 0,
                                                     inds=inds) for inds in inds_list]

    if segnum is not None:
        # Write out center of pointing information
        center_pointing_path = os.path.join(out_path, 'center_pointing_{}_G{}.txt'.format(root, guider))
        utils.write_cols_to_file(center_pointing_path, labels=['segnum'], cols=[segnum], log=LOGGER)

    # Write guiding selections file(s)
    current_dirs = sorted([int(d.split('guiding_config_')[-1]) for d in os.listdir(out_path)
                           if os.path.isdir(os.path.join(out_path, d)) if 'guiding_config' in d])
    new_config_numbers = np.arange(max(current_dirs)+1, max(current_dirs)+1+len(inds_list))
    guiding_selections_path_list = []
    for (i, cols) in zip(new_config_numbers, cols_list):
        guiding_selections_path = os.path.join(out_path, 'guiding_config_{}'.format(i),
                                               'unshifted_guiding_selections_{}_G{}_config{}.txt'.format(
                                                   root, guider, i))
        utils.write_cols_to_file(guiding_selections_path,
                                 labels=['y', 'x', 'countrate'],
                                 cols=cols, log=LOGGER)
        guiding_selections_path_list.append(guiding_selections_path)


    # Shift and write out FSW files
    for i, guiding_selections_file in enumerate(guiding_selections_path_list):
        # Change out_dir to write data to guiding_config_#/ sub-directory next to the selections file
        if 'guiding_config' in guiding_selections_file:
            out_dir_fsw = os.path.join(out_path, 'guiding_config_{}'.format(
                guiding_selections_file.split('guiding_config_')[1].split('/')[0]))

        if shifted:
            fgs_im_fsw, guiding_selections_file_fsw, psf_center_file_fsw = buildfgssteps.shift_to_id_attitude(
                fgs_im, root, guider, out_dir_fsw, guiding_selections_file=guiding_selections_file,
                all_found_psfs_file=unshifted_all_psfs, psf_center_file=psf_center_file,
                crowded_field=crowded_field, logger_passed=True)
        else:
            fgs_im_fsw = fgs_im
            guiding_selections_file_fsw = guiding_selections_file
            psf_center_file_fsw = psf_center_file

        # Rewrite CECIL proc file
        for step in ['ID', 'ACQ1', 'ACQ2', 'TRK']:
            fgs_files_obj = buildfgssteps.BuildFGSSteps(
                fgs_im_fsw, guider, root, step, out_dir=out_dir_fsw, threshold=threshold,
                logger_passed=True, guiding_selections_file=guiding_selections_file_fsw,
                psf_center_file=psf_center_file_fsw, shift_id_attitude=shifted,
            )

            filename_root = '{}_G{}_{}'.format(root, guider, step)
            fgs_files_obj.filename_root = filename_root
            if step in ['ID', 'ACQ1']:
                write_files.write_prc(fgs_files_obj)
            if step != 'ID':
                write_files.write_image(fgs_files_obj)
        LOGGER.info("*** Finished FSW File Writing for Selection #{} ***".format(i + 1))

    LOGGER.info("*** FSW File Writing: COMPLETE ***")
