"""Find all the relevant PSFs in the image, be it manually or using a
file, and generate a guiding_selections*.txt and all_found_psfs*.txt.

Analyze image data using the photutils find_peaks function to identify
the locations and count rates of all stars (or segments) in the image.
Either prompt the use to select which stars to use as the guide and
reference stars, or read the guide and reference stars from a
pre-existing guiding_selections*.txt. Generate an all_found_psfs*.txt file that lists all the
segments in the image, and a guiding_selections*.txt file that lists just the guide
and reference stars.

Authors
-------
    - Keira Brooks
    - Lauren Chambers

Use
---
    This module can be used as such:
    ::
        from jwst_magic.star_selector import select_psfs
        select_psfs.create_reg_file(fgs_im, root, guider)

    Required arguments:

"""

# Standard Library Imports
import io
from collections import OrderedDict
import glob
import logging
import os
import random
import shutil
import yaml

# Third Party Imports
from astropy.io import ascii as asc
from astropy.io import fits
import matplotlib
JENKINS = '/home/developer/workspace/' in os.getcwd()
if matplotlib.get_backend() != 'Qt5Agg' and not JENKINS:
    matplotlib.use('Qt5Agg')  # Make sure that we are using Qt5
from matplotlib import rcParams
from matplotlib.colors import LogNorm
import matplotlib.pyplot as plt
import numpy as np
from photutils import find_peaks
from scipy import ndimage

# Local Imports
if not JENKINS:
    from jwst_magic.star_selector import SelectStarsGUI
from jwst_magic.utils import utils

# Adjust matplotlib parameters
rcParams['image.origin'] = 'upper'
rcParams['font.family'] = 'serif'
rcParams['font.weight'] = 'light'
rcParams['mathtext.bf'] = 'serif:normal'

# Paths
SELECT_PSFS_PATH = os.path.dirname(os.path.realpath(__file__))
PACKAGE_PATH = os.path.split(SELECT_PSFS_PATH)[0]
OUT_PATH = os.path.split(PACKAGE_PATH)[0]  # Location of out/ and logs/ directory
DATA_PATH = os.path.join(PACKAGE_PATH, 'data')

# Start logger
LOGGER = logging.getLogger(__name__)


def plot_centroids(data, coords, root, guider, out_dir):
    """Plot the identified segment locations over the data

    Parameters
    ----------
    data : 2-D numpy array
        Image data
    coords : list
        List of tuples of x and y coordinates of all identified PSFs
    root : str, optional
        Name used to generate output folder and output filenames.
    guider : int
        Guider number (1 or 2)
    out_dir : str
        Where output files will be saved.
    """
    pad = 300

    # Determine x and y limits that encompass all PSFS
    xarray, yarray = [x for (x, y) in coords], [y for (x, y) in coords]
    x_mid = (min(xarray) + max(xarray)) / 2
    y_mid = (min(yarray) + max(yarray)) / 2
    x_range = max(xarray) - min(xarray)
    y_range = max(yarray) - min(yarray)
    ax_range = max(x_range, y_range)  # Choose the larger of the dimensions
    ax_range += 100  # Make sure not to clip off the edge of border PSFS

    plt.figure(figsize=(17, 17))
    plt.imshow(data, cmap='Greys', norm=LogNorm())
    for i in range(len(coords)):
        plt.scatter(coords[i][0], coords[i][1])
        plt.annotate('({}, {})'.format(int(coords[i][0]), int(coords[i][1])),
                     (coords[i][0] - (pad * 0.05),
                      coords[i][1] + (pad * 0.05)))
    plt.title('Centroids found for {}'.format(root))
    plt.xlim(max(0, x_mid - ax_range / 2), min(2048, x_mid + ax_range / 2))
    plt.ylim(min(2048, y_mid + ax_range / 2), max(0, y_mid - ax_range / 2))
    plt.savefig(os.path.join(out_dir, 'all_found_psfs_centroids_{}_G{}.png'.format(root, guider)))

    plt.close()


def parse_in_file(in_file):
    """Get the position and countrates from a provided guiding_selections*.txt

    Determines if the input file contains x, y, and countrate data. If
    so, extracts the locations and countrates of the stars accordingly.
    Recognizes columns named "x" or "xreal"; "y" or "yreal";
    "countrate", "count rate", or "ctot".

    Parameters
    ----------
    in_file : str
        File containing locations and count rates of selected segments

    Returns
    -------
    out_cols : list
        List of positions and countrates of selected segments
    coords : list
        List of tuples with X and Y positions of selected PSFs
    nref : int
        The number of selected reference stars

    Raises
    ------
    TypeError
        Incompatible file provided to in_file
    """

    in_table = asc.read(in_file)
    colnames = in_table.colnames

    # Handle old regfiles where countrate column is titled "count rate" and thus
    # astropy can't automatically identify it as a single column name
    if 'col1' in colnames:
        # If astropy doesn't know how to read the column names...
        fo = open(in_file)
        raw_columns = fo.readlines()[0]
        fo.close()

        if 'x' and 'y' and 'count rate' in raw_columns:
            # If this is formatted like an old regfile, just fix it and continue
            raw_columns = str.replace(raw_columns, 'count rate', 'countrate')
            fix_colnames = raw_columns.split()

            if '#' in fix_colnames:
                # Remove commenting mark
                fix_colnames.pop(fix_colnames.index('#'))

            if len(fix_colnames) != len(in_table.columns):
                # Fix didn't work
                err_message = 'Unknown in_file format: {}. Expecting columns named \
                               "x"/"xreal", "y"/"yreal", and \
                               "count rate"/"countrate"/"ctot". Found columns \
                               named {}. Please rename columns.'.format(in_file,
                                                                        fix_colnames)
                raise TypeError(err_message)

            for old_col, fix_col in zip(colnames, fix_colnames):
                # Assign fixed column names to table
                in_table.rename_column(old_col, fix_col)

            colnames = fix_colnames

        else:
            # If not just an old regfile, it is something unfamiliar
            err_message = 'Unknown in_file format: {}. Expecting columns named \
                           "x"/"xreal", "y"/"yreal", and \
                           "count rate"/"countrate"/"ctot". Found columns \
                           named {}. Please rename columns.'.format(in_file,
                                                                    raw_columns.split())
            raise TypeError(err_message)

    # Make sure all the necessary columns are present
    x_check = 'x' in colnames or 'xreal' in colnames
    y_check = 'y' in colnames or 'yreal' in colnames
    countrate_check = 'countrate' in colnames or 'ctot' in colnames

    if not (x_check and y_check and countrate_check):
        err_message = 'Unknown in_file format: {}. Expecting columns named \
                       "x"/"xreal", "y"/"yreal", and \
                       "count rate"/"countrate"/"ctot". Found columns \
                       named {}. Please rename columns.'.format(in_file, colnames)
        raise TypeError(err_message)

    # Passed all the checkpoints! Move on to process the file.
    LOGGER.info('Star Selection: Checking stars from input file {}'.format(in_file))

    # Rename relevant old columns, if necessary:
    for col in colnames:
        if col == 'xreal':
            in_table.rename_column(col, 'x')
        if col == 'yreal':
            in_table.rename_column(col, 'y')
        if col == 'countrate':
            in_table.rename_column(col, 'ctot')

    # If an incat, only use stars in the catalog
    if 'cat' in colnames:
        in_table = in_table[in_table['cat'] == 1]

    # Set data types
    in_table['x'] = in_table['x'].astype(int)  # Need to be an integer
    in_table['y'] = in_table['y'].astype(int)
    in_table['ctot'] = in_table['ctot'].astype(int)

    out_cols = in_table['y', 'x', 'ctot']  # Make sure to flip x and y!!
    coords = [(x, y) for y, x in zip(in_table['y'], in_table['x'])]
    nref = len(in_table['x']) - 1

    return out_cols, coords, nref


def copy_psfs_files(guiding_selections_file_list, output_file, root, guider, out_dir):
    """By parsing the name of an input guiding_selections*.txt file,
    identify a corresponding file. See "output_file" parameter for information
    on what files are accepted.

    Parameters
    ----------
    guiding_selections_file_list : list of str
        Path(s) to unshifted_guiding_selections*.txt file
    output_file : str
        The type of file to copy and return the path of.
        Eg 'all_found_psfs', 'psf_center', 'center_pointing',
        'all_selections_yaml'
    root : str
        Name used to generate output folder and output file names.
    guider : int
        Guider number (1 or 2)
    out_dir : str
        Directory in which to store the copied all_found_psfs file.

    Returns
    -------
    copied_all_found_psfs : str
        Path to the copied all_found_psfs*.txt file
    """
    # Determine the root of the imported guiding selections file name
    filenames, roots= [], []
    for file in guiding_selections_file_list:
        filename = os.path.basename(file)
        if 'guiding_selections_' in filename:
            if 'unshifted' in filename:
                imported_root = filename.split('unshifted_guiding_selections_')[-1].split('.txt')[0]
            else:
                imported_root = filename.split('guiding_selections_')[-1].split('.txt')[0]
            if 'config' in imported_root:
                imported_root = imported_root.split('_config')[0]
        elif '_regfile' in filename:
            imported_root = filename.split('_regfile.txt')[0]
        else:
            imported_root = None
            LOGGER.warning(
                'Star Selection: Could not parse root from provided guiding selections file ({}). '.format(filename) +
                'Not able to copy over a corresponding {}*.txt file.'.format(output_file)
            )
        filenames.append(filename)
        roots.append(imported_root)

    # Remove any cases where the root is None
    for i, rt in enumerate(roots):
        if rt is None:
            roots.pop(i)
            guiding_selections_file_list.pop(i)
            filenames.pop(i)

    # If there are any guiding files with a known root left
    if len(roots) > 0:
        dirs = [os.path.dirname(file) for file in guiding_selections_file_list]
        dirs = [dr.split('guiding_config')[0] if 'guiding_config_' in dr else dr for dr in dirs]
        rootdirs = ['{}_{}'.format(root, dir) for root, dir in zip(roots, dirs)]

        # If the roots or basepaths differ, choose the first one, but log it
        if len(guiding_selections_file_list) > 1 and len(set(rootdirs)) != 1:
            LOGGER.warning('Star Selection: The multiple guiding selections files chosen do not have matching roots '
                           'and/or basepaths. We will search for the all_found_psf*.txt file using the '
                           'root/basepath of the first guiding selections file. If this is not the '
                           'intended search location, re-arrange the order of the guiding selections files.')

        # Use the info from the first guiding selections file in the list to search for the all found psfs file
        filename = filenames[0]
        imported_root = roots[0]
        dir_to_look = dirs[0]

        # Try to copy the corresponding *all_found_psfs*.txt or *psf_center*.txt file, if present.
        if imported_root is not None:
            txt_files = glob.glob(os.path.join(dir_to_look, "*.txt"))
            yaml_files = glob.glob(os.path.join(dir_to_look, "*.yaml"))
            if output_file == 'psf_center':
                acceptable_files = [
                    os.path.join(dir_to_look, 'unshifted_psf_center_{}.txt'.format(imported_root)),  # newest
                    os.path.join(dir_to_look, 'psf_center_{}.txt'.format(imported_root))]  # oldest
                ftype = 'txt'
            elif output_file == 'all_found_psfs':
                acceptable_files = [
                    os.path.join(dir_to_look, 'unshifted_all_found_psfs_{}.txt'.format(imported_root)),
                    os.path.join(dir_to_look, 'all_found_psfs_{}.txt'.format(imported_root)),
                    os.path.join(dir_to_look, '{}_ALLpsfs.txt'.format(imported_root))]
                ftype = 'txt'
            elif output_file == 'center_pointing':
                acceptable_files = [
                    os.path.join(dir_to_look, 'center_pointing_{}.txt'.format(imported_root))
                ]
                ftype = 'txt'
            elif output_file == 'all_selections_yaml':
                acceptable_files = [
                    os.path.join(dir_to_look, 'all_guiding_selections.yaml')
                ]
                ftype = 'yaml'
            if ftype == 'txt':
                try:
                    file_to_copy = [f for f in acceptable_files if f in txt_files][0]
                except IndexError:
                    file_to_copy = None
            elif ftype == 'yaml':
                try:
                    file_to_copy = [f for f in acceptable_files if f in yaml_files][0]
                except IndexError:
                    file_to_copy = None

        if output_file == 'all_selections_yaml':
            copied_psfs_file = os.path.join(out_dir, 'all_guiding_selections.yaml')
            guiding_selections_file_list = copy_all_selections_yaml(file_to_copy, copied_psfs_file,
                                                                    guiding_selections_file_list, out_dir)
            return file_to_copy, guiding_selections_file_list

        if file_to_copy is not None:
            try:
                if output_file == 'center_pointing':
                    copied_psfs_file = os.path.join(out_dir, '{}_{}_G{}.txt'.format(output_file, root, guider))
                    shutil.copy(file_to_copy, copied_psfs_file)
                else:
                    copied_psfs_file = os.path.join(out_dir, 'unshifted_{}_{}_G{}.txt'.format(output_file, root, guider))
                    shutil.copy(file_to_copy, copied_psfs_file)
                LOGGER.info('Star Selection: Copying over {}'.format(file_to_copy))
                LOGGER.info('Star Selection: Successfully wrote: {}'.format(copied_psfs_file))
                return copied_psfs_file
            except shutil.SameFileError:
                LOGGER.info("Star Selection: {} file already exists in the right location.".format(output_file))
                return file_to_copy
        else:
            if output_file == 'psf_center':
                LOGGER.info(
                    'Star Selection: Could not find a corresponding unshifted_{}*.txt file for '
                    'the provided guiding selections file ({}) in {}. Assuming it is not relevant '
                    'to this image type'.format(output_file, filename, dir_to_look)
                )
                return file_to_copy
            else:
                if output_file == 'center_pointing':
                    psf_file = os.path.join(out_dir, '{}_{}_G{}.txt'.format(output_file, root, guider))
                elif output_file == 'all_found_psfs':
                    psf_file = os.path.join(out_dir, 'unshifted_{}_{}_G{}.txt'.format(output_file, root, guider))

                if os.path.isfile(psf_file):
                    LOGGER.warning(
                        'Star Selection: Could not find a corresponding {}*.txt file for '
                        'the provided guiding selections file ({}) in {} to copy over. Instead '
                        'MAGIC will use the file that is already present: {}'.format(output_file, filename,
                                                                                     dir_to_look, psf_file)
                    )
                    return psf_file

                else:
                    LOGGER.warning('MAGIC cannot find a(n) {} file, either in the location of the '
                                 'input guiding_selections file or in the root directory. This '
                                 'file may be missing, or it may have an incompatible name. '
                                 'Looking for a file named {}. Fix '
                                 'this file before continuing.'.format(output_file, psf_file))


def copy_all_selections_yaml(file_to_copy, final_file, guiding_selections_file_list, out_dir):
    """Handle copying the config information to the yaml file"""
    utils.setup_yaml()

    # Open yaml file to pull config info from
    try:
        with open(file_to_copy, 'r') as stream:
            data_loaded = yaml.safe_load(stream)
    except TypeError:
        data_loaded = None

    # Open or create a yaml file to copy to
    if os.path.exists(final_file):
        try:
            with open(final_file, 'r') as stream:
                final_data = OrderedDict(yaml.safe_load(stream))
            config = int(sorted(final_data.keys(), key=utils.natural_keys)[-1].split('_')[-1]) + 1
        except TypeError:
            final_data = OrderedDict()
            config = 1
    else:
        final_data = OrderedDict()
        config = 1

    # For each guiding selections file passed in
    for i, file in enumerate(guiding_selections_file_list):
        # If we can pull the config number and have a yaml file to check the config number against
        if '/guiding_config_' in file and data_loaded is not None:

            # Pull the config number
            old_config = file.split('/guiding_config_')[-1].split('/')[0]

            # Pull the corresponding config from the yaml file
            config_data = data_loaded['guiding_config_{}'.format(old_config)]

            # Check if the data matches an existing config (rs order doesn't matter)
            chosen_gs = config_data[0]
            chosen_ref = config_data[1:]

            # Pull keys where guide star matches
            key_matching_gs = [key for key, value in final_data.items() if len(value) != 0 and value[0] == chosen_gs]

            # If there's no existing config with a matching guide star, it's a new config
            if len(key_matching_gs) == 0:
                final_data['guiding_config_{}'.format(config)] = config_data
                config += 1
                continue

            # Check these keys for matching ref stars (order doesn't matter)
            l = [True if set(final_data[key][1:]) == set(chosen_ref) else False for key in key_matching_gs]

            # If the current selection has not already been made
            if True not in l:
                final_data['guiding_config_{}'.format(config)] = config_data
                config += 1
            else:
                # Pull key that matches where the config is already saved
                key = key_matching_gs[l.index(True)]

                # Overwrite the guiding_selections_file_list with the path to the matching config
                repeated_file = glob.glob(os.path.join(out_dir, key, "unshifted_guiding_selections*.txt"))[0]
                guiding_selections_file_list[i] = repeated_file
                LOGGER.warning('Guiding selections from file {} match the selections already made and stored in {}/. '
                            'This file will not be used.'.format(file, key))
                pass

        else:
            final_data['guiding_config_{}'.format(config)] = []
            config += 1
            LOGGER.info('Cannot parse {} for a config number. This configuration will not be '
                        'added to all_guiding_selections.yaml file.'.format(file))

    # Write to file
    with io.open(final_file, 'w', encoding="utf-8") as f:
        yaml.dump(final_data, f, default_flow_style=False, allow_unicode=True)

    LOGGER.info('Star Selection: Copying over {}'.format(file_to_copy))
    LOGGER.info('Star Selection: Successfully wrote: {}'.format(final_file))

    # Return new list with any duplicate files deleted
    return guiding_selections_file_list


def manual_star_selection(data, all_found_psfs_path, guider,
                          out_dir, choose_center=False, testing=False, masterGUIapp=None):
    """Launches a GUI to prompt the user to click-to-select guide and
    reference stars.

    Algorithmically find and locate all PSFs in image using
    photutils.find_peaks; prompt user to select guide and reference stars
    using the GUI.

    Parameters
    ----------
    data : 2-D numpy array
        Image data
    all_found_psfs_path : str
        Path to the unshifted_all_found_psfs_{root}_G{guider}.txt file, which
        the locations and count rates of all segments found in the data.
    out_dir : str
        Where output files will be saved. If not provided, the
        image(s) will be saved within the repository at
        jwst_magic/. This path is the level outside the out/root/ dir
    choose_center : bool
        Automatically choose the one, highly-smoothed, PSF found in the image
    testing : bool, optional
        Generates indices randomly (for running pytests)
    masterGUIapp : qApplication, optional
        qApplication instance of parent GUI

    Returns
    -------
    cols_list : list of lists
        List of positions and countrates of selected segments, 1 sub-list per guiding config
    coords : list
        List of tuples with X and Y positions of found PSFs
    nref : list of int
        The number of selected reference stars, 1 int per guiding config
    all_cols : list
        List of positions and countrates of all segments in the image

    Raises
    ------
    ValueError
        The user closed the GUI without selecting any stars.
    """
    read_table = asc.read(all_found_psfs_path)
    x = read_table['x']
    y = read_table['y']
    countrate = read_table['countrate']
    num_psfs = len(x)
    coords = list(zip(x, y))

    # Find the minimum distance between PSFs
    if len(coords) < 2:
        # For cases where we only have star, we assume that we are sufficiently
        # isolated from other stars, but also that the guide star's PSF may be
        # distorted enough that it might appear quite large on the detector
        dist = 20
    else:
        dist = np.floor(np.min(utils.find_dist_between_points(coords))) - 1.

    # Call the GUI to pick PSF indices
    if not testing and not choose_center:
        gui_data = data.copy()
        gui_data[data == 0] = 1  # Alter null pixel values for LogNorm imshow
        inds_list, center_of_pointing = SelectStarsGUI.run_SelectStars(gui_data, x, y, dist, guider,
                                                                       out_dir=out_dir,
                                                                       print_output=False,
                                                                       masterGUIapp=masterGUIapp)

        # Print indices of each guiding configuration
        for i in range(len(inds_list)):
            ind = inds_list[i]
            LOGGER.info('Star Selection: Guiding Configuration {} - GS = {}, RS = {}'.format(i+1, ind[0],
                        ', '.join([str(c) for c in ind[1:]])))

    # If in testing mode, just make a random list of indices
    else:
        # Make random list of inds
        n_select = min(11, num_psfs)
        center_of_pointing = 0
        LOGGER.info('Star Selection: Testing mode; selecting {} PSFs at random.'.format(n_select))
        inds_list = [random.sample(range(num_psfs), n_select)] # a single guiding config

    nref_list = [len(inds) - 1 for inds in inds_list]
    if len(inds_list) == 0:
        raise ValueError('Star Selection: No guide star and no reference stars selected')
    else:
        for i, nref in enumerate(nref_list):
            LOGGER.info('Star Selection: Config {}: 1 guide star and {} reference stars selected'.format(i+1, nref))

    cols_list = [utils.create_cols_for_coords_counts(x, y, countrate, inds=inds) for inds in inds_list]

    return cols_list, coords, nref_list, center_of_pointing


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# MAIN FUNCTION
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


def select_psfs(data, root, guider, all_found_psfs_path, guiding_selections_file_list=None,
                psf_center_path=None, smoothing='default', choose_center=False,
                testing=False, out_dir=None, masterGUIapp=None, logger_passed=False):
    """Select guide and reference segments.

    Locate all of the segments in the provided data, then either parse a
    catalog of sources or prompt the user with a GUI to select which
    segments as the guide and reference segments. Generates one file:
        unshifted__guiding_selections_{root}_G{guider}.txt
            Lists the locations and count rates of the selected
            guide and reference segments

    Parameters
    ----------
    data : 2-D numpy array or str
        Image data
    root : str
        Name used to generate output folder and output filenames.
    guider : int
        Guider number (1 or 2)
    all_found_psfs_path : str
        Path to the unshifted_all_found_psfs_{root}_G{guider}.txt file,
        which the locations and count rates of all segments found in the
        data. File was written out in the convert_image section of MAGIC
    guiding_selections_file_list : list of str, optional
        List of files containing locations and count rates of selected segments
    psf_center_path : str
        Path to psf center file. Only needed for MIMF cases.
    smoothing: str, optional
        Options are "low" for minimal smoothing (e.g. MIMF), "high" for large
        smoothing (e.g. GA), or "default" for medium smoothing for other cases
    choose_center : bool, optional
        Automatically choose the one PSF found in the image
    testing : bool, optional
        Randomly select guide and reference stars (for use with pytests)
    out_dir : str, optional
        Where output files will be saved. If not provided, the
        image(s) will be saved within the repository at
        jwst_magic/. This path is the level outside the out/root/ dir
    masterGUIapp : qApplication, optional
        qApplication instance of parent GUI
    logger_passed : bool, optional
        Denotes if a logger object has already been generated.

    Returns
    -------
    guiding_selections_path : list of str
        List of paths to the guiding_config_#/unshifted_guiding_selections_{root}_G{guider}.txt
        file, which contains the locations and count rates of the selected guide and
        reference segments
    all_found_psfs_path : str
        Path to the unshifted_all_found_psfs_{root}_G{guider}.txt file, which
        the locations and count rates of all segments found in the data
    center_pointing_path : str
        Path to center_pointing_{root}_G{guider}.txt which contains the infomration on the
        center of pointing, either an int (0=mean, #=seg number) or a (y,x) location
    psf_center_path : str
        Path to unshifted_psf_center_{root}_G{guider}.txt file, which is only written for
        smoothing='low', MIMF case.
    """
    if not logger_passed:
        utils.create_logger_from_yaml(__name__, root=root, level='DEBUG')

    try:
        out_dir = utils.make_out_dir(out_dir, OUT_PATH, root)
        utils.ensure_dir_exists(out_dir)

        # Read in image (check if it is a filename)
        if isinstance(data, str):
            data = fits.getdata(data)

        # Pull current config numbers in root
        current_dirs = sorted([int(d.split('guiding_config_')[-1]) for d in os.listdir(out_dir)
                               if os.path.isdir(os.path.join(out_dir, d)) if 'guiding_config' in d])

        if guiding_selections_file_list:  # will be a list of strings
            # Determine the kind of in_file and parse out the PSF locations and
            # countrates accordingly
            LOGGER.info(
                "Star Selection: Reading guide and reference star positions from {}"
                    .format(', '.join(guiding_selections_file_list))
            )
            # Check for duplicate configs and write out yaml file
            yaml_path, new_guiding_selections = copy_psfs_files(guiding_selections_file_list, 'all_selections_yaml',
                                                                 root, guider, out_dir)

            # Remove duplicate files before saving data
            guiding_selections_file_list = [f for f in new_guiding_selections
                    if len([g for g in current_dirs if os.path.join(out_dir, 'guiding_config_{}/'.format(g))
                    in f]) == 0]

            # Return early if all loaded files are duplicates of existing files
            if len(guiding_selections_file_list) == 0:
                all_found_psfs_path = os.path.join(out_dir, 'unshifted_all_found_psfs_{}_G{}.txt'.format(root, guider))
                center_pointing_path = os.path.join(out_dir, 'center_pointing_{}_G{}.txt'.format(root, guider))
                if smoothing == 'low':
                    psf_center_path = os.path.join(out_dir, 'unshifted_psf_center_{}_G{}.txt'.format(root, guider))
                else:
                    psf_center_path = None
                return new_guiding_selections, all_found_psfs_path, center_pointing_path, psf_center_path

            cols_list, nref_list = [], []  # coords will be overwritten, but they should include the same data each time
            for file in guiding_selections_file_list:
                cols, coords, nref = parse_in_file(file)
                cols_list.append(cols)
                nref_list.append(nref)

            # Copy over corresponding all_found_psfs, psf_center, and center_pointing file, if possible.
            all_found_psfs_path = copy_psfs_files(guiding_selections_file_list, 'all_found_psfs', root, guider, out_dir)
            _, all_coords, _ = parse_in_file(all_found_psfs_path)

            psf_center_path = copy_psfs_files(guiding_selections_file_list, 'psf_center', root, guider, out_dir)

            center_of_pointing = None
            center_pointing_path = copy_psfs_files(guiding_selections_file_list, 'center_pointing', root, guider, out_dir)

            # Determine if the guiding file loaded already exists in the right dir - no need for a new config #
            old_configs = [True if os.path.join(out_dir, 'guiding_config_') in path else False for path in
                           guiding_selections_file_list]

        else:
            # If no .incat or reg file provided, create reg file with manual
            # star selection using the SelectStarsGUI
            cols_list, all_coords, nref_list, center_of_pointing = manual_star_selection(data,
                                                                                         all_found_psfs_path,
                                                                                         guider,
                                                                                         out_dir,
                                                                                         choose_center,
                                                                                         testing,
                                                                                         masterGUIapp)
            old_configs = [False] * len(cols_list)

        # Save PNG of image and all PSF locations in out_dir
        if not JENKINS:
            plot_centroids(data, all_coords, root, guider, out_dir)  # coords are in (x,y)

        if center_of_pointing is not None:
            # Write out center of pointing information
            center_pointing_path = os.path.join(out_dir, 'center_pointing_{}_G{}.txt'.format(root, guider))
            utils.write_cols_to_file(center_pointing_path, labels=['center_of_pointing'], cols=[center_of_pointing],
                                     log=LOGGER)

        # Determine config numbers for selections (may re-use old config if loading a file from this out_dir/root)
        new_config_numbers = []
        j = 1
        for i, config in enumerate(old_configs):
           if config is True:
               num = int(guiding_selections_file_list[i].split('guiding_config_')[-1].split('/')[0])
               new_config_numbers.append(num)
           else:
               num = len(current_dirs) + j
               j += 1
               new_config_numbers.append(num)

        # Write catalog of selected PSFs
        guiding_selections_path_list = []
        for (i, cols) in zip(new_config_numbers, cols_list):
            for j in range(len(cols)):
                if j == 0:
                    LOGGER.info("Star Selection: PSF Config {}: Guide Star at "
                                "({},{}) has 3x3 Count Rate of {}".format(i, cols[j][0], cols[j][1], cols[j][2]))
                else:
                    LOGGER.info("Star Selection: PSF Config {}: Ref Star at "
                                "({},{}) has 3x3 Count Rate of {}".format(i, cols[j][0], cols[j][1], cols[j][2]))

            guiding_selections_path = os.path.join(out_dir, 'guiding_config_{}'.format(i),
                                                   'unshifted_guiding_selections_{}_G{}_config{}.txt'.format(
                                                       root, guider, i))
            utils.write_cols_to_file(guiding_selections_path,
                                     labels=['y', 'x', 'countrate'],
                                     cols=cols, log=LOGGER)
            guiding_selections_path_list.append(guiding_selections_path)

    except Exception as e:
        LOGGER.exception(e)
        raise

    return guiding_selections_path_list, all_found_psfs_path, center_pointing_path, psf_center_path

