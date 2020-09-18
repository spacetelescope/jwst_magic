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
import logging
import os
import random
import shutil
import string

# Third Party Imports
from astropy.io import ascii as asc
from astropy.io import fits
import matplotlib
JENKINS = 'jenkins' in os.getcwd()
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


def count_psfs(smoothed_data, gauss_sigma, npeaks=np.inf, choose=False):
    """Use photutils.find_peaks to count how many PSFS are present in the data

    Parameters
    ----------
    smoothed_data : 2-D numpy array
        Image data that has been smoothed with a Gaussian filter
    gauss_sigma : float
        The sigma of the Gaussian smoothing filter
    npeaks : int or np.inf
        Number of peaks to choose with photutils.find_peaks
    choose : bool, optional
        Prompt the user to choose which method to use to select the
        threshold

    Returns
    -------
    num_psfs : int
        The number of PSFs found in the smoothed data
    coords : list
        List of tuples of x and y coordinates of all identified PSFs
    threshold : float
        The threshold used with photutils.find_peaks
    """

    if choose:
        num_psfs, coords, threshold = choose_threshold(smoothed_data, gauss_sigma)

    else:
        # Perform statistics
        median = np.median(smoothed_data)
        std = np.std(smoothed_data)

        # Find PSFs
        threshold = median + (3 * std)  # Used to be median + 3 * std

        sources = find_peaks(smoothed_data, threshold, box_size=gauss_sigma, npeaks=npeaks)
        num_psfs = len(sources)
        if num_psfs == 0:
            raise ValueError("You have no sources in your data.")
        coords = sources['x_peak', 'y_peak']
        coords = [(x, y) for [x, y] in coords]

        LOGGER.info('Star Selection: {} PSFs detected in Gaussian-smoothed data \
            (threshold = {}; sigma = {})'.format(num_psfs, threshold, gauss_sigma))

    return num_psfs, coords, threshold


def choose_threshold(smoothed_data, gauss_sigma):
    """Prompt the user to choose which method to use to select the
    threshold

    Parameters
    ----------
    smoothed_data : 2-D numpy array
        Image data that has been smoothed with a Gaussian filter
    gauss_sigma : float
        The sigma of the Gaussian smoothing filter

    Returns
    -------
    num_psfs : int
        The number of PSFs found in the smoothed data
    coords : list
        List of tuples of x and y coordinates of all identified PSFs
    threshold : float
        The threshold used with photutils.find_peaks

    Raises
    ------
    ValueError
        User did not accept either of the threshold options.
    """
    # Perform statistics
    mean = np.mean(smoothed_data)
    std = np.std(smoothed_data)

    # Run find_peaks with two different threshold options
    thresholds = [3 * std, mean]

    sources_std = find_peaks(smoothed_data, thresholds[0], box_size=gauss_sigma)
    sources_mean = find_peaks(smoothed_data, thresholds[1], box_size=gauss_sigma)

    # Show plots of each for user to choose between
    plt.ion()
    smoothed_data[smoothed_data == 0] = 0.1  # Allow LogNorm plotting
    fig, [ax1, ax2] = plt.subplots(1, 2, figsize=(16, 8))
    fig.subplots_adjust(top=.95, left=.05, bottom=.05)

    ax1.imshow(smoothed_data, cmap='bone', interpolation='nearest',
               clim=(0.1, 100), norm=LogNorm())
    ax1.scatter(sources_std['x_peak'], sources_std['y_peak'], c='r', marker='+')
    ax1.set_title('Threshold = 3$\sigma$ ({} sources found)'.format(len(sources_std)))

    ax2.imshow(smoothed_data, cmap='bone', interpolation='nearest',
               clim=(0.1, 100), norm=LogNorm())
    ax2.scatter(sources_mean['x_peak'], sources_mean['y_peak'], c='r', marker='+')
    ax2.set_title('fThreshold = Mean ({} sources found)'.format(len(sources_mean)))

    plt.get_current_fig_manager().window.raise_()
    plt.show()

    # Prompt user to choose
    choice = input('''
                   Examine the two options presented. To use the stars \
                   selected with a 3 standard deviation threshold, \
                   type "S". To use the stars selected with a mean \
                   threshold, type "M". To use neither and cancel the \
                   program, press enter.

                   Choice: ''')

    plt.close()

    if choice == 'S':
        num_psfs = len(sources_std)
        coords = [(x, y) for [x, y] in sources_std['x_peak', 'y_peak']]
        return num_psfs, coords, thresholds[0]
    if choice == 'M':
        num_psfs = len(sources_mean)
        coords = [(x, y) for [x, y] in sources_mean['x_peak', 'y_peak']]
        return num_psfs, coords, thresholds[1]
    else:
        LOGGER.error('Star Selection: User rejection of identified PSFs.')
        raise ValueError('User rejection of identified PSFs.')


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


def count_rate_total(data, objects, num_objects, x, y, countrate_3x3=True):
    """Get the count rates within each object from a segmentation image.

    Parameters
    ----------
    data : 2-D numpy array
        Image data
    objects : 2-D numpy array
        Segmentation of image data
    num_objects : int
        Number of individual objects in the segmentation data
    x : list
        List of x-coordinates of identified PSFs
    y : list
        List of y-coordinates of identified PSFs
    countrate_3x3 : bool, optional
        Calculate the value of the 3x3 square (True), or of the entire
        object (False)

    Returns
    -------
    countrate : list
        List of count rates of each segmentation object
    val : list
        List of number of pixels within each segmentation object
    """

    countrate = []
    val = []
    for i in range(1, num_objects + 1):
        im = np.copy(objects)
        im[objects != i] = False
        im[objects == i] = True

        if countrate_3x3:
            countrate.append(utils.countrate_3x3(x[i - 1], y[i - 1], np.array(data)))
        else:
            countrate.append(np.sum(im * data))
        val.append(np.sum(im * 1.))  # Number of pixels in object

    return countrate, val


def create_cols_for_coords_counts(x, y, countrate, val, labels=None, inds=None):
    """Format position and count rate data to be written to file.

    Create an array of columns of y, x, and countrate of each PSF to be
    written out. Use the inds returned from pick_stars based on user
    input. If no inds are given, put the PSF with the most compact PSF
    first in the list to make it the guide star.


    Parameters
    ----------
    x : list
        List of x-coordinates of identified PSFs
    y : list
        List of y-coordinates of identified PSFs
    countrate : list
        List of count rates of identified PSFs
    val : list
        List of the number of pixels in each PSF's segmentation object
    labels : list, optional
        Denotes whether the PSF alphabetic labels should be included as
        a column to write out
    inds : list, optional
        List of the indices of the guide and reference stars

    Returns
    -------
    cols : list
        List of segment positions, count rates, and maybe labels for
        each selected segments
    """
    if labels is not None:
        # NOTE: these coordinates are y, x
        cols = [[ll, '{:.4f}'.format(yy),
                 '{:.4f}'.format(xx),
                 '{:.4f}'.format(co)] for ll, yy, xx, co in zip(labels, y, x, countrate)]
    else:
        # NOTE: these coordinates are y, x
        cols = [[yy, xx, co] for yy, xx, co in zip(y, x, countrate)]

    if inds is None:
        min_ind = np.where(val == np.min(val))[0][0]  # Find most compact PSF
        cols.insert(0, cols.pop(min_ind))  # Move most compact PSF to top of the list
    else:
        cols = [cols[i] for i in inds]

    return cols


def match_psfs_to_segments(x, y, smoothing):
    """Match PSFs found in the image to their alphabetic label (between A and R)

    Parameters
    ----------
    x : list
        List of x-coordinates of identified PSFs
    y : list
        List of y-coordinates of identified PSFs
    smoothing: str
        Options are "low" for minimal smoothing (e.g. MIMF), "high" for large
        smoothing (e.g. GA), or "default" for medium smoothing for other cases

    Returns
    -------
    matched_labels : list
        List of alphabetic labels for each identified PSF
    """
    labels = string.ascii_uppercase[:18]

    # Determine boundaries of array
    x_min = min(x)
    x_max = max(x)
    y_min = min(y)
    y_max = max(y)

    if (x_max - x_min) > (y_max - y_min):
        # Horizontal orientation
        x_list = [1, 2, 2, 3, 3, 3, 4, 4, 5, 5, 6, 6, 7, 7, 7, 8, 8, 9]
        y_list = [3, 2, 4, 1, 3, 5, 2, 4, 1, 5, 2, 4, 1, 3, 5, 2, 4, 3]

        x_coords = np.linspace(x_min, x_max, 9)
        y_coords = np.linspace(y_min, y_max, 5)[::-1]
    else:
        # Vertical orientation
        x_list = [3, 2, 4, 1, 3, 5, 2, 4, 1, 5, 2, 4, 1, 3, 5, 2, 4, 3]
        y_list = [9, 8, 8, 7, 7, 7, 6, 6, 5, 5, 4, 4, 3, 3, 3, 2, 2, 1]

        x_coords = np.linspace(x_min, x_max, 5)
        y_coords = np.linspace(y_min, y_max, 9)[::-1]

    seg_coords = np.array([[x_coords[i_x - 1],
                            y_coords[i_y - 1]] for i_x, i_y in zip(x_list, y_list)])

    # Match actual blob coordinates to segment name
    matched_labels = []
    for x_pos, y_pos in zip(x, y):
        seg_distance = 2048
        for i_sc, sc in enumerate(seg_coords):
            x_distance = x_pos - sc[0]
            y_distance = y_pos - sc[1]
            distance = (x_distance**2 + y_distance**2)**0.5
            if distance < seg_distance:
                seg_distance = distance
                i_seg = i_sc
        matched_labels.append(labels[i_seg])

    if len(set(matched_labels)) != len(matched_labels) and smoothing == 'high':
        LOGGER.warning('Could not accurately map labels to segments. It will not '
                       'be possible to run fsw_file_writer.rewrite_prc using the '
                       'all_found_psfs*.txt file generated here.')

    return matched_labels


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
                LOGGER.error(err_message)

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
            LOGGER.error(err_message)

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
        LOGGER.error(err_message)

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


def copy_psfs_files(guiding_selections_file, output_file, root, guider, out_dir):
    """By parsing the name of an input guiding_selections*.txt file,
    identify a corresponding all_found_psfs*.txt file.

    Parameters
    ----------
    guiding_selections_file : str
        Path to unshifted_guiding_selections*.txt file
    output_file : str
        The type of file to copy and return the path of.
        Eg 'all_found_psfs' or 'psf_center'
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
    filename = os.path.basename(guiding_selections_file)
    if 'guiding_selections_' in filename:
        if 'unshifted' in filename:
            imported_root = filename.split('unshifted_guiding_selections_')[-1].split('.txt')[0]
        else:
            imported_root = filename.split('guiding_selections_')[-1].split('.txt')[0]
    elif '_regfile' in filename:
        imported_root = filename.split('_regfile.txt')[0]
    else:
        imported_root = None
        LOGGER.warning(
            'Could not parse root from provided guiding selections file ({}). '.format(filename) +
            'Not able to copy over a corresponding {}*.txt file.'.format(output_file)
        )

    # Try to copy the corresponding *all_found_psfs*.txt or *psf_center*.txt file, if present.
    if imported_root is not None:
        dir_to_look = os.path.dirname(guiding_selections_file)
        psfs_file = os.path.join(dir_to_look,
                                           'unshifted_{}_{}.txt'.format(output_file, imported_root))
        psfs_file_old = os.path.join(dir_to_look,
                                           '{}_{}.txt'.format(output_file, imported_root))
        if output_file == 'all_found_psfs':
            all_found_psfs_file_old2 = os.path.join(dir_to_look,
                                                   '{}_ALLpsfs.txt'.format(imported_root))
        file_to_copy = None

        if os.path.exists(psfs_file):
            file_to_copy = psfs_file
        elif os.path.exists(psfs_file_old):
            file_to_copy = psfs_file_old
        elif os.path.exists(all_found_psfs_file_old2) and output_file == 'all_found_psfs':
            file_to_copy = all_found_psfs_file_old2
        else:
            LOGGER.warning(
                'Could not find a corresponding *{}*.txt file for '
                'the provided guiding selections file ({}) in {}'.format(output_file, filename, dir_to_look)
            )

        if file_to_copy is not None:
            copied_psfs_file = os.path.join(out_dir, 'unshifted_{}_{}_G{}.txt'.format(output_file, root, guider))
            shutil.copy(file_to_copy, copied_psfs_file)
            LOGGER.info('Copying over {}'.format(file_to_copy))
            LOGGER.info('Successfully wrote: {}'.format(copied_psfs_file))

            return copied_psfs_file


def manual_star_selection(data, smoothing, choose_center=False, testing=False, masterGUIapp=None):
    """Launches a GUI to prompt the user to click-to-select guide and
    reference stars.

    Algorithmically find and locate all PSFs in image using
    photutils.find_peaks; prompt user to select guide and reference stars
    using the GUI.

    Parameters
    ----------
    data : 2-D numpy array
        Image data
    smoothing: str, optional
        Options are "low" for minimal smoothing (e.g. MIMF), "high" for large
        smoothing (e.g. GA), or "default" for medium smoothing for other cases
    choose_center : bool
        Automatically choose the one, highly-smoothed, PSF found in the image
    testing : bool, optional
        Generates indices randomly (for running pytests)
    masterGUIapp : qApplication, optional
        qApplication instance of parent GUI

    Returns
    -------
    cols : list
        List of positions and countrates of selected segments
    coords : list
        List of tuples with X and Y positions of selected PSFs
    nref : int
        The number of selected reference stars
    all_cols : list
        List of positions and countrates of all segments in the image

    Raises
    ------
    ValueError
        The user closed the GUI without selecting any stars.
    """
    if smoothing == 'high':
        gauss_sigma = 26
        npeaks = np.inf
    elif smoothing == 'low':
        gauss_sigma = 1
        npeaks = 1
    elif choose_center:
        gauss_sigma = 26
        npeaks = 1
    elif smoothing == 'default':
        gauss_sigma = 5
        npeaks = np.inf

    data = data.astype(float)

    smoothed_data = ndimage.gaussian_filter(data, sigma=gauss_sigma)

    # Use photutils.find_peaks to locate all PSFs in image
    num_psfs, coords, threshold = count_psfs(smoothed_data, gauss_sigma, npeaks=npeaks,
                                             choose=False)
    x, y = map(list, zip(*coords))

    # Use labeling to map locations of objects in array
    # (Kept for possible alternate countrate calculations; see count_rate_total)
    objects = ndimage.measurements.label(smoothed_data > threshold)[0]
    # NOTE: num_objects might not equal num_psfs

    # Find the minimum distance between PSFs
    if len(coords) < 2:
        # For cases where we only have star, we assume that we are sufficiently
        # isolated from other stars, but also that the guide star's PSF may be
        # distorted enough that it might appear quite large on the detector
        dist = 20
    else:
        dist = np.floor(np.min(utils.find_dist_between_points(coords))) - 1.

    # Calculate count rate
    countrate, val = count_rate_total(data, objects, num_psfs, x, y, countrate_3x3=True)

    # Call the GUI to pick PSF indices
    if not testing and not choose_center:
        gui_data = data.copy()
        gui_data[data == 0] = 1  # Alter null pixel values for LogNorm imshow
        inds = SelectStarsGUI.run_SelectStars(gui_data, x, y, dist,
                                              print_output=False,
                                              masterGUIapp=masterGUIapp)
    # Skip the GUI and choose the 0th PSF found (should only use this case when you'll only find 1 PSF)
    elif choose_center:
        inds = [0]
    # If in testing mode, just make a random list of indices
    else:
        # Make random list of inds
        n_select = min(11, num_psfs)
        LOGGER.info('Star Selection: Testing mode; selecting {} PSFs at random.'.format(n_select))
        inds = random.sample(range(num_psfs), n_select)

    nref = len(inds) - 1
    if len(inds) == 0:
        raise ValueError('Star Selection: No guide star and no reference stars selected')
    else:
        LOGGER.info('Star Selection: 1 guide star and {} reference stars selected'.format(nref))

    segment_labels = match_psfs_to_segments(x, y, smoothing)
    all_cols = create_cols_for_coords_counts(x, y, countrate, val,
                                             labels=segment_labels,
                                             inds=range(len(x)))
    cols = create_cols_for_coords_counts(x, y, countrate, val, inds=inds)

    return cols, coords, nref, all_cols


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# MAIN FUNCTION
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


def select_psfs(data, root, guider, guiding_selections_file=None,
                smoothing='default', choose_center=False,
                testing=False, out_dir=None, masterGUIapp=None, logger_passed=False):
    """Select guide and reference segments.

    Locate all of the segments in the provided data, then either parse a
    catalog of sources or prompt the user with a GUI to select which
    segments as the guide and reference segments. Generates two file:
        unshifted_all_found_psfs_{root}_G{guider}.txt
            Lists the locations and count rates of all segments found
            in the data
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
    guiding_selections_file : str, optional
        File containing locations and count rates of selected segments
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
        jwst_magic/
    masterGUIapp : qApplication, optional
        qApplication instance of parent GUI
    logger_passed : bool, optional
        Denotes if a logger object has already been generated.

    Returns
    -------
    guiding_selections_path : str
        Path to the unshifted_guiding_selections_{root}_G{guider}.txt file, which
        contains the locations and count rates of the selected guide and
        reference segments
    all_found_psfs_path : str
        Path to the unshifted_all_found_psfs_{root}_G{guider}.txt file, which
        the locations and count rates of all segments found in the data
    """
    if not logger_passed:
        utils.create_logger_from_yaml(__name__, root=root, level='DEBUG')

    try:
        out_dir = utils.make_out_dir(out_dir, OUT_PATH, root)
        utils.ensure_dir_exists(out_dir)

        # Read in image (check if it is a filename)
        if isinstance(data, str):
            data = fits.getdata(data)

        if guiding_selections_file:
            # Determine the kind of in_file and parse out the PSF locations and
            # countrates accordingly
            LOGGER.info(
                "Star Selection: Reading guide and reference star positions from {}"
                    .format(guiding_selections_file)
            )
            cols, _, _ = parse_in_file(guiding_selections_file)

            # Copy over corresponding all_found_psfs and psf_center file, if possible.
            all_cols = None
            try:
                all_found_psfs_path = copy_psfs_files(guiding_selections_file, 'all_found_psfs',root, guider, out_dir)
            except shutil.SameFileError:
                all_found_psfs_path = guiding_selections_file.replace('guiding_selections', 'all_found_psfs')
            _, all_coords, _= parse_in_file(all_found_psfs_path)

            psf_center_path = None
            try:
                psf_center_path = copy_psfs_files(guiding_selections_file, 'psf_center', root, guider, out_dir)
            except:
                pass

        else:
            # If no .incat or reg file provided, create reg file with manual
            # star selection using the SelectStarsGUI
            cols, all_coords, nref, all_cols = manual_star_selection(data,
                                                                 smoothing,
                                                                 choose_center,
                                                                 testing,
                                                                 masterGUIapp)
            all_found_psfs_path = None
            psf_center_path = None

        # Save PNG of image and all PSF locations in out_dir
        if not JENKINS:
            plot_centroids(data, all_coords, root, guider, out_dir)  # coords are in (x,y)

        if all_cols:
            all_found_psfs_path = os.path.join(out_dir, 'unshifted_all_found_psfs_{}_G{}.txt'.format(root, guider))
            # Write catalog of all identified PSFs
            utils.write_cols_to_file(all_found_psfs_path,
                                     labels=['label', 'y', 'x', 'countrate'],
                                     cols=all_cols, log=LOGGER)

        # Write catalog of selected PSFs
        guiding_selections_path = os.path.join(out_dir, 'unshifted_guiding_selections_{}_G{}.txt'.format(root, guider))
        utils.write_cols_to_file(guiding_selections_path,
                                 labels=['y', 'x', 'countrate'],
                                 cols=cols, log=LOGGER)

        # Calculate and write out center of PSF information for trk file if smoothing is low
        if smoothing == 'low' and psf_center_path is None:
            LOGGER.info(
                "Star Selection: No smoothing chosen so re-running star selection to also calculate PSF center")
            cols_center, _, _, _ = manual_star_selection(data,
                                                         smoothing='default',
                                                         choose_center=True,
                                                         testing=testing,
                                                         masterGUIapp=masterGUIapp)

            LOGGER.info(
                "Star Selection: PSF center information {} vs Guiding knot information {}".format(
                    cols_center, cols))
            psf_center_path = os.path.join(out_dir, 'unshifted_psf_center_{}_G{}.txt'.format(root, guider))
            utils.write_cols_to_file(psf_center_path,
                                     labels=['y', 'x', 'countrate'],
                                     cols=cols_center, log=LOGGER)

    except Exception as e:
        LOGGER.exception(e)
        raise

    return guiding_selections_path, all_found_psfs_path, psf_center_path
