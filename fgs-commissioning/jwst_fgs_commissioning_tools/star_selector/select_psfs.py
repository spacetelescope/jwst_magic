""" Find all the relevant PSFs in the image, be it manually or using a file."""
# STDLIB
import os
import sys
from inspect import currentframe, getframeinfo
import warnings
import string
import random

# Third Party
import matplotlib
if matplotlib.get_backend() != 'Qt5Agg':
    matplotlib.use('Qt5Agg')  # Make sure that we are using Qt5
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from matplotlib import rcParams

from astropy.io import fits
from astropy.io import ascii as asc
from astropy.stats import sigma_clipped_stats
import numpy as np
from photutils import find_peaks
from scipy import ndimage, signal


# LOCAL
from jwst_fgs_commissioning_tools import log, utils
from jwst_fgs_commissioning_tools.star_selector import SelectStarsGUI

# Adjust matplotlib origin
rcParams['image.origin'] = 'upper'

# Make plots pretty
rcParams['font.family'] = 'serif'
rcParams['font.weight']='light'
rcParams['mathtext.bf'] = 'serif:normal'

SS_PATH = os.path.dirname(os.path.realpath(__file__))
PACKAGE_PATH = os.path.split(SS_PATH)[0]
OUT_PATH = os.path.split(PACKAGE_PATH)[0]  # Location of out/ and logs/ directory

def count_psfs(smoothed_data, gauss_sigma, choose=False):
    """Use photutils.find_peaks to count how many PSFS are present in the data
    """

    if choose:
        num_psfs, coords, threshold = choose_threshold(smoothed_data, gauss_sigma)

    else:
        # Perform statistics
        mean, median, std = sigma_clipped_stats(smoothed_data, sigma=0, iters=0)

        # Find PSFs
        threshold = median + (3 * std)  # Used to be median + 3 * std

        sources = find_peaks(smoothed_data, threshold, box_size=gauss_sigma)
        num_psfs = len(sources)
        coords = sources['x_peak', 'y_peak']
        coords = [(x, y) for [x, y] in coords]

        log.info('Star Selection: {} PSFs detected in Gaussian-smoothed data \
            (threshold = {}; sigma = {})'.format(num_psfs, threshold, gauss_sigma))

    return num_psfs, coords, threshold


def choose_threshold(smoothed_data, gauss_sigma):
    # Perform statistics
    mean, median, std = sigma_clipped_stats(smoothed_data, sigma=0, iters=0)

    # Run find_peaks with two different threshold options
    thresholds = [3 * std, mean]
    sources_std = find_peaks(smoothed_data, thresholds[0], box_size=gauss_sigma)
    sources_mean = find_peaks(smoothed_data, thresholds[1], box_size=gauss_sigma)

    # Show plots of each for user to chooose between
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
    choice = raw_input('''
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
        log.error('Star Selection: User rejection of identified PSFs.')
        raise StandardError('User rejection of identified PSFs.')

def plot_centroids(data, coords, root, guider, out_dir):
    pad = 300

    # Determine x and y limits that encompass all PSFS
    xarray, yarray = [x for (x, y) in coords], [y for (x, y) in coords]
    x_mid = (min(xarray) + max(xarray)) / 2
    y_mid = (min(yarray) + max(yarray)) / 2
    x_range = max(xarray) - min(xarray)
    y_range = max(yarray) - min(yarray)
    ax_range = max(x_range, y_range) # Choose the larger of the dimensions
    ax_range += 100 # Make sure not to clip off the edge of border PSFS

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
    plt.savefig(os.path.join(out_dir, '{}_G{}_centers.png'.format(root, guider)))

    plt.close()


def count_rate_total(data, objects, num_objects, x, y, counts_3x3=True):
    """
    Get the x,y, and counts for each psf in the image

    The threshold default of 150, assumes that your data type is uint16. If not,
    a good start is 40.
    """

    counts = []
    val = []
    for i in range(1, num_objects + 1):
        im = np.copy(objects)
        im[objects != i] = False
        im[objects == i] = True

        if counts_3x3:
            counts.append(utils.countrate_3x3(x[i - 1], y[i - 1], np.array(data)))
        else:
            counts.append(np.sum(im * data))
        val.append(np.sum(im * 1.))  # Number of pixels in object

    return counts, val


def create_cols_for_coords_counts(x, y, counts, val, labels=None, inds=None):
    """
    Create an array of columns of y, x, and counts of each PSF to be written out.
    Use the inds returned from pick_stars based on user input.

    If no inds are given, put the PSF with the most compact PSF first in the list to make it the
    Guide Star. **This method is not fool-proof, use at own risk***
    """
    if labels:
        cols = [[ll, '{:.4f}'.format(yy),
                 '{:.4f}'.format(xx), '{:.4f}'.format(co)] for ll, yy, xx, co in zip(labels, y, x, counts)]# these coordinates are y,x
    else:
        cols = [[yy, xx, co] for yy, xx, co in zip(y, x, counts)]# these coordinates are y,x

    if inds is None:
        min_ind = np.where(val == np.min(val))[0][0]  # Find most compact PSF
        cols.insert(0, cols.pop(min_ind))  # Move most compact PSF to top of the list
    else:
        cols = [cols[i] for i in inds]

    return cols

def match_psfs_to_segments(x, y):
    labels = string.ascii_uppercase[:18]
    x_list = [3, 2, 4, 1, 3, 5, 2, 4, 1, 5, 2, 4, 1, 3, 5, 2, 4, 3]
    y_list = [9, 8, 8, 7, 7, 7, 6, 6, 5, 5, 4, 4, 3, 3, 3, 2, 2, 1]

    # Determine boundaries of array
    x_min = min(x)
    x_max = max(y)
    y_min = min(y)
    y_max = max(y)

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

    return matched_labels

def parse_in_file(in_file):
    '''Determines if the input file contains x, y, and countrate data. If so,
    extracts the locations and countrates of the stars accordingly. Recognizes
    columns named "x" or "xreal"; "y" or "yreal"; "countrate", "count rate", or
    "ctot".'''

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
                log.error(err_message)
                return

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
            log.error(err_message)
            return

    # Make sure all the necessary columns are present
    x_check = 'x' in colnames or 'xreal' in colnames
    y_check = 'y' in colnames or 'yreal' in colnames
    counts_check = 'countrate' in colnames or 'ctot' in colnames

    if not (x_check and y_check and counts_check):
        err_message = 'Unknown in_file format: {}. Expecting columns named \
                       "x"/"xreal", "y"/"yreal", and \
                       "count rate"/"countrate"/"ctot". Found columns \
                       named {}. Please rename columns.'.format(in_file, colnames)
        raise TypeError(err_message)
        log.error(err_message)
        return

    # Passed all the checkpoints! Move on to process the file.
    log.info('Star Selection: Selecting stars from input file {}'.format(in_file))

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
    coords = [(y, x) for x, y in zip(in_table['x'], in_table['y'])]
    nref = len(in_table['x']) - 1

    return out_cols, coords, nref

def manual_star_selection(data, global_alignment, testing=False):
    '''Algorithmically find and locate all PSFs in image using
    photutils.find_peaks; prompt user to select guide and reference stars
    using GUI.'''
    if global_alignment:
        gauss_sigma = 26
    else:
        gauss_sigma = 5

    smoothed_data = ndimage.gaussian_filter(data, sigma=gauss_sigma)

    # Use photutils.find_peaks to locate all PSFs in image
    num_psfs, coords, threshold = count_psfs(smoothed_data, gauss_sigma,
                                             choose=False)
    x, y = map(list, zip(*coords))

    # Use labeling to map locations of objects in array
    # (Kept for possible alternate countrate calculations; see count_rate_total)
    objects = ndimage.measurements.label(smoothed_data > threshold)[0]
    # NOTE: num_objects might not equal num_psfs

    # Find the minimum distance between PSFs
    if len(coords) < 2:
        # For cases where we only have star, we assume that we are sufficiently
        #isolated from other stars, but also that the guide star's PSF may be
        #distorted enough that it might appear quite large on the detector
        dist = 20
    else:
        dist = np.floor(np.min(utils.find_dist_between_points(coords))) - 1.

    # Calculate count rate
    counts, val = count_rate_total(data, objects, num_psfs, x, y, counts_3x3=True)

    # Call the GUI to pick PSF indices
    if not testing:
        gui_data = data.copy()
        gui_data[data == 0] = 0.1  # Alter null pixel values for LogNorm imshow
        inds = SelectStarsGUI.run_SelectStars(gui_data, x, y, dist,
                                              print_output=False)
    # If in testing mode, just make a random list of indices
    else:
        # Make random list of inds
        n_select = min(11, num_psfs)
        log.info('Testing mode; selecting {} PSFs at random.'.format(n_select))
        inds = random.sample(range(num_psfs), n_select)

    nref = len(inds) - 1
    if len(inds) == 0:
        log.info('Star Selection: No guide star and no reference stars selected')
    else:
        log.info('Star Selection: 1 guide star and {} reference stars selected'.format(nref))

    segment_labels = match_psfs_to_segments(x, y)
    ALL_cols = create_cols_for_coords_counts(x, y, counts, val,
                                             labels=segment_labels,
                                             inds=range(len(x)))
    cols = create_cols_for_coords_counts(x, y, counts, val, inds=inds)

    return cols, coords, nref, ALL_cols


def create_reg_file(data, root, guider, in_file=None,
                    global_alignment=False, return_nref=False, testing=False,
                    out_dir=None):

    out_dir = utils.make_out_dir(out_dir, OUT_PATH, root)
    utils.ensure_dir_exists(out_dir)

    # Any value above 65535 or below 0 will wrap when converted to uint16
    data = utils.correct_image(data, upper_threshold=65535, upper_limit=65535)

    if in_file:
        # Determine the kind of in_file and parse out the PSF locations and
        # countrates accordingly
        cols, coords, nref = parse_in_file(in_file)
        ALL_cols = None

    else:
        # If no .incat or reg file provided, create reg file with manual
        # star selection using the SelectStarsGUI
        cols, coords, nref, ALL_cols = manual_star_selection(data, global_alignment, testing)

    # Save PNG of image and all PSF locations in out_dir
    plot_centroids(data, coords, root, guider, out_dir)

    if ALL_cols:
        # Write out file of ALL identified PSFs
        utils.write_cols_to_file(out_dir,
                                 filename='{0}_G{1}_ALLpsfs.txt'.format(root, guider),
                                 labels=['label', 'y', 'x', 'countrate'],
                                 cols=ALL_cols)


    # Write out regfile of selected PSFs
    utils.write_cols_to_file(out_dir,
                             filename='{0}_G{1}_regfile.txt'.format(root, guider),
                             labels=['y', 'x', 'countrate'],
                             cols=cols)

    if return_nref:
        return cols, nref
