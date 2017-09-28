# STDLIB
import os
import sys
from inspect import currentframe, getframeinfo
import warnings

# Third Party
from astropy.io import fits
from astropy.io import ascii as asc
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from matplotlib import rcParams
from scipy import ndimage, signal
from astropy.stats import sigma_clipped_stats
from photutils import find_peaks

# LOCAL
import utils
import log
import SelectStarsGUI

# Adjust matplotlib origin
rcParams['image.origin'] = 'upper'

# Make plots pretty
rcParams['font.family'] = 'serif'
rcParams['font.weight']='light'
rcParams['mathtext.bf'] = 'serif:normal'

# rcParams['xtick.labelsize'] = 14
# rcParams['ytick.labelsize'] = 14


def count_psfs(smoothed_data, gauss_sigma, choose=False):
    """Use photutils.find_peaks to count how many PSFS are present in the data
    """

    if choose:
        num_psfs, coords, threshold = choose_threshold(smoothed_data, gauss_sigma)

    else:
        # Perform statistics
        mean, median, std = sigma_clipped_stats(smoothed_data, sigma=0, iters=0)

        # Find PSFs
        threshold = (3 * std)  # Used to be median + 3 * std
        sources = find_peaks(smoothed_data, threshold, box_size=gauss_sigma)
        num_psfs = len(sources)
        coords = sources['x_peak', 'y_peak']
        coords = [(x, y) for [x, y] in coords]

        log.info('{} PSFs detected in Gaussian-smoothed data \
            (threshold = {}; sigma = {})'.format(num_psfs, threshold, gauss_sigma))

    return num_psfs, coords, threshold


def choose_threshold(smoothed_data, gauss_sigma):
    # Perform statistics
    mean, median, std = sigma_clipped_stats(smoothed_data, sigma=0, iters=0)

    # Run find_peaks with two different threshold options
    thresholds = [3*std, mean]
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
    ax2.set_title('Threshold = Mean ({} sources found)'.format(len(sources_mean)))

    # fig = plt.gcf()
    # fig.canvas.manager.window.raise_()
    plt.get_current_fig_manager().window.raise_()

    plt.show()

    # Prompt user to choose
    choice = raw_input('''
Examine the two options presented. To use the stars selected with a 3 \
standard deviation threshold, type "S". To use the stars selected with a mean \
threshold, type "M". To use neither and cancel the program, press enter.

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
        log.error('User rejection of identified PSFs.')
        raise StandardError('User rejection of identified PSFs.')


def countrate_3x3(coords, data):
    """
    Using the coordinates of each PSF,place a 3x3 box around center pixel, and sum
    the counts of the pixels in this box.
    """
    xx = int(coords[1])
    yy = int(coords[0])

    counts = np.sum(data[yy - 1:yy + 2, xx - 1:xx + 2])
    return counts


def plot_centroids(data, coords, root, guider, out_dir):
    # if compact:
    #     pad = 300
    # else:
    #     pad = 500
    pad=300

    # Determine x and y limits that encompass all PSFS
    xarray, yarray = [x for (x,y) in coords], [y for (x,y) in coords] # Backwards... :^(
    x_mid = (min(xarray) + max(xarray)) / 2
    y_mid = (min(yarray) + max(yarray)) / 2
    x_range = max(xarray) - min(xarray)
    y_range = max(yarray) - min(yarray)
    ax_range = max(x_range, y_range) # Choose the larger of the dimensions
    ax_range += 100 # Make sure not to clip off the edge of border PSFS

    plt.figure(figsize=(17,17))
    plt.imshow(data,cmap='Greys',norm=LogNorm())
    for i in range(len(coords)):
        plt.scatter(coords[i][1],coords[i][0])
        plt.annotate('({},{})'.format(int(coords[i][0]),int(coords[i][1])),
                     (coords[i][1]-(pad*0.05), coords[i][0]+(pad*0.05)))
    plt.title('Centroids found for {}'.format(root))
    plt.ylim(min(2048, x_mid + ax_range/2), max(0, x_mid - ax_range/2))
    plt.xlim(max(0, y_mid - ax_range/2), min(2048, y_mid + ax_range/2))
    plt.savefig(os.path.join(out_dir,'{}_G{}_centers.png'.format(root,guider)))

    plt.close()


def count_rate_total(data, objects, num_objects, gs_points, counts_3x3=True):
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
            counts.append(countrate_3x3(gs_points[i-1], data))
        else:
            counts.append(np.sum(im*data))
        val.append(np.sum(im*1.))

    return counts, val


def create_cols_for_coords_counts(y, x, counts, val, inds=None):
    """
    Create an array of columns of y, x, and counts of each PSF to be written out.
    Use the inds returned from pick_stars based on user input.

    If no inds are given, put the PSF with the most compact PSF first in the list to make it the
    Guide Star. **This method is not fool-proof, use at own risk***
    """

    cols = []
    for i, (yy,xx,co) in enumerate(zip(y,x,counts)):
        cols.append([yy,xx,co]) ## these coordinates are y,x

    if inds is None:
        min_ind = np.where(val == np.min(val))[0][0] #Find most compact PSF
        cols.insert(0, cols.pop(min_ind)) #Move most compact PSF to top of the list
    else:
        cols = [cols[i] for i in inds]

    return cols


def create_reg_file(data, root, guider, out_dir, return_nref=False,
                    global_alignment=False, incat=None):
   
    # If no .incat file provided, create reg file with manual star selection in GUI
    if incat == None:
        if isinstance(data, str):
            data = utils.read_fits(data)[1]

        if global_alignment:
            gauss_sigma = 26
        else:
            gauss_sigma = 5

        smoothed_data = ndimage.gaussian_filter(data, sigma = gauss_sigma)

        # Use photutils.find_peaks to locate all PSFs in image
        num_psfs, coords, threshold = count_psfs(smoothed_data, gauss_sigma, choose=True)
        x, y = map(list, zip(*coords))

        # Use labeling to map locations of objects in array
        objects, num_objects = ndimage.measurements.label(smoothed_data > threshold)
        # NOTE: num_objects might not equal num_psfs

        if len(coords)<2:
            log.error('Less than two objects have been found. Cannot proceed. Exiting')
            raise ValueError('cannot guide on < 2 objects')

        #find the minimum distance between PSFs
        dist = np.floor(np.min(utils.find_dist_between_points(coords))) - 1. 

        plot_centroids(data, coords, root, guider, out_dir) # Save pretty PNG in out dir
        counts, val = count_rate_total(data, objects, num_psfs, coords, counts_3x3=True) # Calculate count rate
        # print(y, x, counts,val)

        # Call the GUI
        dataToShow = data.copy()
        dataToShow[data == 0] = 0.1 # Alter null pixel values to appear as black in LogNorm image
        inds = SelectStarsGUI.run_SelectStars(dataToShow, x, y, dist, printOutput=False)
        nref = len(inds)-1
        log.info('1 guide star and {} reference stars selected'.format(nref))

        cols = create_cols_for_coords_counts(y, x, counts, val, inds=inds)
   
    # If .incat file provided, create reg file with provided information
    else:
        # Read in .incat file
        log.info('Selecting stars from .incat file {}'.format(incat))
        incat = asc.read(incat, names=['x', 'y', 'ctot', 'inimg', 'incat'])
        incat['x'] = incat['x'].astype(int) # Need to be an integer
        incat['y'] = incat['y'].astype(int)
        incat['ctot'] = incat['ctot'].astype(int)

        # Only use stars in the catalog
        incat = incat[incat['incat']==1]
        # Only use relevant columns
        cols = incat['y', 'x', 'ctot'] # Make sure to flip x and y!!
        inds = incat['x']
        nref = len(inds)-1

        coords = [(y,x) for x, y in zip(incat['x'], incat['y'])]

        plot_centroids(data, coords, root, guider, out_dir)  # Save pretty PNG

    utils.write_cols_to_file(out_dir,
                           filename='{0}_G{1}_regfile.txt'.format(root, guider),
                           labels=['y','x','count rate'],
                           cols=cols)

    if return_nref:
        return nref
