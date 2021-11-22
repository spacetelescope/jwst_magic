import os

from astropy.io import fits
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from matplotlib.colors import from_levels_and_colors
from matplotlib import cm
import numpy as np
import pandas as pd
from scipy import ndimage

from jwst_magic import convert_image


# DO NOT CHANGE THESE LOCATIONS
GA_PSF_LOCATIONS = {'A1-1':[1037, 599], 'A2-2':[1391, 808], 'A3-3':[1375, 1218],
                    'A4-4':[1020, 1426], 'A5-5':[670, 1204], 'A6-6':[660, 783],
                    'B1-7':[1047, 177], 'B2-9':[1756, 608], 'B3-11':[1771, 1459],
                    'B4-13':[1010, 1789], 'B5-15':[297, 1409], 'B6-17':[327, 589],
                    'C1-8':[1394, 404], 'C2-10':[1748, 1024], 'C3-12':[1370, 1637],
                    'C4-14':[650, 1630], 'C5-16':[309, 987], 'C6-18':[690, 391]
                    }


def correct_pixel_values(image):
    """
    Make sure there are no 0/negative numbers in the image for plotting purposes, and make
    sure all hot pixels (over 100000) are reset (they confuse the PSF finding algorithm)
    """
    image_copy = image.copy()
    image_copy[image_copy > 100000] = 0
    image_copy[image_copy <= 0] = 3e-20

    return image_copy


def distance_to_target_center(x, y, truth_x, truth_y):
    '''
    calculate the distance from the PSF in OTE-01 to it's location in OTE-06/LOS-02
    truth_x and truth_y are the location of the boresight
    '''

    return np.sqrt((truth_x - x)**2 + (truth_y - y)**2)


def get_position_from_magic(image, smoothing='default', npeaks=np.inf):
    """ Get x and y postitions of each identifiable PSF in the MOSAIC. The identified PSFs
    might not all be from the target star and some from the target star might be missing
    """

    if smoothing == 'high':
        gauss_sigma = 26
    elif smoothing == 'default':
        gauss_sigma = 5
    else:
        raise ValueError("Only 'high' and 'default' are supported.")

    data = image.copy()
    data = data.astype(float)
    smoothed_data = ndimage.gaussian_filter(data, sigma=gauss_sigma)

    # Use photutils.find_peaks to locate all PSFs in image
    _, coords, _ = convert_image.convert_image_to_raw_fgs.count_psfs(smoothed_data,
                                                                     gauss_sigma,
                                                                     npeaks=npeaks,
                                                                     choose=False)
    x_list, y_list = map(list, zip(*coords))

    return x_list, y_list


def get_image_information(image, info_from_shadow=None, params_to_include=None,
                          target_location=None, smoothing='high'):
    """
    Take in the pick_log.txt file from Shadow to load in PSF characteristics. If no pick_log.txt
    is provided, the locations of the PSFs are found with a MAGIC function and no other PSF
    characteristics are measured.
    """
    if info_from_shadow is not None:
        print(f'Creating the information table from {info_from_shadow}')
        info_df = read_shadow_log(info_from_shadow, params_to_include=params_to_include)
    else:
        x_list, y_list = get_position_from_magic(image, smoothing=smoothing)
        print(f'Creating the information table from the X and Y values calculated using MAGIC')
        print(f'{len(x_list)} PSFs found')
        if x_list:
            info_df = pd.DataFrame(data={'x': x_list, 'y': y_list})
        else:
            info_df = None

    # Make a larger dictionary that organizes the information by segment
    if target_location is None:
        target_x = np.median(info_df['x'].values)
        target_y = np.median(info_df['y'].values)
        target_location = (int(target_x), int(target_y))
    else:
        target_x, target_y = target_location
    info_df['distance_to_target'] = distance_to_target_center(info_df['x'].values,
                                                              info_df['y'].values,
                                                              target_x,
                                                              target_y)

    return info_df, target_location


def convert_df_to_dictionary(info_df):
    """
    """
    if info_df is not None:
        info_dictionary = info_df.to_dict(orient='index')
    else:
        info_dictionary = None

    return info_dictionary


def read_shadow_log(filename, params_to_include=['x', 'y', 'fwhm', 'fwhm_x', 'fwhm_y', 'ellipse']):
    """
    From the txt file provided from Shadow, pull out helpful information.

    The default is to pull out the x & y locations of each identified PSF, the average FWHM,
    the FWHM in x and y, and the ellipticity.
    """
    labels = ['ra_txt', 'dec_txt', 'equinox', 'x', 'y', 'fwhm', 'fwhm_x', 'fwhm_y', 'starsize',
              'ellipse', 'background', 'skylevel', 'brightness', 'time_local', 'time_ut', 'ra_deg',
              'dec_deg']

    info = pd.read_csv(filename, delimiter=' ', skiprows=[0, 1], names=labels)

    info_squeezed = info[params_to_include].copy()
    return info_squeezed


def match_psf_params_to_segment(info_dictionary, matching_dictionary):
    """
    Match the location and other PSF parameters to the segment based on the matching
    dictionary. Also calculate the distance from each segment to the boresight
    """

    seg_location_dictionary = {}
    seg_location_dictionary = {matching_dictionary[ind]:info_dictionary[ind] \
                                    for ind in matching_dictionary.keys()}

    return seg_location_dictionary


def add_segment_to_dict(info_dictionary, matching_dictionary):
    """
    Match the location and other PSF parameters to the segment based on the matching
    dictionary. Also calculate the distance from each segment to the boresight
    """
    for k in info_dictionary.keys():
        try:
            info_dictionary[k]['seg_id'] = matching_dictionary[k]
        except KeyError:
            info_dictionary[k]['seg_id'] = None

    return info_dictionary


def grab_vmin_vmax(param):
    """ Give a vmin & vmax depending on what parameter we are plotting
    """
    if 'fwhm' in param:
        vmin = 10 # TBD
        vmax = 20 # TBD
        cmap = cm.RdYlGn_r
    elif param == 'ellipse':
        vmin = 0
        vmax = 1
        cmap = cm.RdYlGn
    else:
        vmin = None
        vmax = None
        cmap = cm.RdYlGn_r

    return vmin, vmax, cmap


def plot_mosaic_with_psfs(mosaic, car, x_list=None, y_list=None, est_target_location=None,
                          xlim=None, ylim=None, label_color='white'):
    """
    Plot out the mosaic image with the identified PSFs and estimated target location plotted
    """
    plt.figure(figsize=(14, 14))

    plt.imshow(mosaic, norm=LogNorm(vmin=1, vmax=1000), origin='lower')
    if est_target_location is not None:
        plt.scatter(est_target_location[0], est_target_location[1], s=500, marker='*', c='white',
                    edgecolor='C1', label='Estimated Target Location')
        plt.legend()
    if x_list is not None and y_list is not None:
        for j, (x, y) in enumerate(zip(x_list, y_list)):
            plt.annotate(j, xy=(x, y), xytext=(x-75, y), c=label_color, fontsize=14, weight='bold')
        if xlim is not None:
            plt.xlim(xlim)
        if ylim is not None:
            plt.ylim(ylim)

    plt.title(f"{car} Mosaic")
    plt.show()


def plot_each_identified_segment_psf(image, x_list, y_list, window_size=200, num_psf_per_row=6,
                                     seg_list=None):
    """
    Cut each segment out of the mosaic/image provided based on the x and y locations given with
    x_list, and y_list
    """
    radius = window_size//2

    rows = len(x_list)//num_psf_per_row+1
    columns = len(x_list) if len(x_list) < num_psf_per_row else num_psf_per_row
    fig, ax = plt.subplots(rows, columns, figsize=(6*columns, 6*rows))

    i = 0
    for j in range(rows):
        for k in range(columns):
            try:
                cutout = image[int(y_list[i])-radius:int(y_list[i])+radius,
                               int(x_list[i])-radius:int(x_list[i])+radius]
                ax[j, k].imshow(cutout, origin='lower', norm=LogNorm(vmin=1))
                ind = seg_list[i] if seg_list else i
                ax[j, k].set_title(f'PSF {ind}')
                i += 1
            except IndexError:
                ax[j, k].axis('off')

    plt.show()


#TODO: remove?
def plot_parameters(xs, ys, values, title, cmap=cm.RdYlGn_r):
    """
    Plot the PSF parameters
    """
    center_x = np.median(xs)
    center_y = np.median(ys)

    # plot out new image with distance color coded by how far from boresight
    plt.figure(figsize=(10, 10))
    plt.title(title)
    plt.scatter(xs, ys, marker='o', s=100, c=values, cmap=cmap)
    plt.colorbar()
    plt.ylim(center_y+1024, center_y-1024)
    plt.xlim(center_x-1024, center_x+1024)
    plt.show()


def plot_multiple_parameters(parameters, info_df, xs=None, ys=None, xlim=None, ylim=None):
    """
    Give a list of parameters that you want to plot information for. This plot will show the
    value of each at the location of the PSF that is provided in the info_dictionary, but if the
    parameters xs and ys are provided, those will be used instead. When we are plotting the PSFs
    in their GA location, we must provide the GA xs and ys in order to overide the position of
    each PSF in the original image.

    """
    # Check that all parameters are in info_dictionary
    parameters = list(set(parameters) & set(list(info_df.columns)))
    if not parameters:
        print('Requested parameters are not in the povided data frame. Nothing to plot')
    else:
        fig, ax = plt.subplots(1, len(parameters), figsize=(6*len(parameters), 8))
        try:
            segs = info_df['seg_id'].values
        except KeyError:
            segs = np.arange(len(info_df))
        if xs is None and ys is None:
            xs = info_df['x'].values
            ys = info_df['y'].values
        for i, param in enumerate(parameters):
            values = info_df[param].values
            for seg, x, y in zip(segs, xs, ys):
                try:
                    ax[i].annotate(seg, (x, y), (x+15, y+15))
                except TypeError:
                    ax.annotate(seg, (x, y), (x+15, y+15))
            vmin, vmax, cmap = grab_vmin_vmax(param)
            try:
                pl = ax[i].scatter(xs, ys, marker='o', s=100, c=values, cmap=cmap,
                                   vmin=vmin, vmax=vmax)
                fig.colorbar(pl, ax=ax[i], orientation='horizontal')
                ax[i].set_title(f"{param} for each PSF in the image")
            except TypeError:
                pl = ax.scatter(xs, ys, marker='o', s=100, c=values, cmap=cmap,
                                vmin=vmin, vmax=vmax)
                fig.colorbar(pl, ax=ax, orientation='horizontal')
                ax.set_title(f"{param} for each PSF in the image")
            if xlim and ylim:
                plt.xlim(xlim)
                plt.ylim(ylim)


def get_nrc_data_from_list(nrc_file_list):
    """
    From a list of nircam A or B images, create a list of the data and filenames
    """
    data_list = []
    name_list = []
    for det_file in nrc_file_list:
        data = fits.getdata(det_file)
        data_list.append(correct_pixel_values(data))
        name_list.append(det_file.split('/')[-1])

    return data_list, name_list


def plot_nrca_images(data_list, name_list):
    """
    Plot the 4 SW NRCA images in their respecitve postions
    """

    a1, a2, a3, a4 = data_list
    a1_name, a2_name, a3_name, a4_name = name_list

    fig, ax = plt.subplots(2, 2, figsize=(12, 12))
    ax[0, 0].imshow(a2, norm=LogNorm(vmin=1, vmax=1000), origin='lower')
    ax[0, 0].set_title(f'A2: {a2_name}')
    ax[0, 1].imshow(a4, norm=LogNorm(vmin=1, vmax=1000), origin='lower')
    ax[0, 1].set_title(f'A4: {a4_name}')
    ax[1, 0].imshow(a1, norm=LogNorm(vmin=1, vmax=1000), origin='lower')
    ax[1, 0].set_title(f'A1: {a1_name}')
    ax[1, 1].imshow(a3, norm=LogNorm(vmin=1, vmax=1000), origin='lower')
    ax[1, 1].set_title(f'A3: {a3_name}')

    plt.show()


def plot_nrcb_images(data_list, name_list):
    """
    Plot the 4 SW NRCB images in their respecitve postions
    """

    b1, b2, b3, b4 = data_list
    b1_name, b2_name, b3_name, b4_name = name_list

    fig, ax = plt.subplots(2, 2, figsize=(12, 12))
    ax[0, 0].imshow(b3, norm=LogNorm(vmin=1, vmax=1000), origin='lower')
    ax[0, 0].set_title(f'B3: {b3_name}')
    ax[0, 1].imshow(b1, norm=LogNorm(vmin=1, vmax=1000), origin='lower')
    ax[0, 1].set_title(f'B1: {b1_name}')
    ax[1, 0].imshow(b4, norm=LogNorm(vmin=1, vmax=1000), origin='lower')
    ax[1, 0].set_title(f'B4: {b4_name}')
    ax[1, 1].imshow(b2, norm=LogNorm(vmin=1, vmax=1000), origin='lower')
    ax[1, 1].set_title(f'B2: {b2_name}')

    plt.show()


def create_basic_mosaic(nrca_data_list, nrcb_data_list):
    """
    Create a mosaic array from 4 nrca and 4 nrcb images.
    """
    a1, a2, a3, a4 = nrca_data_list
    b1, b2, b3, b4 = nrcb_data_list

    big_image = np.zeros((2*2048, 4*2048))

    # NRCA
    big_image[:2048, :2048] = a1
    big_image[:2048, 2048:4096] = a3

    big_image[2048:4096, :2048] = a2
    big_image[2048:4096, 2048:4096] = a4

    #NRCB
    big_image[:2048, 4096:6144] = b4
    big_image[:2048, 6144:8192] = b2

    big_image[2048:4096, 4096:6144] = b3
    big_image[2048:4096, 6144:8192] = b1

    return big_image


def create_image_array(image, seg_location_dictionary, window_size=60):
    """
    Cut out postage stamps of an image and using the locations of each PSF and it's assumed
    segment ID, place it in a large image array.
    """
    # Create an array of zeros the size of a standard dector size
    large_image_array = np.zeros([2048, 2048])

    ga_xs = []
    ga_ys = []
    for seg in seg_location_dictionary.keys():
        mosaic_x = seg_location_dictionary[seg]['x']
        mosaic_y = seg_location_dictionary[seg]['y']
        psf = image[mosaic_y-window_size:mosaic_y+window_size,
                    mosaic_x-window_size:mosaic_x+window_size]
        ga_x, ga_y = GA_PSF_LOCATIONS[seg]
        ga_xs.append(ga_x)
        ga_ys.append(ga_y)
        large_image_array[ga_y-window_size: ga_y+window_size,
                          ga_x-window_size:ga_x+window_size] = psf

    large_image_array = correct_pixel_values(large_image_array)

    return large_image_array, ga_xs, ga_ys


def list_good_psfs(info_dictionary, fwhm_limit, ellipse_limit):
    """
    """
    good_inds = []
    for key in info_dictionary.keys():
        if info_dictionary[key]['fwhm'] < fwhm_limit and \
           info_dictionary[key]['ellipse'] > ellipse_limit:
            good_inds.append(key)

    return good_inds
