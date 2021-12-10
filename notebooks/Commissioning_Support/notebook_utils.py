import os

from astropy.io import fits
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from matplotlib.colors import from_levels_and_colors
from matplotlib import cm
import numpy as np
import pandas as pd
import poppy
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


def get_image_information(image, params_to_include=None,
                          target_location=None, smoothing='high'):
    """
    Take in the pick_log.txt file from Shadow to load in PSF characteristics. If no pick_log.txt
    is provided, the locations of the PSFs are found with a MAGIC function and no other PSF
    characteristics are measured.
    """
    x_list, y_list = get_position_from_magic(image, smoothing=smoothing)
    fwhm_x, fwhmy, ee = get_psf_characteristics(image, x_list, y_list)
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


def add_segment_to_df(info_df, matching_dictionary):
    """
    Add the Segment ID to the dataframe based on the matching in matching_dictionary
    """
    segments = []

    for i in range(len(info_df)):
        try:
            segments.append(matching_dictionary[i])
        except KeyError:
            segments.append('Unknown')

    info_df['segment'] = segments
    return info_df


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


def plot_mosaic_with_psfs(mosaic, car, info_df, est_target_location=None,
                          xlim=None, ylim=None, label_color='white'):
    """
    Plot out the mosaic image with the identified PSFs and estimated target location plotted
    """
    x_list = info_df['x']
    y_list = info_df['y']
    try:
        seg_list = info_df['segment']
    except KeyError:
        seg_list = None

    plt.figure(figsize=(14, 14))
    plt.imshow(mosaic, norm=LogNorm(vmin=1, vmax=1000), origin='lower')
    if est_target_location is not None:
        plt.scatter(est_target_location[0], est_target_location[1], s=500, marker='*', c='white',
                    edgecolor='C1', label='Estimated Target Location')
        plt.legend()
    if x_list is not None and y_list is not None:
        for j, (x, y) in enumerate(zip(x_list, y_list)):
            name = seg_list[j] if seg_list is not None else j
            plt.annotate(name, xy=(x, y), xytext=(x-75, y), c=label_color, fontsize=14,
                         weight='bold')
        if xlim is not None:
            plt.xlim(xlim)
        if ylim is not None:
            plt.ylim(ylim)

    plt.title(f"{car} Mosaic")
    plt.show()


def plot_each_identified_segment_psf(image, info_df, window_size=200, num_psf_per_row=6):
    """
    Cut each segment out of the mosaic/image provided based on the x and y locations given with
    x_list, and y_list
    """
    x_list = info_df['x'].values
    y_list = info_df['y'].values
    try:
        seg_list = info_df['segment'].values
    except KeyError:
        seg_list = np.arange(len(x_list))

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
                ax[j, k].set_title(f'PSF {seg_list[i]}')
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
            segs = info_df['segment'][info_df['segment'] != 'Unknown'].values
        except KeyError:
            segs = np.arange(len(info_df))

        if xs is None and ys is None:
            xs = info_df['x'].values
            ys = info_df['y'].values

        for i, param in enumerate(parameters):
            ax = ax if len(parameters) == 1 else ax[i]

            for seg, x, y in zip(segs, xs, ys):
                ax.annotate(seg, (x, y), (x+15, y+15))

            vmin, vmax, cmap = grab_vmin_vmax(param)

            try:
                values = info_df[param][info_df['segment'] != 'Unknown'].values
            except KeyError:
                values = info_df[param].values

            pl = ax.scatter(xs, ys, marker='o', s=100, c=values, cmap=cmap,
                            vmin=vmin, vmax=vmax)
            fig.colorbar(pl, ax=ax, orientation='horizontal')
            ax.set_title(f"{param} for each PSF in the image")

            if xlim and ylim:
                plt.xlim(xlim)
                plt.ylim(ylim)


def separate_nircam_images(all_images_list):
    nrca_images = []
    nrcb_images = []
    for fi in all_images_list:
        instrument = fi.split('/')[-1].split('_')[3]
        if 'a' in instrument:
            nrca_images.append(fi)
        elif 'b' in instrument:
            nrcb_images.append(fi)

    return sorted(nrca_images), sorted(nrcb_images)


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


def create_image_array(image, info_df, window_size=60):
    """
    Cut out postage stamps of an image and using the locations of each PSF and it's assumed
    segment ID, place it in a large image array.
    """
    # Create an array of zeros the size of a standard dector size
    large_image_array = np.zeros([2048, 2048])

    ga_xs = []
    ga_ys = []
    x_list = info_df['x'].values
    y_list = info_df['y'].values
    try:
        seg_list = info_df['segment'].values

        for mosaic_x, mosaic_y, seg in zip(x_list, y_list, seg_list):
            if seg != 'Unknown':
                psf = image[mosaic_y-window_size:mosaic_y+window_size,
                            mosaic_x-window_size:mosaic_x+window_size]
                ga_x, ga_y = GA_PSF_LOCATIONS[seg]
                ga_xs.append(ga_x)
                ga_ys.append(ga_y)
                large_image_array[ga_y-window_size: ga_y+window_size,
                                  ga_x-window_size:ga_x+window_size] = psf

        large_image_array = correct_pixel_values(large_image_array)
        return large_image_array, ga_xs, ga_ys

    except KeyError:
        print('No segment list. Cannot make pseudo large image array')
        return


def list_good_psfs(info_df, fwhm_limit, ellipse_limit, seg_list=None):
    """
    Determine which PSFs are "good" based on their FWHM and ellipticity values
    """
    good_inds = []
    for i, (fwhm, ellipse) in enumerate(zip(info_df['fwhm'], info_df['ellipse'])):
        if fwhm < fwhm_limit and ellipse > ellipse_limit:
            if seg_list:
                good_inds.append(seg_list[i])
            else:
                good_inds.append(i)

    return good_inds

# Got this from https://grit.stsci.edu/wfsc/tools/-/blob/master/ote-commissioning/pre-mimf/ote28_psf_analysis.py
def measure_fwhm(array):
    """Fit a Gaussian2D model to a PSF and return the fitted PSF
    the FWHM is x and y can be found with fitted_psf.x_fwhm, fitted_psf.y_fwhm

    Parameters
    ----------
    array : numpy.ndarray
        Array containing PSF

    Returns
    -------
    x_fwhm : float
        FWHM in x direction in units of pixels

    y_fwhm : float
        FWHM in y direction in units of pixels
    """
    yp, xp = array.shape
    y, x, = np.mgrid[:yp, :xp]
    p_init = models.Gaussian2D(amplitude = array.max(), x_mean=xp*0.5,y_mean=yp*0.5)
    fit_p = fitting.LevMarLSQFitter()
    fitted_psf = fit_p(p_init, x, y, array)
    return fitted_psf.x_fwhm, fitted_psf.y_fwhm

def measure_ee(array):
    '''Wrapper function around poppy's measure_ee
    '''
    hdu1 = fits.PrimaryHDU(array)
    new_hdul = fits.HDUList([hdu1])
    return poppy.measure_ee(new_hdul, normalize='None')

def get_psf_characteristics(data, x_list, y_list, radius):
        '''

        '''
        # results = peaks.evaluate_peaks(psf_peaks, data, fwhm_radius=fwhm_radius,
        #                                fwhm_method='gaussian', ee_total_radius=ee_total_radius)
        for x, y in zip(x_list, y_list):
            cutout = data[int(y_list[i])-radius:int(y_list[i])+radius,
                           int(x_list[i])-radius:int(x_list[i])+radius]
            # Measure FWHM
            fwhm_x, fwhm_y = measure_fwhm(cutout)
            # Measure EE
            ee_fn = measure_ee(cutout)
            if ee_fn is None:
                print(f'Could not measure EE for {segment}. Recording EE of None.')
                d['encircled_energy_1a'] = None # the encircled energy at 1 arcsec
                d['encircled_energy_.5a'] = None # the encircled energy at half an arcsec
            else:
                ee_at_1arcsec = ee_fn(arcsec)
                ee_at_half_arcsec = ee_fn(arcsec/2)
                d['encircled_energy_1a'] = float(ee_at_1arcsec) # the encircled energy at 1 arcsec
                d['encircled_energy_.5a'] = float(ee_at_half_arcsec) # the encircled energy at half an arcsec

        return fwhm_x, fwhmy, ee
