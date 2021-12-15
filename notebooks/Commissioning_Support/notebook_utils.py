"""
This module supports JWST MAGIC notebooks, in particular the Pre-LOS-02 Guiding Analysis notebook.

Included in this notebook are functions and a class that supports guiding analysis for CARs prior
to LOS-02.
"""

import os

from astropy.io import fits
from astropy.modeling import models, fitting
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
# These values are determined by the Get simulated NIRCam locations of LOS-02 PSFs section of the Pre-LOS-02 notebook
GA_PSF_LOCATIONS = {'A1-1':[1037, 599], 'A2-2':[1391, 808], 'A3-3':[1375, 1218],
                    'A4-4':[1020, 1426], 'A5-5':[670, 1204], 'A6-6':[660, 783],
                    'B1-7':[1047, 177], 'B2-9':[1756, 608], 'B3-11':[1771, 1459],
                    'B4-13':[1010, 1789], 'B5-15':[297, 1409], 'B6-17':[327, 589],
                    'C1-8':[1394, 404], 'C2-10':[1748, 1024], 'C3-12':[1370, 1637],
                    'C4-14':[650, 1630], 'C5-16':[309, 987], 'C6-18':[690, 391]
                    }
NIRCAM_SW_PIXELSCALE = 0.031 # arcsec/pixel

class PsfAnalysis():
    def __init__(self, car, image=None, header=None, pixelscale=None, target_location=None):
        """
        The PsfAnalysis class takes in an image and associated information in order to perform PSF
        analysis on each PSF in that image.

        Parameters:
        -----------
        car: string
            The name of the CAR - this must match how it is labelled in the out path
        image: 2D array
            This is generally expected to be a mosaic image
        header: fits header object
            A header that is associated with the image, it does not need to be exactly for that
            image if that image is made up of many other images.
        target_location: tuple
            The estimated location of the target in the image if measured. If not known, the
            average of the found PSFs will be taken as this target_location

        Attributes:
        -----------
        car:string
            The name of the CAR - this must match how it is labelled in the out path
        header: fits header object
            A header that is associated with the image. It can be the header associated with one
            of the images that makes up the mosaic.
        pixelscale: float
            The pixelscale associated with this image (if known)
        target_location: tuple
            The estimated location of the target in the image
        corrected_image: 2D array
            A negative value-corrected version of the input image that will be used for all analysis
        info_df: pandas DataFrame object
            A table of the segment PSF characteristics for each PSF found in the input image
        psfs: list of 2D arrays
            A list of the cutout PSFs that have been identified in the input image
        psf_window_radius: int
            Half the size of one side of the box used to cut out the PSFs. This will change based
            on the pixelscale
        ee_fns: list of scipy interpolation objects
            A list of the enclided energy scipy interpolation objects for each PSF that has been
            identified.
        """
        self.car = car
        self.header = header
        if pixelscale is None:
            self.pixelscale = NIRCAM_SW_PIXELSCALE
        else:
            self.pixelscale = pixelscale

        self.target_location = target_location

        if self.header is not None:
            self.check_header_for_pixelscale()
        if image is not None:
            # Correct any negative values
            self.corrected_image = correct_pixel_values(image)
        else:
            self.corrected_image = None

        # Set up attributes that will be used in the other methods
        self.info_df = None # The data frame holding the PSF information
        self.psfs = None # List of cutout PSFs
        self.psf_window_radius = None # The radius of the box used for cutting out PSFs
        self.ee_fns = None # list of EE interpolation functions

    def check_header_for_pixelscale(self):
        """Make sure that if the header does not already include the PIXELSCL keyword, add it.
        This header keyword is required in order to measure the enircled energy
        """
        if 'PIXELSCL' not in self.header.keys() and self.pixelscale is not None:
            self.header['PIXELSCL'] = self.pixelscale
        elif 'PIXELSCL' in self.header.keys() and self.pixelscale is None:
            self.pixelscale = self.header['PIXELSCL']

    def get_image_information(self, smoothing='high', psf_window_radius=40):
        """
        In the input image, find all PSFs, cut out the PSFs and make a list of them, and
        measure the FWHM in x and y, distance from target, and create a list of EE interpolation
        function. Add the PSF locations, FWHM, and distance to target information to an
        information data frame. The PSFs and EE interpolation functions are returned as lists.
        """
        # Save the radius used to cut out the postage stamp as it will be used later
        self.psf_window_radius = psf_window_radius
        # Find each PSF in the input image
        x_list, y_list = get_position_from_magic(self.corrected_image, smoothing=smoothing)
        print(f'{len(x_list)} PSFs found')

        psfs = []
        for x, y in zip(x_list, y_list):
            psfs.append(cut_out_psf(x, y, self.corrected_image, psf_window_radius))
        self.create_info_df(x_list, y_list, psfs)

    def create_info_df(self, x_list, y_list, psfs=None):
        """
        For each PSF, measure the FWHM, distance from target, and create a EE function,
        then add to a information data frame.
        """
        # Get PSF characteristis for each PSF
        if x_list:
            fwhm_xs, fwhm_ys, self.ee_fns = self.get_psf_characteristics(psfs)
            self.get_target_location(x_list, y_list)
            distances_to_target = self.get_distances_to_target(x_list, y_list)
            # Create the data frame
            self.info_df = pd.DataFrame(data={'x': x_list, 'y': y_list,
                                              'fwhm_x': fwhm_xs, 'fwhm_y': fwhm_ys,
                                              'distance_to_target': distances_to_target})
        else:
            self.info_df = None

    def get_psf_characteristics(self, psfs):
        '''
        Get the FWHM in x and y and the ee function for each identified PSF
        '''
        # Make sure that if we don't have the pixelscale value at this point, that the
        #  user knows that we can't measure the EE
        if 'PIXELSCL' not in self.header.keys():
            print('EE cannot be measured because there is no pixelscale information.')
        fwhm_xs = []
        fwhm_ys = []
        ee_fns = []
        for cutout in psfs:
            # Measure FWHM
            fwhm_x, fwhm_y = measure_fwhm(cutout)
            fwhm_xs.append(fwhm_x)
            fwhm_ys.append(fwhm_y)
            # Measure EE
            try:
                ee_fn = measure_ee(cutout, self.header)
            except KeyError:
                ee_fn = None
            ee_fns.append(ee_fn)

        return fwhm_xs, fwhm_ys, ee_fns

    def get_target_location(self, x_list, y_list):
        """ If we don't have a target location in the image, calculate it from the
        locations of segment PSFs that we do have
        """
        if self.target_location is None:
            self.target_location = (int(np.median(x_list)), int(np.median(y_list)))

    def get_distances_to_target(self, x_list, y_list):
        """Calculate the distance to the target location for each known PSF
        """
        target_x, target_y = self.target_location
        distances_to_target = distance_to_target_center(x_list, y_list,
                                                        target_x, target_y)
        return distances_to_target

    def plot_mosaic_with_psfs(self, xlim=None, ylim=None,
                              label_color='white', legend_location='best'):
        """
        Plot out the mosaic image with the identified PSFs and estimated target location plotted
        The info_df attribute cannot be None before trying to run this function.
        """
        if self.info_df is None:
            print("Cannot plot anything without the info_df being populated.")
            return
        x_list = self.info_df['x']
        y_list = self.info_df['y']
        if xlim is None:
            xlim = (min(x_list)-100, max(x_list)+100)
        if ylim is None:
            ylim = (min(y_list)-100, max(y_list)+100)
        try:
            seg_list = self.info_df['segment']
        except KeyError:
            seg_list = None
        # Plot it!
        plt.figure(figsize=(14, 14))
        plt.imshow(self.corrected_image, norm=LogNorm(vmin=1, vmax=1000), origin='lower')
        plt.scatter(self.target_location[0], self.target_location[1], s=500, marker='*', c='white',
                    edgecolor='C1', label='Estimated Target Location')
        plt.legend(loc=legend_location)
        if x_list is not None and y_list is not None:
            for j, (x, y) in enumerate(zip(x_list, y_list)):
                name = seg_list[j] if seg_list is not None else j
                plt.annotate(name, xy=(x, y), xytext=(x-100, y), c=label_color, fontsize=14,
                             weight='bold')
            plt.xlim(xlim)
            plt.ylim(ylim)

        plt.title(f"{self.car} Mosaic")
        plt.show()

    def plot_each_identified_segment_psf(self, num_psf_per_row=6, labels=None):
        """
        Cut each segment out of the mosaic/image provided based on the x and y locations given with
        x_list, and y_list
        """
        if self.psfs is None:
            print("Cannot plot anything without the PSFs list being populated.")
            return

        if labels is None:
            labels = np.arange(len(self.psfs))

        rows = len(self.psfs)//num_psf_per_row+1
        columns = len(self.psfs) if len(self.psfs) < num_psf_per_row else num_psf_per_row
        _, ax = plt.subplots(rows, columns, figsize=(4*columns, 4*rows))

        i = 0
        for j in range(rows):
            for k in range(columns):
                try:
                    ax[j, k].imshow(self.psfs[i], origin='lower', norm=LogNorm(vmin=1))
                    ax[j, k].set_title(f'PSF {labels[i]}')
                    i += 1
                except IndexError:
                    ax[j, k].axis('off')
        plt.tight_layout()
        plt.show()

    def plot_multiple_parameters(self, parameters, xs=None, ys=None, xlim=None, ylim=None):
        """
        Give a list of parameters that you want to plot information for. This plot will show the
        value of each at the location of the PSF that is provided in the info_dictionary, but if the
        parameters xs and ys are provided, those will be used instead. When we are plotting the PSFs
        in their GA location, we must provide the GA xs and ys in order to overide the position of
        each PSF in the original image.
        """
        if self.info_df is None:
            print("Cannot plot anything without the info_df being populated.")
            return
        # Check that all parameters are in info_dictionary
        parameters = list(set(parameters) & set(list(self.info_df.columns)))
        if not parameters:
            print('Requested parameters are not in the povided data frame. Nothing to plot')
        else:
            fig, ax = plt.subplots(1, len(parameters), figsize=(6*len(parameters), 8))
            try:
                segs = self.info_df['segment'][self.info_df['segment'] != 'Unknown'].values
            except KeyError:
                segs = np.arange(len(self.info_df))

            if xs is None and ys is None:
                xs = self.info_df['x'].values
                ys = self.info_df['y'].values

            for i, param in enumerate(parameters):
                axs = ax if len(parameters) == 1 else ax[i]
                for seg, x, y in zip(segs, xs, ys):
                    axs.annotate(seg, (x, y), (x+30, y+30), fontsize=16)
                try:
                    values = self.info_df[param][self.info_df['segment'] != 'Unknown'].values
                except KeyError:
                    values = self.info_df[param].values

                pl = axs.scatter(xs, ys, marker='o', s=100, c=values, cmap=cm.RdYlGn_r)
                fig.colorbar(pl, ax=axs, orientation='horizontal')
                axs.set_title(f"{param} for each PSF in the image", fontsize=16)

                if xlim and ylim:
                    axs.set_xlim(xlim)
                    axs.set_ylim(ylim)

    def plot_ee(self, num_psf_per_row=3, labels=None):
        '''If we have EE measurements, plot them
        '''
        if self.pixelscale is not None:
            ee_pixel_radius = self.psf_window_radius
            ee_arcsec_radius = ee_pixel_radius * self.pixelscale
            xs = np.linspace(0, ee_arcsec_radius)

            num_data_points = len(self.ee_fns)
            if labels is None:
                labels = np.arange(num_data_points)

            # Just double check that we have an EE function for this PSF
            rows = num_data_points//num_psf_per_row+1
            columns = num_data_points if num_data_points < num_psf_per_row else num_psf_per_row
            _, ax = plt.subplots(rows, columns, figsize=(4*columns, 4*rows))

            i = 0
            for j in range(rows):
                for k in range(columns):
                    try:
                        ys_total = self.ee_fns[i](xs)
                        ys_percent = ys_total/ys_total.max()*100
                        ax[j, k].plot(xs, ys_percent)
                        ax[j, k].set_title(f'PSF {labels[i]}')
                        ax[j, k].set_xlabel("Radius [arcsec]")
                        ax[j, k].set_ylabel("Encircled Energy (percent)")
                        ax[j, k].set_ylim(None, 110)
                        i += 1
                    except IndexError:
                        ax[j, k].axis('off')
            plt.tight_layout()
            plt.show()


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
    Calculate the distance from a segment PSF to the expected location of the target in the sky.
    truth_x and truth_y are the location of the boresight
    '''
    return np.sqrt((truth_x - np.asarray(x))**2 + (truth_y - np.asarray(y))**2)


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


def cut_out_psf(x, y, image, psf_window_radius):
    """ Cut out each found PSF """
    cutout = image[int(y)-psf_window_radius:int(y)+psf_window_radius,
                   int(x)-psf_window_radius:int(x)+psf_window_radius]
    return cutout


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
    p_init = models.Gaussian2D(amplitude=array.max(), x_mean=xp*0.5, y_mean=yp*0.5)
    fit_p = fitting.LevMarLSQFitter()
    fitted_psf = fit_p(p_init, x, y, array)
    return fitted_psf.x_fwhm, fitted_psf.y_fwhm


def measure_ee(array, header):
    '''Wrapper function around poppy's measure_ee
    '''
    hdu1 = fits.PrimaryHDU(data=array, header=header)
    new_hdul = fits.HDUList([hdu1])
    return poppy.measure_ee(new_hdul, normalize='None')


def create_combined_pointing_df(info_dfs, matching_dictionary, ee_fn_lists):
    """Given multiple info_dfs, combine to include the information for all identified
    PSFs. This also means creating a new list of encircled energy functions.
    `info_dfs` is a list of the info_dfs in order of pointing (pointing1 first)
    `ee_fn_lists` is a list of the ee_fns in order of pointing (pointing1 first)
    """
    new_dfs = []
    new_ee_fns = []
    for info_df, pointing, ee_fns in zip(info_dfs, matching_dictionary.keys(), ee_fn_lists):
        pointing_dict = matching_dictionary[pointing]
        selected_indicies = list(pointing_dict.keys())
        new_df = info_df.loc[selected_indicies]
        # Make the new EE fns list
        new_ee_fns += [ee_fns[i] for i in selected_indicies]
        # Grab the right list of segments and add to the new dataframe
        segments = matching_dictionary[pointing].values()
        new_df['segment'] = segments
        # Make sure the pointing is recorded in the new dataframe
        new_df['pointing'] = np.repeat(pointing, len(segments))
        new_dfs.append(new_df)

    combined_info_df = pd.concat(new_dfs)
    combined_info_df = combined_info_df[['segment', 'pointing', 'x', 'y', 'fwhm_x', 'fwhm_y',
                                         'distance_to_target']]

    return combined_info_df, new_ee_fns

def get_psfs_from_all_pointings(combined_df, pointings_images, psf_window_radius):
    """
    Create a new list of PSF cutouts from each pointingfor a full list of each segment PSF
    """
    x_list = combined_df['x']
    y_list = combined_df['y']
    pointing_list = combined_df['pointing']

    psfs = []
    for x, y, pointing in zip(x_list, y_list, pointing_list):
        pointing_num = pointing.split('_')[-1]
        image = pointings_images[int(pointing_num)-1]
        psfs.append(cut_out_psf(x, y, image, psf_window_radius))

    return psfs


def separate_nircam_images(all_images_list):
    """
    Separate a list of NIRCam images to the NIRCam A and NIRCam B images
    """
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
    From a list of NIRCam A or B images, create a list of the data and filenames
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
    Plot the 4 SW NIRCam A images in their respecitve postions
    """

    a1, a2, a3, a4 = data_list
    a1_name, a2_name, a3_name, a4_name = name_list

    _, ax = plt.subplots(2, 2, figsize=(12, 12))
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
    Plot the 4 SW NIRCam B images in their respecitve postions
    """

    b1, b2, b3, b4 = data_list
    b1_name, b2_name, b3_name, b4_name = name_list

    _, ax = plt.subplots(2, 2, figsize=(12, 12))
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
    Create a mosaic array from 4 NIRCam A and 4 NIRCam B images.
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


def create_image_array(psfs, segment_list, window_size):
    """
    Cut out postage stamps of an image and using the locations of each PSF and it's assumed
    segment ID, place it in a large image array.
    """
    # Create an array of zeros the size of a standard dector size
    large_image_array = np.zeros([2048, 2048])

    ga_xs = []
    ga_ys = []
    for psf, seg in zip(psfs, segment_list):
        if seg != 'Unknown':
            ga_x, ga_y = GA_PSF_LOCATIONS[seg]
            ga_xs.append(ga_x)
            ga_ys.append(ga_y)
            large_image_array[ga_y-window_size: ga_y+window_size,
                              ga_x-window_size:ga_x+window_size] = psf

    large_image_array = correct_pixel_values(large_image_array)
    return large_image_array, ga_xs, ga_ys


def save_out_image_for_magic(image, header, filename, out_dir):
    """
    Write data to a fits file. Can take in an array/header
    or lists of arrays/headers
    """
    outfile = os.path.join(out_dir, filename)

    # Create a fake DQ array that treats every pixel as good
    dq_array = np.zeros_like(image)

    # Make a list of the data and header
    data = [image, dq_array]
    header['DETECTOR'] = 'NRCA3'
    headers = [header, None]

    header_list = header if isinstance(header, list) else [header]
    for hdr in header_list:
        if not any([isinstance(hdr, fits.header.Header), hdr is None]):
            raise TypeError(f'Header to be written out in {outfile} is not either "None" or of type fits.header.Header')

    hdu_list = []
    for i, (dat, hdr) in enumerate(zip(data, headers)):
        if i == 0:
            hdu = fits.PrimaryHDU(data=dat, header=hdr)
        else:
            hdu = fits.ImageHDU(data=dat, header=hdr, name='DQ')
        hdu_list.append(hdu)
    hdul = fits.HDUList(hdu_list)

    hdul.writeto(outfile, overwrite=True)
