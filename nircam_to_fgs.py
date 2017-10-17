# STDLIB
import os

# Third Party
import numpy as np
from astropy.io import fits
from glob import glob
from scipy import ndimage, signal
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import matplotlib
from matplotlib import cycler
matplotlib.rcParams['image.cmap'] = 'viridis'
matplotlib.rcParams['axes.prop_cycle'] = cycler(u'color',['#1f77b4',
            '#ff7f0e',
            '#2ca02c',
            '#d62728',
            '#9467bd',
            '#8c564b',
            '#e377c2',
            '#7f7f7f',
            '#bcbd22',
            '#17becf'])
matplotlib.rcParams['image.origin'] = 'upper'

# LOCAL
import counts_to_jmag
import log
import select_psfs
import utils


'''
Python Img Re-Binning Tool for FGS comissioning data

Input NIRCam images, rebin to FGS platescale and renormalize to desired FGS

NIRCam FOV: 1.09' by 1.09'
NIRCam pixel scale: 0.032"/pixel

FGS FOV: 2.4' by 2.4'
FGS pixel scale: 0.070"/pixel


For development of this code, used mock-up NIRCam files from Ball Aerospace's ITM tool,
simulating Global Alignment using the short wavelength channel

(Sherie was told they include detector noise and bias offsets between the readout channels)
'''


#-------------------------------------------------------------------------------
def bad_pixel_correction(data, bp_thresh):
    '''Finds and smooths out bad pixels with a median filter'''
    #apply median filter
    smooth = signal.medfilt(data, 3)

    #set negative values to zero
    j = smooth.copy()
    j[j < 0] = 0

    #difference between image and smoothed image; leaves the background behind
    # so we can filter out the bad pixels
    delta = data-smooth

    #Locating the bad pixels. If there are still bpix in masked image, fiddle
    # with delta threshold
    j = np.where(delta > bp_thresh)

    #using location of the bad pixels, replace the bpix value with median value
    # of the smoothed image
    #also get rid of any negative numbers
    data[j] = np.median(smooth)
    data[data < 0] = 0

    #recast as unsigned integers
    data = np.int_(data)

    #clip any over saturated/hot pixels left, replace with integer form of
    # median value of smoothed image
    data[data > 50000] = np.int_(np.median(smooth))

    return data


def rotate_nircam_image(im, fgs_guider, header, nircam_mod):
    '''
    Given NIRCAM module A or B (given by the header in your original NIRCAM image),
    rotate/flip to put in correct orientation for FGS 1 and 2.
    '''
    # The Dectector keyword retruns 'NRCA*' or 'NRCB*' so to simplify matters
    # I just pull out the 4th character in the string
    if nircam_mod is not None:
        module = nircam_mod
    else:
        module = header['DETECTOR'][3]

    if module == 'A':
        ## NIRCAM Module A
        if fgs_guider == 1:
            ## FGS guider = 1; Perform a Left-Right flip
            im = np.fliplr(im)# equivalent to im[:,::-1]
        else:
            ## FGS guider = 2; Perform a 180 degree rotation
            im = np.rot90(im, k=20)

    elif module == 'B':
        ## NIRCAM Module B
        if fgs_guider == 1:
            ## FGS guider = 1; Perform a Up-Down flip
            im = np.flipud(im)# equivalent to im[::-1,...]
        else:
            ## FGS guider = 2; No change necessary!
            pass
    else:
        log.error('Check the header keyword "DETECTOR" for the NIRCAM module, \
              then re-run using the "nircam_mod" keyword to bypass the header query.')
    return im


def pad_data(data, padding):
    """
    Pad data with median of data with Poisson noise
    """
    size = np.shape(data)[0]
    print(size)

    # Create an array of size binned data + 2*padding
    padded_size = size + 2 * (padding)
    print(padded_size)
    background = np.zeros((padded_size, padded_size))

    # Remove NIRCam pedestals
    ped_size = size / 4
    peddata = np.zeros((size, size))
    for i in range(4):
        ped_start = i * ped_size
        ped_stop = (i + 1) * ped_size
        ped_strip = data[:, ped_start:ped_stop]
        pedestal = np.median(ped_strip)

        # Subtract median from each pedestal strip
        peddata[:, ped_start:ped_stop] = data[:, ped_start:ped_stop] - pedestal
        print('Removing pedestal {} value: {}'.format(i + 1, pedestal))

    # Add Poisson noise
    padded_data = np.random.poisson(lam=5, size=background.shape)

    # Replace center of array with real data
    padded_data[padding:padding + size, padding:padding + size] = data

    # Add in FGS pedestal
    fgs_ped = np.fix(15 * np.random.random_sample(size=4)).astype(int)  # Randomly generate values
    ped_size = padded_size / 4
    for i in range(4):
        ped_start = i * ped_size
        ped_stop = (i + 1) * ped_size
        pedestal = fgs_ped[i]

        # Add randomly generated pedestal to each padded pedestal strip
        padded_data[:, ped_start:ped_stop] += pedestal

    return padded_data

def resize_nircam_image(data, nircam_scale, fgs_pix, fgs_plate_size):
    cropped = data[4:-4, 4:-4]  # crop 4pixel zero-padding
    binned_pix = int(round((data.shape[0] * nircam_scale * fgs_pix) / (fgs_plate_size * 60)))

    data_resized = utils.resize_array(cropped, binned_pix, binned_pix)

    padding = int((cropped.shape[0] - binned_pix) / 2)
    data_pad = pad_data(data_resized, padding)
    fgs_data = np.pad(data_pad, 4, 'constant')  # Add back reference pixels

    return fgs_data


def normalize_data(data, fgs_counts, threshold=5):
    '''
    Threshold of 5 assumes background is very low. *This will need to be automated
    later.*
    '''
    mask = data > threshold
    data_norm = np.copy(mask*data.astype(np.float64))
    data_norm *= (fgs_counts/data_norm.sum()) #renormalize by sum of non-masked data
    data_norm[mask == 0] = data[mask == 0] #background is not normalized

    return data_norm



def add_bias_to_data(bias_data_path, fgs_data, root, guider='', output_path='',
                     save_to_fits=True):
    """
    OUT OF DATE - 6/21/17
    Adds in the bias from the guider (two seperate files) to FGS
    (or simulated FGS image that has been padded and normalized)

    This assumes that the format of the guider bias filename is
       job<ID>_g<guider_number>bias.fits.
    If you have a different filename format, pass in the guider_name.

    Default guider bias files are found in the following location:
       guider1: "g1bias.fits"
       guider2: "g2bias.fits"
    """
    bias_data = utils.read_fits(bias_data_path)[1]

    binned_pad_norm_bias = fgs_data + bias_data

    if save_to_fits:
        utils.ensure_dir_exists(os.path.join(output_path, 'bin_norm_bias_imgs'))

        if guider is None:
            guider = bias_data_path.split('/')[-1].split('.')[0][-6:]
        out_path =  os.path.join(output_path, 'bin_norm_bias_imgs',
                      '{}_G{}_binned_pad_norm.fits'.format(root,guider))
        utils.write_fits(out_path, binned_pad_norm_bias)

    return binned_pad_norm_bias



#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
def convert_im(input_im, guider, fgs_counts=None, jmag=None, nircam_mod=None,
               return_im=True, output_path=None):
    '''
    Takes NIRCam image and turns it into an FGS-like image, gets count rate and location of
    each guide star in each image

    Parameters
    ==========
    input_im: str,list of strings
        This can be the path to a file, the path to a directory where all .fits
        files are images you want to convert, or a list of paths (i.e. use glob to
        get this list)
        input_path = '/user/kbrooks/itar/FGS/ga/subset_060617'
    guider: int
        1 or 2
    fgs_counts: int
        The FGS counts in the star. If set to 'None', pass in the jmag. SET TO 65353 FOR NOW
    jmag: float
        The J magnitude of the star. If set to 'None' and fgs_counts set to 'None',
        will defaul to 11.
    nircam_mod: str
        The NIRCAM module, otherwise the header will be parsed
    '''

    # Establish paths and necessary files
    # 'local_path' is the path where this script exists
    local_path = os.path.dirname(os.path.realpath(__file__))
    # Data path is the directory that includes *.fits files (ie newmagicHdrImg,bias0, etc)
    data_path = os.path.join(local_path, 'data')
    # Guider-dependent files
    header_file = os.path.join(data_path, 'newG{}magicHdrImg.fits'.format(guider))
    # bias_data_path = os.path.join(data_path, 'g{}bias0.fits'.format(guider))

    # ---------------------------------------------------------------------
    # Constants
    nircam_scale = 0.032  # NIRCam pixel scale
    fgs_pix = 2048  # FGS image size in pixels
    fgs_plate_size = 2.4  # FGS image size in arcseconds
    # Constants to change
    bp_thresh = 2000  # Bad pixel threshold

    # ---------------------------------------------------------------------
    # Find FGS counts to be used for normalization
    if fgs_counts is None:
        if jmag is None:
            log.warning('No counts or J magnitude given, setting to default')
            jmag = 11
        fgs_counts = counts_to_jmag.jmag_to_fgs_counts(jmag, guider)
    else:
        jmag = counts_to_jmag.fgs_counts_to_jmag(fgs_counts, guider)

    log.info('J magnitude = {:.1f}'.format(jmag))

    # ---------------------------------------------------------------------
    # Get list of images from input path: can take file,list,dir
    if isinstance(input_im, list):
        im_list = input_im
    elif os.path.isfile(input_im):
        im_list = [input_im]
    elif os.path.isdir(input_im):
        im_list = (glob(os.path.join(input_im, '*.fits')))
    else:
        log.error("Input format not recognized. Exiting.")
        return

    # ---------------------------------------------------------------------
    # For the images requested, convert to FGS images
    for im in im_list:
        basename = os.path.basename(im)
        root = basename.split('.')[0]
        log.info('Beginning to create FGS image from {}'.format(root))

        if output_path is None:
            output_path = os.path.join(local_path, 'out', root)
            utils.ensure_dir_exists(output_path)

        header, data = utils.read_fits(im)

        # ---------------------------------------------------------------------
        # Create FGS image
        # Mask out bad pixels
        data_masked = bad_pixel_correction(data, bp_thresh)
        # Rotate the NIRCAM image into FGS frame
        data_rot = rotate_nircam_image(data_masked, guider, header, nircam_mod)
        # Pad image
        data_pad = resize_nircam_image(data_rot, nircam_scale, fgs_pix, fgs_plate_size)
        # Normalize image
        data_norm = normalize_data(data_pad, fgs_counts)

        out_path = os.path.join(output_path, 'FGS_imgs',
                                '{}_G{}_binned_pad_norm.fits'.format(root, guider))
        # Any value about 65535 will wrap when converted to uint16
        data_norm[data_norm >= 65535] = 65535
        hdr = utils.read_fits(header_file)[0]
        utils.write_fits(out_path, np.uint16(data_norm), header=hdr)

        print("Finished for {}, Guider = {}".format(root, guider))

        if return_im:
            return data_norm
