# STDLIB
from glob import glob
import os

# THIRD PARTY
import numpy as np
from astropy.io import fits
from scipy import signal

# LOCAL
from jwst_fgs_commissioning_tools.nircam_to_fgs import counts_to_jmag
from jwst_fgs_commissioning_tools import log, utils


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

# Constants
NIRCAM_SW_SCALE = 0.031  # NIRCam SW pixel scale (arcsec/pixel)
NIRCAM_LW_SCALE = 0.063  # NIRCam LW pixel scale (arcsec/pixel)
FGS_PIXELS = 2048  # FGS image size in pixels
FGS_PLATE_SIZE = 2.4  # FGS image size in arcseconds

# Constants to change
BAD_PIXEL_THRESH = 2000  # Bad pixel threshold

# -------------------------------------------------------------------------------
def bad_pixel_correction(data, bp_thresh):
    '''Finds and smooths out bad pixels with a median filter'''
    # apply median filter
    smooth = signal.medfilt(data, 3)

    # set negative values to zero
    j = smooth.copy()
    j[j < 0] = 0

    # difference between image and smoothed image; leaves the background behind
    # so we can filter out the bad pixels
    delta = data - smooth

    # Locating the bad pixels. If there are still bpix in masked image, fiddle
    # with delta threshold
    j = np.where(delta > bp_thresh)

    # using location of the bad pixels, replace the bpix value with median value
    # of the smoothed image
    # also get rid of any negative numbers
    data[j] = np.median(smooth)
    data[data < 0] = 0

    # recast as unsigned integers
    data = np.int_(data)

    # clip any over saturated/hot pixels left, replace with integer form of
    # median value of smoothed image
    data[data > 50000] = np.int_(np.median(smooth))

    return data


def rotate_nircam_image(image, fgs_guider, header, nircam_det):
    '''
    Given NIRCAM module A or B (given by the header in your original NIRCAM image),
    rotate/flip to put in correct orientation for FGS 1 and 2.
    '''
    # The Dectector keyword retruns 'NRCA*' or 'NRCB*' so to simplify matters
    # I just pull out the 4th character in the string
    if nircam_det is not None:
        detector = nircam_det
    else:
        detector = header['DETECTOR'][3:].strip()

    # Determine whether the NIRCam image is short- or long-wave to determine
    # the pixel scale
    if '5' in detector:
        # Longwave
        nircam_scale = NIRCAM_LW_SCALE
    else:
        # Shortwave
        nircam_scale = NIRCAM_SW_SCALE

    # Based on the specific detector frame, rotate to match the raw FGS frame
    if detector in ['A2', 'A4', 'B1', 'B3', 'B5']:
        if fgs_guider == 1:
            # FGS guider = 1; Perform a Left-Right flip and swap axes
            image = np.fliplr(image)  # equivalent to image[:,::-1]
            image = np.swapaxes(image, 0, 1)
        elif fgs_guider == 2:
            # FGS guider = 2; Perform a 180 degree rotation and swap axes
            image = np.rot90(image, k=2)
            image = np.swapaxes(image, 0, 1)

    elif detector in ['A1', 'A3', 'A5', 'B2', 'B4']:
        if fgs_guider == 1:
            # FGS guider = 1; Perform a Up-Down flip and swap axes
            image = np.flipud(image)  # equivalent to image[::-1,...]
            image = np.swapaxes(image, 0, 1)
        elif fgs_guider == 2:
            # FGS guider = 2; Swap axes!
            image = np.swapaxes(image, 0, 1)

    else:
        log.error('Unfamiliar NIRCam detector provided. Check the header keyword' +
                  ' "DETECTOR" for the NIRCAM module, then re-run using the ' +
                  '"nircam_det" keyword to bypass the header query.')

    return nircam_scale, image


def pad_data(data, padding, fgs_pix):
    """
    Pad data with mean of data
    """
    size = np.shape(data)[0]

    # Remove NIRCam pedestals
    ped_size = size / 4
    noped_data = np.zeros(np.shape(data))
    for i in range(4):
        ped_start = int(i * ped_size)
        ped_stop = int((i + 1) * ped_size)
        ped_strip = data[:, ped_start:ped_stop]
        pedestal = np.median(ped_strip)

        # Subtract median from each pedestal strip
        noped_data[:, ped_start:ped_stop] = data[:, ped_start:ped_stop] - pedestal
        print('Removing pedestal {} value: {}'.format(i + 1, pedestal))

    # Create an array of size (binned data + 2*padding), filled with the mean data value
    padded_size = size + 2 * (padding)
    if padded_size != fgs_pix - 8:
        # If just a +1 error from odd size of image
        if padded_size == 2039:
            padded_size = 2040
        # If something else is going on....
        else:
            raise ValueError('Padded image not of proper size (should be 2040): {}'.format(padded_size))
    avg_signal = np.mean(noped_data)
    padded_data = np.full((padded_size, padded_size), avg_signal)

    # Replace center of array with real data
    padded_data[padding:padding + size, padding:padding + size] = noped_data

    # Correct high or low pixels
    padded_data[padded_data < 0] = 0.
    padded_data[padded_data > 65000] = 0.

    return padded_data

def resize_nircam_image(data, nircam_scale, fgs_pix, fgs_plate_size):
    '''
    Resize the passed in NIRCam image to match the expected FGS size
    '''
    cropped = data[4:-4, 4:-4]  # crop 4pixel zero-padding
    binned_pix = int(round((data.shape[0] * nircam_scale * fgs_pix) / (fgs_plate_size * 60)))
    data_resized = utils.resize_array(cropped, binned_pix, binned_pix)

    padding = int((cropped.shape[0] - binned_pix) / 2)
    data_pad = pad_data(data_resized, padding, fgs_pix)
    fgs_data = np.pad(data_pad, 4, 'constant')  # Add back reference pixels

    return fgs_data


def normalize_data(data, fgs_counts, threshold=5):
    '''
    Threshold of 5 assumes background is very low. *This will need to be automated
    later.*
    '''
    mask = data > threshold
    data_norm = np.copy(mask * data.astype(np.float64))
    data_norm *= (fgs_counts / data_norm.sum())  # renormalize by sum of non-masked data
    data_norm[mask == 0] = data[mask == 0]  # background is not normalized

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
    bias_data = fits.getdata(bias_data_path)

    binned_pad_norm_bias = fgs_data + bias_data

    if save_to_fits:
        utils.ensure_dir_exists(os.path.join(output_path, 'bin_norm_bias_imgs'))

        if guider is None:
            guider = bias_data_path.split('/')[-1].split('.')[0][-6:]
        biasout_path = os.path.join(output_path, 'bin_norm_bias_imgs',
                                '{}_G{}_binned_pad_norm.fits'.format(root, guider))
        utils.write_fits(biasout_path, binned_pad_norm_bias)

    return binned_pad_norm_bias


# -------------------------------------------------------------------------------
# -------------------------------------------------------------------------------
def convert_im(input_im, guider, fgs_counts=None, jmag=None, nircam_det=None,
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
    nircam_det: str
        The NIRCAM module, otherwise the header will be parsed
    '''

    # Establish paths and necessary files
    local_path = os.path.dirname(os.path.realpath(__file__))  # where this script exists
    package_path = os.path.split(local_path)[0]  # where the package exists
    out_path = os.path.split(package_path)[0]  # where the out/ dir goes
    data_path = os.path.join(package_path, 'data')  # Includes data/*.fits files (ie newmagicHdrImg,bias0, etc)
    header_file = os.path.join(data_path, 'newG{}magicHdrImg.fits'.format(guider))  # Guider-dependent files

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
    all_ims = []
    for image in im_list:
        basename = os.path.basename(image)

        root = basename.split('.')[0]
        log.info('Beginning to create FGS image from {}'.format(root))

        if output_path is None:
            output_path_save = os.path.join(out_path, 'out', root)
            utils.ensure_dir_exists(output_path_save)
        else:
            output_path_save = output_path

        data = fits.getdata(image, header=False)
        header = fits.getheader(image, ext=0)

        # ---------------------------------------------------------------------
        # Create FGS image
        # Mask out bad pixels
        data_masked = bad_pixel_correction(data, BAD_PIXEL_THRESH)
        # Rotate the NIRCAM image into FGS frame
        nircam_scale, data_rot = rotate_nircam_image(data_masked, guider, header, nircam_det)
        # Pad image
        data_pad = resize_nircam_image(data_rot, nircam_scale, FGS_PIXELS, FGS_PLATE_SIZE)
        # Normalize image
        data_norm = normalize_data(data_pad, fgs_counts)

        fgsout_path = os.path.join(output_path_save, 'FGS_imgs',
                                '{}_G{}_binned_pad_norm.fits'.format(root, guider))
        # Any value about 65535 will wrap when converted to uint16
        data_norm[data_norm >= 65535] = 65535
        dummy_data, hdr = fits.getdata(header_file, header=True)
        # print hdr, dummy_data
        utils.write_fits(fgsout_path, np.uint16(data_norm), header=hdr)

        print("Finished for {}, Guider = {}".format(root, guider))

        all_ims.append(data_norm)

    if return_im:
        return all_ims
