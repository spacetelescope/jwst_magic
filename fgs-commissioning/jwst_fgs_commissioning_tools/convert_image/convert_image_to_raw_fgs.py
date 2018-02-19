# STDLIB
from glob import glob
import itertools
import os

# THIRD PARTY
import numpy as np
from astropy.io import fits
from scipy import signal

# LOCAL
from jwst_fgs_commissioning_tools.convert_image import counts_to_jmag
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
# Paths
FSW_PATH = os.path.dirname(os.path.realpath(__file__))
PACKAGE_PATH = os.path.split(FSW_PATH)[0]
OUT_PATH = os.path.split(PACKAGE_PATH)[0]  # Location of out/ and logs/ directory
DATA_PATH = os.path.join(PACKAGE_PATH, 'data')

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
    delta = data - smooth

    # Locating the bad pixels.
    j = np.where(delta > bp_thresh)

    # using location of the bad pixels, replace the bpix value with median value
    # of the smoothed image
    data[j] = np.median(smooth)

    # clip any over saturated/hot pixels left, replace with integer form of
    # median value of smoothed image
    data = utils.correct_image(data, upper_threshold=50000, upper_limit=np.median(smooth))
    # recast as unsigned integers
    data = np.uint16(data)

    return data

def correct_nircam_dq(image, dq_array, bit_arr=None):
    '''
    Based on a conversation with Alicia Canipe on the NIRCam team on 02/09/2018,
    try focusing on these flags
    Bit Value Name      Description
    0	1	    DO_NOT_USE	 Bad pixel. Do not use.
    9	512	    NON_SCIENCE	 Pixel not on science portion of detector
    10	1024	DEAD	     Dead pixel
    13	8192	LOW_QE	     Low quantum efficiency
    16	65536	NONLINEAR	 Pixel highly nonlinear
    19	524288	NO_GAIN_VALUE	Gain cannot be measured
    20	1048576	NO_LIN_CORR	 Linearity correction not available
    21	2097152	NO_SAT_CHECK Saturation check not available
    '''
    # Convert bits into values
    if bit_arr is None:
        bit_arr = [0, 9, 13, 16, 19, 20, 21]

    flags = [2**x for x in bit_arr]
    # Find all combinations of bits
    flag_combinations = [seq for i in range(len(flags), 0, -1) for seq in itertools.combinations(flags, i)]

    inds = []
    for flag_comb in flag_combinations:
        # Add all combinations together to find all values we care about
        flag = (np.sum(flag_comb))
        # Now create DQ array that just flags the pixels with values that we care about
        ind = list(zip(*np.where(dq_array == flag)))
        if ind:
            inds.extend(ind)

    # Fix bad pixel by taking median value in 3x3 box
    im_copy = np.copy(image)
    for ind in inds:
        im_copy[ind] = np.median(image[ind[0]-2:ind[0]+2, ind[1]-2:ind[1]+2])

    return im_copy

def fgs_add_dq(image, guider):
    '''
    Add FGS bad pixels to image
    Currently, we only have a map of all flagged pixels but no indication as to
    why they are flagged. For now, we set all flagged pixels to saturation.
    '''
    if guider == 1:
        dq_arr = fits.getdata(os.path.join(DATA_PATH, 'fgs_dq_G1.fits'))
    elif guider == 2:
        dq_arr = fits.getdata(os.path.join(DATA_PATH, 'fgs_dq_G2.fits'))

    # Apply dq_arr to image
    # FIXME for now, set all flagged pixels to saturation
    image[dq_arr == 1] = 65535

    return image


def nircam_raw_to_fgs_raw(image, nircam_detector, fgs_guider):
    '''
    rotate image from NIRCam detector (raw) coordinate frame to FGS raw

    (See 'notebooks/Convert from NIRCam to FGS coordinate frames.ipynb')
    '''
    # Based on the specific detector frame, rotate to match the raw FGS frame
    if nircam_detector in ['A2', 'A4', 'B1', 'B3', 'B5']:
        if fgs_guider == 1:
            # FGS guider = 1; Perform 270 degree rotation
            image = np.rot90(image, k=3)
        elif fgs_guider == 2:
            # FGS guider = 2; Perform a 180 degree rotation and swap axes
            image = np.rot90(image, k=2)
            image = np.swapaxes(image, 0, 1)

    elif nircam_detector in ['A1', 'A3', 'A5', 'B2', 'B4']:
        if fgs_guider == 1:
            # FGS guider = 1; One 90 degree rotation
            image = np.rot90(image, k=1)
        elif fgs_guider == 2:
            # FGS guider = 2; Swap axes!
            image = np.swapaxes(image, 0, 1)

    else:
        log.error('Unfamiliar NIRCam detector provided. Check the header keyword' +
                  ' "DETECTOR" for the NIRCAM module, then re-run using the ' +
                  '"nircam_det" keyword to bypass the header query.')

    return image

def sci_to_fgs_raw(image, fgs_guider):
    '''
    Rotate image from DMS science coordinate (same for NIRCam and FGS) frame to FGS raw
    ** This is the expected frame for output DMS images **

    (See 'notebooks/Convert from NIRCam to FGS coordinate frames.ipynb' and
    'notebooks/Convert FGS coordinate frames.ipynb')
    '''
    if fgs_guider == 1:
        # FGS guider = 1; Swap axes
        image = np.swapaxes(image, 0, 1)
    elif fgs_guider == 2:
        # FGS guider = 2; Perform 90 degree rotation
        image = np.rot90(image, k=1)

    return image

def rotate_nircam_image(image, fgs_guider, header, nircam_det,
                        nircam_coord_frame='sci'):
    '''
    Given NIRCAM module A or B (given by the header in your original NIRCAM image),
    rotate/flip to put in correct orientation for FGS 1 and 2.

    Parameters:
    -----------
    image: array-like
        NIRCam image to be rotated into correct FGS frame
    fgs_guider: int
        Guider 1 or 2
    header: .fits header object
        The header of the input NIRCam image
    nircam_det: str
        The NIRCam detector with which the image was taken.
        Expects: A1, A2, A3, A4, A5, B1, B2, B3, B4, or B5
    nircam_coord_frame: str
        The coordinate frame that the input image is in.
        Expects: 'sci' for the science frame, or 'raw' or 'det' for the raw
        detector frame

    Returns:
    --------
    nircam_scale: float
        Detector scale factor depending on if the input image is from a long- or
        shortwave detector
    image: array-like
        The rotated NIRCam image
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

    if nircam_coord_frame == 'sci':
        image = sci_to_fgs_raw(image, fgs_guider)
    elif nircam_coord_frame == 'raw' or nircam_coord_frame == 'det':
        image = nircam_raw_to_fgs_raw(image, detector, fgs_guider)
    else:
        raise ValueError('Unrecognized coordinate frame name.')

    return nircam_scale, image


def pad_data(data, padding, fgs_pix):
    """
    Pad data with mean of data
    """
    size = np.shape(data)[0]

    # Remove NIRCam pedestals
    ped_size = size // 4
    noped_data = np.zeros(np.shape(data))
    for i in range(4):
        ped_start = i * ped_size
        ped_stop = (i + 1) * ped_size
        ped_strip = data[:, ped_start:ped_stop]
        pedestal = np.median(ped_strip)

        # Subtract median from each pedestal strip
        noped_data[:, ped_start:ped_stop] = data[:, ped_start:ped_stop] - pedestal
        # print('Removing pedestal {} value: {}'.format(i + 1, pedestal))

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
    padded_data = utils.correct_image(padded_data, 65000, 0)

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

# -------------------------------------------------------------------------------
# -------------------------------------------------------------------------------
def convert_im(input_im, guider, nircam=True, fgs_counts=None, jmag=None, nircam_det=None,
               out_dir=None):
    '''
    Takes NIRCam image and turns it into an FGS-like image, gets count rate and location of
    each guide star in each image

    Parameters
    ==========
    input_im: str
        The path to the input image
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
    # For the images requested, convert to FGS images
    basename = os.path.basename(input_im)

    root = basename.split('.')[0]
    log.info('Beginning to create FGS image from {}'.format(root))

    output_path_save = utils.make_out_dir(out_dir, OUT_PATH, root)
    utils.ensure_dir_exists(output_path_save)

    data = fits.getdata(input_im, header=False)
    header = fits.getheader(input_im, ext=0)

    if len(data.shape) > 2:
        print(data.shape)
        raise TypeError('Expecting a single frame or slope image.')

    # ---------------------------------------------------------------------
    # Create FGS image
    # Mask out bad pixels
    # data = bad_pixel_correction(data, BAD_PIXEL_THRESH)

    if nircam:
        log.info("This is a NIRCam image")

        # Pull out DQ array for this image
        dq_arr = fits.getdata(input_im, extname='DQ')
        if not dq_arr.min() == 1 and not dq_arr.max() == 1:
            data = correct_nircam_dq(data, dq_arr)

        # Rotate the NIRCAM image into FGS frame
        nircam_scale, data = rotate_nircam_image(data, guider, header, nircam_det)
        # Pad image
        data = resize_nircam_image(data, nircam_scale, FGS_PIXELS, FGS_PLATE_SIZE)

    else:
        log.info("This is an FGS image")
        guider = utils.get_guider(header)
        try:
            origin = header['ORIGIN'].strip()
            if origin == 'ITM':
                log.info("Data provided in science/DMS frame; rotating to raw FGS frame.")
                data = sci_to_fgs_raw(data, guider)
        except KeyError:
            pass

    # Normalize image
    data_norm = normalize_data(data, fgs_counts)

    # Any value above 65535 or below 0 will wrap when converted to uint16
    data_norm = utils.correct_image(data_norm, upper_threshold=65535, upper_limit=65535)

    # Load header file
    header_file = os.path.join(DATA_PATH, 'newG{}magicHdrImg.fits'.format(guider))
    fgsout_path = os.path.join(output_path_save, 'FGS_imgs',
                               '{}_G{}.fits'.format(root, guider))

    hdr = fits.getheader(header_file, ext=0)
    utils.write_fits(fgsout_path, np.uint16(data_norm), header=hdr)

    print("Finished for {}, Guider = {}".format(root, guider))

    return data_norm
