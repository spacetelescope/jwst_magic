"""Convert NIRCam or FGS images to raw FGS frame.

This tool takes input NIRCam or FGS images, re-bins to FGS plate scale
if necessary, applies the DQ array if necessary, rotates the image
according to the SIAF coordinate transformations, and re-normalizes to
the desired FGS countrate/magnitude. The development of this code used
mock-up NIRCam files from Ball Aerospace's ITM tool, simulating Global
Alignment using the short wavelength channel.

Authors
-------
    - Keira Brooks
    - Lauren Chambers
    - Shannon Osborne


Use
---
    This module can be executed in a Python shell as such:
    ::
        from jwst_magic.convert_image import convert_image_to_raw_fgs
        convert_image_to_raw_fgs.convert_im(input_im, guider, root):

    Required arguments:
        ``input_image`` - filepath for the input (NIRCam or FGS) image
        ``guider`` - number for guider 1 or guider 2
        ``root`` - will be used to create the output directory, ./out/{root}
    Optional arguments:
        ``nircam`` - denotes if the input_image is an FGS or NIRCam
            image. If True, the image will be converted to FGS format.
        ``nircam_det`` - used to specify the detector of a provided
            NIRCam image. If left blank, the detector will be extracted
            from the header of the NIRCam FITS file.
        ``normalize`` - denotes if the image will be normalized.
        ``norm_value`` and ``norm_unit`` - If the image will be
            normalized, specifies the value to normalize to and the
            units of that value (either FGS Magnitude or FGS countrate).
        ``out_dir`` - where output files will be saved. If not provided,
            the image(s) will be saved within the repository at
            jwst_magic/
        ``coarse_pointing`` - denotes if the image will have a Gaussian
            filter applied to simulate the effects of jitter when the
            observatory is in coarse pointing rather than fine guide.
        ``jitter_rate_arcsec`` - the rate of the spacecraft jitter, in
            arcseconds, that will be used to apply the Gaussian filter
            if coarse_pointing is True.
        ``logger_passed`` - denotes if a logger object has already been
            generated.

References
----------
    See JWST-STScI-001550 (Rev A), "Description and Use of the JWST
    Science Instrument Aperture File," on SOCCER for definition of
    coordinate transformations between NIRCam and FGS frames, and
    between DMS/Science frames and raw frames.

Notes
-----
    NIRCam short wave field-of-view: 1.09' x 1.09'
    NIRCam short wave pixel scale: 0.032"/pixel

    FGS field-of-view: 2.4' x 2.4'
    FGS pixel scale: 0.070"/pixel
"""

# Standard Library Imports
import copy
import itertools
import logging
import os

# Third Party Imports
from astropy.io import fits
import fgscountrate
import numpy as np
import pysiaf
from scipy import signal
from scipy.ndimage.filters import gaussian_filter

# Local Imports
from jwst_magic.utils import utils

# Paths
__location__ = os.path.realpath(os.path.join(os.getcwd(), os.path.dirname(__file__)))
PACKAGE_PATH = os.path.split(__location__)[0]
OUT_PATH = os.path.split(PACKAGE_PATH)[0]  # Location of out/ and logs/ directory
DATA_PATH = os.path.join(PACKAGE_PATH, 'data')

# Constants
NIRCAM_SW_SCALE = 0.031  # NIRCam SW pixel scale (arcsec/pixel)
NIRCAM_LW_SCALE = 0.063  # NIRCam LW pixel scale (arcsec/pixel)
FGS_SCALE = 0.069  # FGS pixel scale (arcsec/pixel)
FGS_PIXELS = 2048  # FGS image size in pixels
FGS_PLATE_SIZE = 2.4  # FGS image size in arcseconds

# Start logger
LOGGER = logging.getLogger(__name__)

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# SUPPORTING FUNCTIONS
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


def apply_coarse_pointing_filter(data, jitter_rate_arcsec, pixel_scale):
    """Apply a Gaussian filter to simulate coarse pointing

    Parameters
    ----------
    data : 2-D numpy array
        Image data
    jitter_rate_arcsec : float
        The rate of jitter of the telescope (arcsec/sec)
    pixel_scale : float
        The pixel scale of the detector (arcsec/pixel)

    Returns
    -------
    data_gauss : 2-D numpy array
        Image data with coarse pointing filter applied
    """
    t_fullframe_read = 10.7  # sec
    jitter_rate = jitter_rate_arcsec / pixel_scale  # pixel/sec
    sigma = jitter_rate * t_fullframe_read / 3  # pixel

    data_gauss = gaussian_filter(data, sigma)

    return data_gauss


def bad_pixel_correction(data, bp_thresh):
    """Finds bad pixels that are above a threshold; smooths them with a
    3-pixel median filter.

    Parameters
    ----------
    data : 2-D numpy array
        Image data
    bp_thresh : int
        Threshold above which pixels are considered "bad" and smoothed

    Returns
    -------
    data : 2-D numpy array
        Image data with pixel correction applied
    """
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
    """Apply a data quality array correction to an input NIRCam image.

    Parameters
    ----------
    image : 2-D numpy array
        Image data
    dq_array : 2-D numpy array
        Data quality array
    bit_arr : list, optional
        List of DQ flags that indicate pixels that should be fixed (see Notes)

    Returns
    -------
    im_copy : 2-D numpy array
        Image data with DQ array corrections applied

    Notes
    -----
    Based on a conversation with Alicia Canipe on the NIRCam team on 02/09/2018,
    try focusing on these flags:

    Bit  Value     Name             Description
    --   -----     ----             -----------
    0    1         DO_NOT_USE       Bad pixel. Do not use.
    9    512       NON_SCIENCE      Pixel not on science portion of detector
    10   1024      DEAD             Dead pixel
    13   8192      LOW_QE           Low quantum efficiency
    16   65536     NONLINEAR        Pixel highly nonlinear
    19   524288    NO_GAIN_VALUE    Gain cannot be measured
    20   1048576   NO_LIN_CORR      Linearity correction not available
    21   2097152   NO_SAT_CHECK     Saturation check not available
    """

    # Convert bits into values
    if bit_arr is None:
        bit_arr = [0, 9, 13, 16, 19, 20, 21]

    flags = [2**x for x in bit_arr]
    # Find all combinations of bits
    flag_combinations = [seq for i in range(len(flags), 0, -1)
                         for seq in itertools.combinations(flags, i)]

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
        im_copy[ind] = np.median(image[ind[0] - 2:ind[0] + 2,
                                       ind[1] - 2:ind[1] + 2])

    return im_copy


def transform_nircam_raw_to_fgs_raw(image, from_nircam_detector, to_fgs_detector):
    """Transform image from NIRCam detector raw/det coordinate frame to FGS
    raw, using transformations as defined by the Science Instrument
    Aperture File (SIAF).

    Parameters
    ----------
    image : 2-D numpy array
        Input image data
    from_nircam_detector : str
        Name of NIRCam detector of the input image (A1, A2, A3, A4, A5,
        B1, B2, B3, B4, or B5)
    to_fgs_detector : int
        Guider number of the desired output image (1 or 2)

    Returns
    -------
    image : 2-D numpy array
        Image data with coordinate frame transformation applied
    """

    # 1) Transform to NIRCam sci
    # Get the Det2Sci angle and parity from the SIAF
    from_nircam_detector = 'NRC' + from_nircam_detector
    nircam_siaf = pysiaf.Siaf('NIRCam')
    aperture = nircam_siaf['{}_FULL'.format(from_nircam_detector)]
    angle = aperture.DetSciYAngle
    parity = aperture.DetSciParity

    # Flip the X axis according to the parity
    if parity == -1:
        image = np.fliplr(image)

    # Rotate the image according to the angle
    n_90_deg_rots = angle // 90
    image = np.rot90(image, k=n_90_deg_rots)

    # 2) Transform to FGS raw
    image = transform_sci_to_fgs_raw(image, to_fgs_detector)

    return image

def transform_sci_to_fgs_raw(image, to_fgs_detector):
    """Rotate NIRCam or FGS image from DMS/science coordinate frame
    (the expected frame for output DMS images) to FGS raw. Note that
    it is not necessary to specify the input image detector because
    the DMS coordinate frame is identical for all FGS and NIRCam
    detectors.

    Parameters
    ----------
    image : 2-D numpy array
        Input image data
    to_fgs_detector : int
        Guider number of the desired output image (1 or 2)

    Returns
    -------
    image : 2-D numpy array
        Image data with coordinate frame transformation applied

    """
    # Get the Det2Sci angle and parity from the SIAF
    to_fgs_detector = 'FGS' + str(to_fgs_detector)
    fgs_siaf = pysiaf.Siaf('FGS')
    aperture = fgs_siaf['{}_FULL'.format(to_fgs_detector)]
    angle = aperture.DetSciYAngle
    parity = aperture.DetSciParity

    # Flip the X axis according to the parity
    if parity == -1:
        image = np.fliplr(image)

    # Rotate the image according to the angle
    n_90_deg_rots = angle // 90
    image = np.rot90(image, k=n_90_deg_rots)

    # Because goal is FGS raw (not det), swap X and Y
    image = np.swapaxes(image, 0, 1)

    return image


def transform_nircam_image(image, to_fgs_detector, from_nircam_detector, header,
                           nircam_coord_frame='sci'):
    """Given NIRCam image and detector, rotate and flip to put in
    correct orientation for FGS 1 or 2.

    Parameters
    ----------
    image : 2-D numpy array
        NIRCam image to be rotated into correct FGS frame
    to_fgs_detector : int
        Guider 1 or 2
    header : astropy.io.fits.Header object
        The header of the input NIRCam image
    from_nircam_detector : str
        The NIRCam detector with which the image was taken.
        Expects: A1, A2, A3, A4, A5, B1, B2, B3, B4, or B5
    nircam_coord_frame : str, optional
        The coordinate frame that the input image is in.
        Expects: 'sci' for the science/DMS frame, or 'raw' or 'det' for
        the raw detector frame

    Returns
    --------
    nircam_scale : float
        Detector pixel scale depending on if the input image is from a
        long- or shortwave detector
    image : 2-D numpy array
        The transformed NIRCam image
    """
    # If the NIRCam detector isn't specified, get it from the header
    if from_nircam_detector in [None, "-Parse from Header-"]:
        # The Detector keyword returns 'NRCA*' or 'NRCB*', so to simplify matters
        # just pull out the 4th & 5th character in the string
        from_nircam_detector = header['DETECTOR'][3:].strip()
    LOGGER.info("Image Conversion: Transforming from NIRCAM Detector = {}".format(from_nircam_detector))

    # Determine whether the NIRCam image is short- or long-wave to determine
    # the pixel scale
    nircam_scale = NIRCAM_LW_SCALE if '5' in from_nircam_detector else NIRCAM_SW_SCALE

    # Perform the transformation
    if nircam_coord_frame == 'sci':
        LOGGER.info("Image Conversion: Input NIRCam image in SCI coordinate frame.")
        image = transform_sci_to_fgs_raw(image, to_fgs_detector)

    elif nircam_coord_frame == 'raw' or nircam_coord_frame == 'det':
        LOGGER.info("Image Conversion: Input NIRCam image in RAW/DET coordinate frame.")
        image = transform_nircam_raw_to_fgs_raw(image, from_nircam_detector, to_fgs_detector)

    else:
        raise ValueError('Unrecognized coordinate frame name.')

    return nircam_scale, image


def pad_data(data, padding, fgs_pix):
    """Pad re-binned NIRCam data with mean of data; effectively placing
    NIRCam data onto an appropriately-sized FGS array with the
    appropriate pixel scale.

    Parameters
    ----------
    data : 2-D numpy array
        Image data
    padding : int
        Width of padded area around data
    fgs_pix : int
        Number of pixels along one side of an FGS image (probably 2048)

    Returns
    -------
    padded_data
        Image data padded to match FGS pixel scale

    Raises
    ------
    ValueError
        If the output FGS image is not going to be 2048x2048 (i.e. if
        the padding width does not match the size of the input image
        data)
    """
    # Determine the size of the data array
    size = np.shape(data)[0]

    # Create an array of size (binned data + 2*padding), filled with the mean data value
    padded_size = size + 2 * padding
    if padded_size != fgs_pix - 8:
        # If just a +1 error from odd size of image
        if padded_size == 2039:
            padded_size = 2040
        # If something else is going on....
        else:
            raise ValueError('Padded image not of proper size (should be 2040): {}'.format(padded_size))
    mean_signal = np.median(data)
    padded_data = np.full((padded_size, padded_size), mean_signal)

    # Replace center of array with real data
    padded_data[padding:padding + size, padding:padding + size] = data

    # Correct high or low pixels
    padded_data = utils.correct_image(padded_data, 65000, 0)

    return padded_data


def resize_nircam_image(data, nircam_scale, fgs_pix, fgs_plate_size):
    """Resize a NIRCam image to the expected FGS size and pixel scale

    Parameters
    ----------
    data : 2-D numpy array
        Image data
    nircam_scale : float
        Pixel scale of NIRCam detector
    fgs_pix : int
        Number of pixels along one side of an FGS image (probably 2048)
    fgs_plate_size : float
        Pixel scale of FGS detector

    Returns
    -------
    fgs_data
        Re-binned and padded image data
    """
    cropped = data[4:-4, 4:-4]  # crop 4pixel zero-padding
    binned_pix = int(round((data.shape[0] * nircam_scale * fgs_pix) / (fgs_plate_size * 60)))
    data_resized = utils.resize_array(cropped, binned_pix, binned_pix)

    padding = int((cropped.shape[0] - binned_pix) / 2)
    data_pad = pad_data(data_resized, padding, fgs_pix)
    fgs_data = np.pad(data_pad, 4, 'constant')  # Add back reference pixels

    return fgs_data


def normalize_data(data, fgs_countrate):
    """Re-normalize data to the desired FGS countrate without
    masking out any of the background. See JWSTFGS-160 for the
    reasoning behind this.

    Parameters
    ----------
    data : 2-D numpy array
        Image data
    fgs_countrate : float
        The FGS countrate value to normalize to

    Returns
    -------
    data_norm
        Normalized image data

    """

    data_norm = np.copy(data.astype(np.float64))
    data_norm *= (fgs_countrate / data_norm.sum())

    return data_norm


def remove_pedestal(data, nircam, itm):
    """Subtract the vertical (for NIRCam) or horizontal (for FGS)
     pedestals/amps from a raw frame. The pedestal is calculated
     as the median of the amps if the image is ITM and the median
     of the reference pixels if it is not.

    Parameters
    ----------
    data : 2-D numpy array
        Image data
    nircam : bool
        True if the data is NIRCam SCI frame data. False if
        the data is FGS SCI frame data.
    itm : bool
        True if the data is an ITM image (which won't have
        accurate reference pixels)

    Returns
    -------
    noped_data : 2-D numpy array
        Image data with pedestals subtracted
    """
    size = np.shape(data)[0]
    ped_size = size // 4
    pedestals = []

    if itm is True:  # Use the median of the amps and subtract it from the whole amp
        noped_data = np.zeros(np.shape(data))
        for i in range(4):
            ped_start = i * ped_size
            ped_stop = (i + 1) * ped_size

            # Subtract median from each pedestal strip
            if nircam:
                ped_strip = data[:, ped_start:ped_stop]
                pedestal = np.median(ped_strip)
                pedestals.append(pedestal)
                noped_data[:, ped_start:ped_stop] = data[:, ped_start:ped_stop] - pedestal
            else:
                ped_strip = data[ped_start:ped_stop, :]
                pedestal = np.median(ped_strip)
                pedestals.append(pedestal)
                noped_data[ped_start:ped_stop, :] = data[ped_start:ped_stop, :] - pedestal

    # For CV3/Other hardware data (that includes reference pixels)
    else:  # Use the median of the refpix and subtract it only from the non-refpix on the amp
        noped_data = copy.deepcopy(data)  # retain the refpix - will overwrite everything else
        for i in range(4):
            ped_start = i * ped_size
            ped_stop = (i + 1) * ped_size

            # Get pedestal start/stop points excluding the reference pixels
            if i == 0:
                start = ped_start + 4  # for left rows of refpix
                stop = ped_stop
            elif i == 3:
                start = ped_start
                stop = ped_stop - 4  # for right rows of refpix
            else:
                start = ped_start
                stop = ped_stop

            # Pull each pedestral strip (FGS vs NIRCam pedestal runs different ways)
            # Subtract median of 4px refpix from each pedestal strip (not including the reference pixels)
            if nircam:
                pedestal = np.median([data[0:4, start:stop], data[-4:, start:stop]])
                pedestals.append(pedestal)
                noped_data[4:-4, start:stop] = data[4:-4, start:stop] - pedestal
            else:
                pedestal = np.median([data[start:stop, 0:4], data[start:stop, -4:]])
                pedestals.append(pedestal)
                noped_data[start:stop, 4:-4] = data[start:stop, 4:-4] - pedestal

    LOGGER.info("Image Conversion: " +
                "Removed pedestal values from image: {} ".
                format(', '.join(['{:.2f}'.format(p) for p in pedestals])))

    return noped_data

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# MAIN FUNCTIONS
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


def convert_im(input_im, guider, root, nircam=True,
               nircam_det=None, normalize=True, norm_value=12.0,
               norm_unit="FGS Magnitude", gs_catalog=None, coarse_pointing=False,
               jitter_rate_arcsec=None, logger_passed=False, itm=False):
    """Takes NIRCam or FGS image and converts it into an FGS-like image.

    Parameters
    ----------
    input_im : str
        Filepath for the input (NIRCam or FGS) image
    guider : int
        Guider number (1 or 2)
    root : str
        Name used to create the output directory, {out_dir}/out/{root}
    nircam : bool, optional
        Denotes if the input_image is an FGS or NIRCam image. If True,
        the image will be converted to FGS format. Unless out_dir is
        specified, the FGS-formatted image will be saved to
        ../out/{root}/FGS_imgs/{root}_binned_pad_norm.fits
    nircam_det : str, optional
        The detector of a provided NIRCam image. If left blank, the
        detector will be extracted from the header of the NIRCam FITS
        file.
    normalize : bool, optional
        Denotes if the image will be normalized. If True, norm_value
        and norm_unit will be used to determine the normalization value
    norm_value : str or float, optional
        Specifies the Guide Star ID or the count rate/magnitude to which
        to normalize.
    norm_unit : str, optional
        Specifies the unit of norm_value ("FGS Magnitude", "FGS countrate",
        or "Guide Star ID")
    gs_catalog : str, optional
        Guide star catalog version to query. E.g. 'GSC242'. None will use
        the default catalog as defined in teh FGS Count Rate Module.
    coarse_pointing : bool, optional
        Denotes if the image will have a Gaussian filter applied to
        simulate the effects of jitter when the observatory is in
        coarse pointing rather than fine guide.
    jitter_rate_arcsec : None, optional
        The rate of the spacecraft jitter, in arcseconds per second,
        that will be used to apply the Gaussian filter if
        coarse_pointing is True.
    logger_passed : bool, optional
        Denotes if a logger object has already been generated.
    itm : bool, optional
        If this image come from the ITM simulator (important for normalization).

    Returns
    -------
    data : 2-D numpy array
        Image formatted like a raw FGS image

    Raises
    ------
    TypeError
        The input filename has more than one frame.
    ValueError
        An input NIRCam file has an obstruction in the pupil.
    """
    # Start logging
    if not logger_passed:
        utils.create_logger_from_yaml(__name__, root=root, level='DEBUG')

    try:
        LOGGER.info("Image Conversion: " +
                    "Beginning image conversion to guider {} FGS image".format(guider))
        LOGGER.info("Image Conversion: Input image is expected to be in units of ADU/sec (countrate)")

        data = fits.getdata(input_im, header=False)
        header = fits.getheader(input_im, ext=0)

        if len(data.shape) > 2:
            raise TypeError('Expecting a single frame or slope image.')

        # Check if this is an ITM image and the itm flag is set correctly (backwards compatibility)
        try:
            origin = header['ORIGIN'].strip()
            if origin == 'ITM':
                try:
                    assert itm is True
                except AssertionError:
                    itm = True
                    LOGGER.warning("Deprecation Warning: This is an ITM image, setting itm flag to 'True'")
        except KeyError:
            origin = None

        # Try to check that the units on the input image are as expected (Dn/s = ADU/s; *_rate.fits)
        try:
            header_sci = fits.getheader(input_im, extname='sci')
        except KeyError:
            header_sci = {}

        for hdr in [header, header_sci]:
            if 'BUNIT' in hdr:
                input_unit = hdr['BUNIT'].lower()
            if 'PHOTMJSR' in hdr:
                photmjsr = hdr['PHOTMJSR']

        try:
            if input_unit == 'mjy/sr':
                convert_to_adu_s = photmjsr
                data /= convert_to_adu_s
                LOGGER.info('Image Conversion: Input is a Cal image. Converting from MJy/sr to ADU/s')
            elif input_unit == 'dn/s':
                LOGGER.info('Image Conversion: Image in correct units of ADU/s.')
        except NameError:
            LOGGER.info("Image Conversion: Can't check image type because of missing "
                        "BUNIT keyword. User should confirm the input image is a rate image.")
            pass

        # Create raw FGS image...

        # Remove pedestal from NIRCam or FGS data
        # pedestal should be taken out in refpix correction - only run if that hasn't been run
        # and if not labeled as test data which is made without a pedestal
        if 'S_REFPIX' in header.keys() and header['S_REFPIX'] == 'COMPLETE':
            LOGGER.info("Image Conversion: Skipping removing pedestal - Reference pixel correction run in pipeline.")
        elif 'TEST' in header.keys() and header['TEST'] == 'True':
            LOGGER.info("Image Conversion: Skipping removing pedestal - Test flag found in image header.")
        else:
            data = remove_pedestal(data, nircam, itm)

        # -------------- From NIRCam --------------
        if nircam:
            LOGGER.info("Image Conversion: This is a NIRCam image")

            # Check that the pupil is clear
            try:
                pupil_keyword = header['PUPIL']
                if pupil_keyword in ['CLEAR', 'Imaging Pupil']:
                    pass
                else:
                    raise ValueError(
                        'NIRCam "PUPIL" header keyword for provided file is {}. '.format(pupil_keyword) +
                        'Only the CLEAR/Imaging Pupil can be used to realistically simulate FGS images.'
                    )
            except KeyError:
                pass

            # Pull out DQ array for this image
            try:
                dq_arr = fits.getdata(input_im, extname='DQ')
                if not dq_arr.min() == 1 and not dq_arr.max() == 1:
                    data = correct_nircam_dq(data, dq_arr)
            except KeyError:
                LOGGER.warning("Image Conversion: No DQ extension found; DQ correction cannot be performed.")

            # Rotate the NIRCAM image into FGS frame
            nircam_scale, data = transform_nircam_image(data, guider, nircam_det, header)
            # Pad image
            data = resize_nircam_image(data, nircam_scale, FGS_PIXELS, FGS_PLATE_SIZE)

        # -------------- From FGS --------------
        else:
            LOGGER.info("Image Conversion: This is an FGS image")
            # Check if header keyword is equal to fgs raw to determine if rotation to raw is needed
            if origin is not None and origin.upper() == 'FGSRAW':  # this is a magic-team keyword that should only be in test files
                LOGGER.info("Image Conversion: Data is already provided in raw frame; no rotation done")
                LOGGER.warning("Assume input guider is same as output guider; no rotation done")
            else:
                LOGGER.info(
                    "Image Conversion: Expect that data provided is in science/DMS frame; rotating to raw FGS frame.")
                data = transform_sci_to_fgs_raw(data, guider)

        # Apply Gaussian filter to simulate coarse pointing
        if coarse_pointing:
            pixel_scale = nircam_scale if nircam else FGS_SCALE

            data = apply_coarse_pointing_filter(data, jitter_rate_arcsec, pixel_scale)
            LOGGER.info("Image Conversion: Applied Gaussian filter to simulate "
                        "coarse pointing with jitter of {:.3f} arcsec/sec".format(jitter_rate_arcsec))

        # Normalize the image, if the "normalize" flag is True
        # The ITM simulations are only created for relative SNR so they need to
        # normalized to one before anything else happens
        if itm:
            LOGGER.info("Image Conversion: This is an ITM image.")
            data -= data.min() #set minimum at 0.
            data /= data.sum()  # set total countrate to 1.
            if not norm_value:
                norm_value = 12
                norm_unit = 'FGS Magnitude'
                LOGGER.warning("Image Conversion: No normalization was specified but is required for an ITM image. Using FGS Magnitude of 12.")

        if normalize or itm:
            # Convert magnitude/countrate to FGS countrate using new count rate module
            # Take norm_value and norm_unit to pass to count rate module
            if norm_unit == 'FGS countrate' and isinstance(norm_value, (float, int)): # If already a count rate, do nothing
                if norm_value < 1000:
                    LOGGER.warning("Image Conversion: Small count rate has been requested. This may be an error.")
                fgs_countrate = norm_value
                fgs_mag = ''
            elif norm_unit == 'FGS countrate' and isinstance(norm_value, str):
                LOGGER.error("Image Conversion: Value type does not match expectation for unit type.")
                raise TypeError("Mismatched normalization value ({}) and expected type ({}) for unit of {}".format(
                    norm_value, "str", norm_unit))
            else: # If a GSID or magnitude, pass back gsid
                # Sneaky string replacement for FGS Magnitude
                norm_unit = 'FGS Magnitude' if norm_unit == 'FGS Mag (10, 11, 12, 13, or 14)' else norm_unit
                gsid = check_norm_value_unit(norm_value, norm_unit)
                LOGGER.info("Image Conversion: Using GSID {} to normalize image.".format(gsid))
                fgs = fgscountrate.FGSCountrate(guide_star_id=gsid, guider=guider)
                fgs_countrate, _, fgs_mag, _ = fgs.query_fgs_countrate_magnitude(catalog=gs_catalog)

            # Normalize the data
            data = normalize_data(data, fgs_countrate)
            LOGGER.info("Image Conversion: Normalizing to {} FGS Countrate (FGS Mag: {})".format(fgs_countrate,
                                                                                                 fgs_mag))

    except Exception as e:
        LOGGER.exception(e)
        raise

    return data

def check_norm_value_unit(norm_value, norm_unit):
    '''
    For the case when the norm_unit is not "FGS countrate", get back the
    guide star ID for the norm_value given.

    Parameters
    ----------
    norm_value : str or float, optional
        Specifies the Guide Star ID or the count rate/magnitude to which to normalize.
    norm_unit : str, optional
        Specifies the unit of norm_value ("FGS Magnitude", "FGS countrate", or "Guide Star ID")

    Returns
    -------
    gsid: str
        The guide star ID to be used to get the FGS count rate

    '''
    if norm_unit == 'Guide Star ID' and isinstance(norm_value, str):
        gsid = norm_value
    elif norm_unit == 'FGS Magnitude' and isinstance(norm_value, (float,int)):
        gsid = query_jmag_mini_library(norm_value)
    else:
        LOGGER.error("Image Conversion: Value type does not match expectation for unit type.")
        raise TypeError("Mismatched normalization value ({}) and expected type ({}) for unit of {}".format(
            norm_value, ["float or int" if norm_unit == 'FGS Magnitude' else 'str'], norm_unit))
    return gsid


def query_jmag_mini_library(fgsmag):
    '''
    Since the count rate module requires a Guide Star ID (GSID), we have
    pre-defined GSIDs for Guide Stars with magnitudes between (approximately)
    FGS Magnitude of 10 and 14.
    For a FGS Mag of 10, 11, 12, 13, or 14, return the following GSID:
    10: N135000314 (FGS Mag (G1) = 10.09767, J Mag = 9.99699)
    11: N13A000006 (FGS Mag (G1) = 11.16569, J Mag = 11.00100)
    12: N13A000158 (FGS Mag (G1) = 12.16066, J Mag = 12.03899)
    13: N13A000125 (FGS Mag (G1) = 13.07528, J Mag = 13.00800)
    14: N13A002729 (FGS Mag (G1) = 13.95341, J Mag = 13.98499)

    Parameters
    ----------
    fgsmag: int
        Can only be int between (inclusive)10 and 14.

    Raises
    ------
    ValueError
        The input is not an integer between (inclusive) 10 and 14.

    '''
    accepted_mags = np.arange(10, 15)
    if fgsmag not in accepted_mags:
        LOGGER.error("Image Conversion: FGS Magnitude not value. Please enter an integer between 10 and 14.")
        raise ValueError("Unacceptable FGS Magnitude value")
    else:
        LOGGER.info("Image Conversion: Using pre-defined Guide Star ID for FGS Magnitude of: {}".format(fgsmag))

    fgsmag_dict = {
        10: 'N135000314',
        11: 'N13A000006',
        12: 'N13A000158',
        13: 'N13A000125',
        14: 'N13A002729'
        }

    return fgsmag_dict[fgsmag]

def write_fgs_im(data, out_dir, root, guider, fgsout_path=None):
    """Writes an array of FGS data to the appropriate file:
    {out_dir}/out/{root}/FGS_imgs/{root}_G{guider}.fits

    Parameters
    ----------
    data : 2-D numpy array
        FGS image data
    out_dir : str, optional
        Where output files will be saved. If not provided, the
        image(s) will be saved within the repository at
        jwst_magic/
    root : str
        Name used to create the output directory, {out_dir}/out/{root}
    guider : int
        Guider number (1 or 2)
    fgsout_path : str, optional
        Alternate directory in which to save the FGS files. If not
        provided, the FGS images will be saved to
        {out_dir}/out/{root}/FGS_imgs/

    Returns
    -------
    fgsout_path : str
        Filepath for the output FGS image
    """
    data = utils.correct_image(data, upper_threshold=65535, upper_limit=65535)

    # Define output path
    output_path_save = utils.make_out_dir(out_dir, OUT_PATH, root)
    utils.ensure_dir_exists(output_path_save)
    if not fgsout_path:
        fgsout_path = os.path.join(output_path_save, 'FGS_imgs')
    fgsout_file = os.path.join(fgsout_path, 'unshifted_{}_G{}.fits'.format(root, guider))

    # Load header file
    header_file = os.path.join(DATA_PATH, 'newG{}magicHdrImg.fits'.format(guider))
    hdr = fits.getheader(header_file, ext=0)

    # Write FITS file
    utils.write_fits(fgsout_file, data, header=hdr, log=LOGGER)

    return fgsout_path
