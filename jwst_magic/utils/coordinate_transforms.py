"""Convert between FGS raw/native frame (pixels), ideal angle frame
(arcsec), and DHAS frame (arcsec).

Note that the input *cal.fits images used in MAGIC are all undistorted
immediately, so whenever the "raw", "det" or "sci" frames are discussed
here, they are not the actual frames from the SIAF definition, but are
instead an undistorted frame where the coodinte system matches the "raw",
"det" or "sci" frames from the SIAF respectively.

Authors
-------
    - Lauren Chambers
    - Shannon Osborne

Use
---
    ::
        from jwst_magic.coordinate_transforms import raw2idl
        x_idl, y_idl = raw2idl(x_raw, y_raw)
"""

from astropy import units as u
import numpy as np
import pysiaf

# Open SIAF with pysiaf
FGS_SIAF = pysiaf.Siaf('FGS')


def nrcpixel_offset_to_v2v3_offset(x_offset, y_offset, detector):
    """Convert a boresight offset from NIRCam pixels to V2/V3 arcsec

    Parameters
    ----------
    x_offset : float
        Boresight offset in NIRCam X pixels
    y_offset : float
        Boresight offset in NIRCam Y pixels
    detector : str
        NIRCam detector to transform. E.g. 'NRCA3'

    Returns
    -------
    v2_offset, v3_offset : tup
        Boresight offset in V2/V3 (arcsec)
    """
    # Get pixel scale
    nrc_siaf = pysiaf.Siaf('NIRCam')
    nrc_det = nrc_siaf[f'{detector}_FULL']
    nircam_x_scale = nrc_det.XSciScale  # arcsec/pixel
    nircam_y_scale = nrc_det.YSciScale  # arcsec/pixel

    # Convert x/y offsets to V2/V3
    v2_offset = x_offset * nircam_x_scale  # arcsec
    v3_offset = y_offset * nircam_y_scale  # arcsec

    return v2_offset, v3_offset


def raw2idl(x_raw, y_raw, guider):
    """Pass in undistorted X and Y pixels in the raw/native coordinate frame
     and get out X Y angles in the ideal frame

    Parameters
    ----------
    x_raw : float
        X pixels in the raw detector frame
    y_raw : float
        Y pixels in the raw detector frame
    guider : int
        Which FGS detector to convert for (1 or 2)

    Returns
    -------
    x_idealangle : float
        X angle (arcsec) in the ideal frame
    y_idealangle : float
        Y angle (arcsec) in the ideal frame
    """
    if int(guider) == 1:
        fgs_full = FGS_SIAF['FGS1_FULL']
    elif int(guider) == 2:
        fgs_full = FGS_SIAF['FGS2_FULL']
    else:
        raise ValueError('Unrecognized guider number: {}'.format(guider))

    # Add +1 to incoming values to convert from python array 0-indexing to pixel coordinates 1-indexing
    x_raw_px = x_raw + 1
    y_raw_px = y_raw + 1

    # Convert from RAW -> SCI (just a coordinate origin change, no distortion change)
    x_sci, y_sci = fgs_full.raw_to_sci(x_raw_px, y_raw_px)

    # Shift to the IDL origin location
    x_idl_pix = x_sci - fgs_full.XSciRef
    y_idl_pix = y_sci - fgs_full.YSciRef

    # Convert from IDL pixels to idl arcseconds
    x_idealangle = x_idl_pix * fgs_full.XSciScale
    y_idealangle = y_idl_pix * fgs_full.YSciScale

    return x_idealangle, y_idealangle


def raw2tel(x_raw, y_raw, guider):
    """Pass in undistorted X and Y pixels in the raw/native coordinate frame
     and get out angles in the V2/V3 frame

    Parameters
    ----------
    x_raw : float
        X pixels in the raw detector frame
    y_raw : float
        Y pixels in the raw detector frame
    guider : int
        Which FGS detector to convert for (1 or 2)

    Returns
    -------
    v2 : float
        Angle (arcsec) in the V2 frame
    v3 : float
        Angle (arcsec) in the V3 frame
    """
    if int(guider) == 1:
        fgs_full = FGS_SIAF['FGS1_FULL']
    elif int(guider) == 2:
        fgs_full = FGS_SIAF['FGS2_FULL']
    else:
        raise ValueError('Unrecognized guider number: {}'.format(guider))

    # Convert from Raw to IDL (1 added to raw value in raw2idl to go from python 0-based to pixel coord 1-based)
    x_idl, y_idl = raw2idl(x_raw, y_raw, guider)

    # Convert from IDL -> TEL (just a coordinate origin change, no distortion change)
    v2, v3 = fgs_full.idl_to_tel(x_idl, y_idl)

    return v2, v3


def idl2dhas(x_idealangle, y_idealangle):
    """Pass in X and Y angles in the ideal frame and get out X and Y angles in the
    frame DHAS requires.

    Parameters
    ----------
    x_idealangle : float
        X angle (arcsec) in the ideal frame
    y_idealangle : float
        Y angle (arcsec) in the ideal frame

    Returns
    -------
    x_dhas : float
        X angle (arcsec) in the DHAS frame
    y_dhas : float
        Y angle (arcsec) in the DHAS frame
    """

    x_dhas = -x_idealangle
    y_dhas = y_idealangle

    return x_dhas, y_dhas


def raw2dhas(x_raw, y_raw, guider):
    """Pass in undistorted X and Y pixels in the raw/native coordinate frame
    and get out X and Y angles in the frame DHAS requires.

    Parameters
    ----------
    x_raw : float
        Undistorted X pixels in the raw coordinate frame
    y_raw : float
        Undistorted Y pixels in the raw coordinate frame
    guider : int
        Which FGS detector to convert for (1 or 2)

    Returns
    -------
    x_dhas : float
        X angle (arcsec) in the DHAS frame
    y_dhas : float
        Y angle (arcsec) in the DHAS frame
    """
    # +1 added to raw value in raw2idl to go from python 0-based to pixel coord 1-based
    x_idealangle, y_idealangle = raw2idl(x_raw, y_raw, guider)
    x_dhas, y_dhas = idl2dhas(x_idealangle, y_idealangle)

    return x_dhas, y_dhas


def raw2sci(x_raw, y_raw, guider):
    """Pass in undistorted X and Y pixels in the raw/native coordinate frame
    and get out undistorted X and Y pixels in the sci coordinate frame

    Parameters
    ----------
    x_raw : float
        Undistorted X pixels in the raw coordinate frame
    y_raw : float
        Undistorted Y pixels in the raw coordinate frame
    guider : int
        Which FGS detector to convert for (1 or 2)

    Returns
    -------
    x_sci : float
        Undistorted X pixels in the sci coordinate frame
    y_sci : float
        Undistorted Y pixels in the sci coordinate frame
    """
    if int(guider) == 1:
        fgs_full = FGS_SIAF['FGS1_FULL']
    elif int(guider) == 2:
        fgs_full = FGS_SIAF['FGS2_FULL']
    else:
        raise ValueError('Unrecognized guider number: {}'.format(guider))

    # Add +1 to incoming values to convert from python array 0-indexing to pixel coordinates 1-indexing
    x_raw_px = x_raw + 1
    y_raw_px = y_raw + 1

    x_sci, y_sci = fgs_full.raw_to_sci(x_raw_px, y_raw_px)

    # Remove 1 from outgoing values to convert back to python array 0-indexing
    x_sci -= 1
    y_sci -= 1

    return x_sci, y_sci


def convert_sky_to_idl(gs_ra, gs_dec, pa, ra_list, dec_list, guider, oss=False):
    """
    Convert
    """
    detector = f'FGS{guider}_FULL_OSS' if oss else f'FGS{guider}_FULL'
    fgs = pysiaf.Siaf('FGS')[detector]

    # Guide segment information
    gs_ra *= u.deg
    gs_dec *= u.deg
    pa *= u.deg

    # The v2/v3 location of the guide segment on the detector - in the middle
    v2 = fgs.V2Ref * u.arcsec
    v3 = fgs.V3Ref * u.arcsec

    attitude = pysiaf.rotations.attitude_matrix(v2, v3, gs_ra, gs_dec, pa)
    fgs.set_attitude_matrix(attitude)

    idl_x, idl_y = fgs.sky_to_idl(ra_list * u.deg, dec_list * u.deg)

    return idl_x.round(decimals=7), idl_y.round(decimals=7)


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
