"""Convert between FGS raw/native frame (pixels), ideal angle frame
(arcsec), and DHAS frame (arcsec).

Note that the images used in MAGIC are all undistorted, so whenever
the "raw", "det" or "sci" frames are discussed here, they are not the
actual frames from the SIAF definition, but are instead an undistorted
frame where the coodinte system matches the "raw", "det" or "sci" frames
from the SIAF respectively.

Authors
-------
    - Lauren Chambers
    - Shannon Osborne

Use
---
    ::
        from jwst_magic.coordinate_transforms import Raw2Idl
        x_idl, y_idl = Raw2Idl(x_raw, y_raw)
"""

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
    nrc_det = nrc_siaf[f'{detector}_FULL_OSS']
    nircam_x_scale = nrc_det.XSciScale  # arcsec/pixel
    nircam_y_scale = nrc_det.YSciScale  # arcsec/pixel

    # Convert x/y offsets to V2/V3
    v2_offset = x_offset * nircam_x_scale  # arcsec
    v3_offset = y_offset * nircam_y_scale  # arcsec

    return v2_offset, v3_offset


def Raw2Idl(x_raw, y_raw, guider):
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

    # Convert from RAW -> SCI (just a coordinate origin change, no distortion change)
    x_sci, y_sci = fgs_full.raw_to_sci(x_raw, y_raw)

    # Shift to the IDL origin location
    x_idl_pix = x_sci - fgs_full.XSciRef  # TODO TBD ON IF USING THIS VALUE IS OKAY
    y_idl_pix = y_sci - fgs_full.YSciRef

    # Convert from IDL pixels to idl arcseconds
    x_idealangle = x_idl_pix * fgs_full.XSciScale
    y_idealangle = y_idl_pix * fgs_full.YSciScale

    return x_idealangle, y_idealangle


def Raw2Tel(x_raw, y_raw, guider):
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

    # Convert from Raw to IDL
    x_idl, y_idl = Raw2Idl(x_raw, y_raw, guider)

    # Convert from IDL -> TEL (just a coordinate origin change, no distortion change)
    v2, v3 = fgs_full.idl_to_tel(x_idl, y_idl)

    return v2, v3


def Idl2DHAS(x_idealangle, y_idealangle):
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


def Raw2DHAS(x_raw, y_raw, guider):
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
    x_idealangle, y_idealangle = Raw2Idl(x_raw, y_raw, guider)
    x_dhas, y_dhas = Idl2DHAS(x_idealangle, y_idealangle)

    return x_dhas, y_dhas
