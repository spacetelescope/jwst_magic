"""Convert between FGS raw/native frame (pixels), ideal angle frame
(arcsec), and DHAS frame (arcsec).

Authors
-------
    - Lauren Chambers

Use
---
    ::
        from jwst_magic.coordinate_transforms import Raw2Det
        x_det, y_det = Raw2Det(x_raw, y_raw)
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


def Raw2Det(x_raw, y_raw):
    """Pass in X Y pixels in the raw/native frame and get out X Y pixels in the
    SIAF detector frame

    Parameters
    ----------
    x_raw : float
        X pixels in the raw detector frame
    y_raw : float
        Y pixels in the raw detector frame

    Returns
    -------
    x_det : float
        X pixels in the SIAF detector frame
    y_det : float
        Y pixels in the SIAF detector frame
    """
    # Flip raw axes to get det axes
    x_det = y_raw
    y_det = x_raw
    return x_det, y_det


def Raw2Idl(x_raw, y_raw, guider):
    """Pass in X Y pixels in the raw/native frame and get out X Y angles in the
    ideal frame

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
    # Flip raw axes to get det axes
    x_det, y_det = Raw2Det(x_raw, y_raw)

    if int(guider) == 1:
        fgs_full = FGS_SIAF['FGS1_FULL']
    elif int(guider) == 2:
        fgs_full = FGS_SIAF['FGS2_FULL']

    # If invalid guider number provided...
    else:
        raise ValueError('Unrecognized guider number: {}'.format(guider))

    # Convert detector frame to ideal frame
    x_idealangle, y_idealangle = fgs_full.det_to_idl(x_det, y_det)

    return x_idealangle, y_idealangle


def Raw2Tel(x_raw, y_raw, guider):
    """Pass in X Y pixels in the raw/native frame and get out angles in the
    V2/V3 frame

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
    # Flip raw axes to get det axes
    x_det, y_det = Raw2Det(x_raw, y_raw)

    if int(guider) == 1:
        fgs_full = FGS_SIAF['FGS1_FULL']
    elif int(guider) == 2:
        fgs_full = FGS_SIAF['FGS2_FULL']

    # If invalid guider number provided...
    else:
        raise ValueError('Unrecognized guider number: {}'.format(guider))

    # Convert detector frame to V2/V3 frame
    v2, v3 = fgs_full.det_to_tel(x_det, y_det)

    # Subtract V2 and V3 references
    v2ref = fgs_full.V2Ref
    v3ref = fgs_full.V3Ref
    v2 -= v2ref  # this is putting the JWST FOV TEL frame, to have the origin on the FGS v2/v3 reference point
    v3 -= v3ref

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

def Idl2Raw(x_idealangle, y_idealangle, guider):
    """Pass in X and Y angles in the ideal frame and get out X Y pixels in the
    raw detector frame.

    Parameters
    ----------
    x_idealangle : float
        X angle (arcsec) in the ideal frame
    y_idealangle : float
        Y angle (arcsec) in the ideal frame
    guider : int
        Which FGS detector to convert for (1 or 2)

    Returns
    -------
    x_raw : float
        X pixels in the raw detector frame
    y_raw : float
        Y pixels in the raw detector frame
    """

    x_dhas, y_dhas = Idl2DHAS(x_idealangle, y_idealangle)
    x_raw, y_raw = DHAS2Raw(x_dhas, y_dhas, guider)

    return x_raw, y_raw


def Raw2DHAS(x_raw, y_raw, guider):
    """Pass in X and Y pixels in the raw/native frame and get out X and Y angles in
    the frame DHAS requires.

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
    x_dhas : float
        X angle (arcsec) in the DHAS frame
    y_dhas : float
        Y angle (arcsec) in the DHAS frame
    """
    x_idealangle, y_idealangle = Raw2Idl(x_raw, y_raw, guider)
    x_dhas, y_dhas = Idl2DHAS(x_idealangle, y_idealangle)

    return x_dhas, y_dhas


def DHAS2Raw(x_dhas, y_dhas, guider):
    """Pass in X and Y angles in the DHAS frame and get out X Y pixels in the
    raw detector frame

    Parameters
    ----------
    x_dhas : float
        X angle (arcsec) in the DHAS frame
    y_dhas : float
        Y angle (arcsec) in the DHAS frame
    guider : int
        Which FGS detector to convert for (1 or 2)

    Returns
    -------
    x_raw : float
        X pixels in the raw detector frame
    y_raw : float
        Y pixels in the raw detector frame
    """
    if int(guider) == 1:
        fgs_full = FGS_SIAF['FGS1_FULL']
    elif int(guider) == 2:
        fgs_full = FGS_SIAF['FGS2_FULL']

    # If invalid guider number provided...
    else:
        raise ValueError('Unrecognized guider number: {}'.format(guider))

    # DHAS to ideal
    x_idealangle = -x_dhas
    y_idealangle = y_dhas

    # Ideal to detector
    x_det, y_det = fgs_full.idl_to_det(x_idealangle, y_idealangle)

    # Detector to raw
    try:
        x_raw = y_det.astype(int)
        y_raw = x_det.astype(int)
    except AttributeError:
        x_raw = int(y_det)
        y_raw = int(x_det)

    return x_raw, y_raw


def write_to_file(xangle, yangle):
    """ Write ideal angles to file

    Parameters
    ----------
    xangle : float
        X angle (arcsec)
    yangle : float
        Y angle (arcsec)
    """
    with open('ideal.tmp', 'w') as f:
        f.write('index xangle yangle\n')
        try:
            len(xangle)
            for i, (xang, yang) in enumerate(zip(xangle, yangle)):
                f.write('{} {} {}\n'.format(i, xang, yang))
        except TypeError:
            f.write('0 {} {}\n'.format(xangle, yangle))
