import pysiaf

# LOCAL
try:
    import log
except ImportError:
    in_tool = False

'''
Convert between FGS raw/native frame (pixels), ideal angle frame (arcsec), and
DHAS frame (arcsec).
'''

# Open SIAF with pysiaf
fgs_siaf = pysiaf.Siaf('FGS')

def Raw2Det(x_raw, y_raw):
    '''
    Pass in X Y pixels in the raw/native frame and get out X Y pixels in the
    SIAF detector frame
    '''
    # Flip raw axes to get det axes
    x_det = y_raw
    y_det = x_raw
    return x_det, y_det

def Raw2Idl(x_raw, y_raw, guider):
    '''
    Pass in X Y pixels in the raw/native frame and get out X Y angles in the
    ideal frame
    '''
    # Flip raw axes to get det axes
    x_det, y_det = Raw2Det(x_raw, y_raw)

    if guider == 1:
        fgs_full = fgs_siaf['FGS1_FULL']
    elif guider == 2:
        fgs_full = fgs_siaf['FGS2_FULL']

    # If invalid guider number provided...
    elif in_tool:
        log.error('Unrecognized guider number: {}'.format(guider))
    else:
        raise ValueError('Unrecognized guider number: {}'.format(guider))

    # Convert detector frame to ideal frame
    x_idealangle, y_idealangle = fgs_full.det_to_idl(x_det, y_det)

    return x_idealangle, y_idealangle

def Raw2Tel(x_raw, y_raw, guider):
    '''
    Pass in X Y pixels in the raw/native frame and get out X Y angles in the
    V2/V3 frame
    '''
    # Flip raw axes to get det axes
    x_det, y_det = Raw2Det(x_raw, y_raw)

    if int(guider) == 1:
        fgs_full = fgs_siaf['FGS1_FULL']
    elif int(guider) == 2:
        fgs_full = fgs_siaf['FGS2_FULL']

    # If invalid guider number provided...
    elif in_tool:
        log.error('Unrecognized guider number: {}'.format(guider))
    else:
        raise ValueError('Unrecognized guider number: {}'.format(guider))

    # Convert detector frame to V2/V3 frame
    v2, v3 = fgs_full.det_to_tel(x_det, y_det)

    # Subtract V2 and V3 references
    v2ref = fgs_full.V2Ref
    v3ref = fgs_full.V3Ref
    v2 -= v2ref
    v3 -= v3ref

    return v2, v3

def Idl2DHAS(x_idealangle, y_idealangle):
    '''
    Pass in X and Y angles in the ideal frame and get out X and Y angles in the
    frame DHAS requires.
    '''

    x_dhas = -x_idealangle
    y_dhas = y_idealangle

    return x_dhas, y_dhas

def Raw2DHAS(x_raw, y_raw, guider):
    '''
    Pass in X and Y pixels in the raw/native frame and get out X and Y angles in
    the frame DHAS requires.
    '''
    x_idealangle, y_idealangle = Raw2Idl(x_raw, y_raw, guider)
    x_dhas, y_dhas = Idl2DHAS(x_idealangle, y_idealangle)

    return x_dhas, y_dhas

def DHAS2Raw(x_dhas, y_dhas, guider):

    if int(guider) == 1:
        fgs_full = fgs_siaf['FGS1_FULL']
    elif int(guider) == 2:
        fgs_full = fgs_siaf['FGS2_FULL']

    # If invalid guider number provided...
    elif in_tool:
        log.error('Unrecognized guider number: {}'.format(guider))
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
    ''' Write ideal angles to file'''
    with open('ideal.tmp', 'w') as f:
        f.write('index xangle yangle\n')
        try:
            len(xangle)
            for i, (xang, yang) in enumerate(zip(xangle, yangle)):
                f.write('{} {} {}\n'.format(i, xang, yang))
        except TypeError:
            f.write('0 {} {}\n'.format(xangle, yangle))
