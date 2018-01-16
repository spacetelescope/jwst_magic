import jwxml
import glob

# LOCAL
try:
    import log
except ImportError:
    in_tool = False

'''
Convert between FGS raw/native frame (pixels), ideal angle frame (arcsec), and
DHAS frame (arcsec).
'''

# Find most recent SIAF .xml file
fgs_siaf_dir = '***REMOVED***/share/SIAF_WG/Instruments/FGS/'
all_xmls = glob.glob(fgs_siaf_dir + '*.xml')
all_xmls.sort
siaf_filename = all_xmls[-1]

# Open with JWXML
fgs_siaf = jwxml.SIAF(filename=siaf_filename)

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
        fgs_full = fgs_siaf.apertures['FGS1_FULL']
    elif guider == 2:
        fgs_full = fgs_siaf.apertures['FGS2_FULL']

    # If invalid guider number provided...
    elif in_tool:
        log.error('Unrecognized guider number: {}'.format(guider))
    else:
        raise ValueError('Unrecognized guider number: {}'.format(guider))

    # Convert detector frame to ideal frame
    x_idealangle, y_idealangle = fgs_full.Det2Idl(x_det, y_det)

    return x_idealangle, y_idealangle


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
