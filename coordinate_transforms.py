'''
Convert between real pixel, ideal angle, and DHAS angle coordinates.
'''
import jwxml
import glob

#LOCAL
try:
    import log
except ImportError:
    in_tool = False

# Find most recent SIAF .xml file
fgs_siaf_dir = '***REMOVED***/share/SIAF_WG/Instruments/FGS/'
all_xmls = glob.glob(fgs_siaf_dir + '*.xml')
all_xmls.sort
siaf_filename = all_xmls[-1]

# Open with JWXML
fgs_siaf = jwxml.SIAF(filename=siaf_filename)

def rptoia(x_realpixel, y_realpixel, guider):
    '''
    Pass in X Y in real pixels and get out X Y in ideal angle
    '''

    if guider == 1:
        fgs_full = fgs_siaf.apertures['FGS1_FULL']
    elif guider == 2:
        fgs_full = fgs_siaf.apertures['FGS2_FULL']
    elif in_tool:
        log.error('Unrecognized guider number: {}'.format(guider))
    else:
        raise ValueError('Unrecognized guider number: {}'.format(guider))

    x_idealangle, y_idealangle = fgs_full.Det2Idl(x_realpixel, y_realpixel)

    return x_idealangle, y_idealangle

def iatoDHAS(x_idealangle, y_idealangle, guider):
    '''
    Pass in X and Y in the ideal angle frame and get out X and Y in the frame
    DHAS requires.
    '''

    if guider == 1:
        # Reverse ideal y
        y_idealangle = -y_idealangle
    elif guider == 2:
        # Reverse ideal x
        x_idealangle = -x_idealangle
    elif in_tool:
        log.error('Unrecognized guider number: {}'.format(guider))
    else:
        raise ValueError('Unrecognized guider number: {}'.format(guider))

    # Flip X and Y axes
    x_dhas, y_dhas = y_idealangle, x_idealangle

    return x_dhas, y_dhas

def iatorp(x_idealangle, y_idealangle, guider):
    '''
    Pass in X Y in ideal angle and get out X Y in real pixel
    '''

    if guider == 1:
        fgs_full = fgs_siaf.apertures['FGS1_FULL']
    elif guider == 2:
        fgs_full = fgs_siaf.apertures['FGS2_FULL']
    elif in_tool:
        log.error('Unrecognized guider number: {}'.format(guider))
    else:
        raise ValueError('Unrecognized guider number: {}'.format(guider))

    x_realpixel, y_realpixel = fgs_full.Idl2Det(x_idealangle, y_idealangle)

    return x_realpixel, y_realpixel

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
