#STDLIB
import os

#THIRD PARTY
from astropy.io import fits
import numpy as np

#LOCAL
from jwst_fgs_commissioning_tools import log, utils, coordinate_transforms
from jwst_fgs_commissioning_tools.fsw_file_writer import mkproc

FSW_PATH = os.path.dirname(os.path.realpath(__file__))
PACKAGE_PATH = os.path.split(FSW_PATH)[0]
OUT_PATH = os.path.split(PACKAGE_PATH)[0]  # Location of out/ and logs/ directory
DATA_PATH = os.path.join(PACKAGE_PATH, 'data')

def write_all(obj):
    '''
    Create **all** the files and images needs for ID & ACQ
    Requires: time_normed_im, step, x, y, countrate, bias, idarr

    ID & ACQ:
    <name>_G<guider>_<step>sky.fits :   fits file of they time_normed_im array

    <name>_G<guider>_<step>.stc :       list of x, y, and countrate for the guide
                                        and reference stars
    <name>_G<guider>_<step>bias.fits :  fits file of the bias cube for ID

    <name>_G<guider>_<step>cds.fits :

    ID Only:
    <name>_G<guider>_<step>.gssscat :   list of x and y coordinates of the guide
                                        and reference stars
    <name>_G<guider><step>ff.fits :     fits cube of the sky with added bias and noise
                                        to simulate the read and ramp cycle for ID
    <name>_G<guider>_<step>strips.fits : fits cube of overlaping strips cut
                                         from arr

    ACQ Only:
    <name>_G<guider>_<step>.cat :   list of x and y coordinates of the guide
                                        and reference stars
    <name>_G<guider><step>.fits :   fits cube of the sky with added bias and noise
                                        to simulate the read and ramp cycle for ID
    '''
    ## STScI only files - mostly just for quick checks of the data
    # Sky imge
    filename_sky = os.path.join(obj.out_dir,
                                'stsci',
                                '{}_G{}_{}sky.fits'.format(obj.root,
                                                           obj.guider,
                                                           obj.step))
    utils.write_fits(filename_sky, obj.time_normed_im)

    # STC files using offset, rotated catalog
    filename_stc = os.path.join(obj.out_dir,
                                'stsci',
                                '{}_G{}_{}.stc'.format(obj.root,
                                                       obj.guider,
                                                       obj.step))
    if obj.step == 'ACQ1' or obj.step == 'ACQ2':
        write_stc(filename_stc, obj.xarr - obj.imgsize // 2,
                  obj.yarr - obj.imgsize // 2, obj.countrate,
                  obj.guider)
    else:
        write_stc(filename_stc, obj.xarr, obj.yarr,
                  obj.countrate, obj.guider)

    # Bias image
    if obj.bias is not None:
        filename_bias = os.path.join(obj.out_dir,
                                     'stsci',
                                     '{}_G{}_{}bias.fits'.format(obj.root,
                                                                 obj.guider,
                                                                 obj.step))
        utils.write_fits(filename_bias, obj.bias)

    # Create CDS image
    if obj.cds is not None:
        filename_cds = os.path.join(obj.out_dir,
                                    'stsci',
                                    '{}_G{}_{}cds.fits'.format(obj.root,
                                                               obj.guider,
                                                               obj.step))
        utils.write_fits(filename_cds, obj.cds)

    if obj.step == 'ID':
        ## STScI only files - mostly just for quick checks of the data
        # Star catalog in real pixs
        filename_starcat = os.path.join(obj.out_dir,
                                        'stsci',
                                        '{}_G{}_{}.gssscat'.format(obj.root,
                                                                   obj.guider,
                                                                   obj.step))
        write_cat(filename_starcat, obj.xarr, obj.yarr)

        # Create ff fits file
        filename_ff = os.path.join(obj.out_dir,
                                   'stsci',
                                   '{}_G{}_{}ff.fits'.format(obj.root,
                                                             obj.guider,
                                                             obj.step))
        utils.write_fits(filename_ff, np.uint16(obj.image))

        ## DHAS file
        # Extract strips from ff img
        filename_id_strips = os.path.join(obj.out_dir,
                                          'dhas',
                                          '{}_G{}_{}strips.fits'.format(obj.root,
                                                                        obj.guider,
                                                                        obj.step))
        # Write to strips to fits file
        filename_hdr = os.path.join(DATA_PATH,
                                    'newG{}magicHdrImg.fits'.format(obj.guider))
        hdr0 = fits.getheader(filename_hdr, ext=0)
        utils.write_fits(filename_id_strips, obj.strips, header=hdr0)
        mkproc.Mkproc(obj.guider, obj.root, obj.xarr, obj.yarr, obj.countrate,
                      step='ID', out_dir=obj.out_dir)

        ## Ground system file
        convert_fits_to_dat(filename_id_strips, obj.step,
                            os.path.join(obj.out_dir, 'ground_system'))

    elif obj.step == 'ACQ1' or obj.step == 'ACQ2':
        ## STScI only files - mostly just for quick checks of the data
        # star catalog in real pixs
        filename_starcat = os.path.join(obj.out_dir,
                                        'stsci',
                                        '{}_G{}_{}.cat'.format(obj.root,
                                                               obj.guider,
                                                               obj.step))
        write_cat(filename_starcat, obj.xarr, obj.yarr)

        ## DHAS file
        # Noisy sky acquisition fits images
        filename_noisy_sky = os.path.join(obj.out_dir,
                                          'dhas',
                                          '{}_G{}_{}.fits'.format(obj.root,
                                                                  obj.guider,
                                                                  obj.step))

        utils.write_fits(filename_noisy_sky, np.uint16(obj.image))

        if obj.step == 'ACQ1':
            mkproc.Mkproc(obj.guider, obj.root, obj.xarr, obj.yarr, obj.countrate,
                          step='ACQ', out_dir=obj.out_dir,
                          acq1_imgsize=obj.acq1_imgsize,
                          acq2_imgsize=obj.acq2_imgsize)

        ## Ground system file
        convert_fits_to_dat(filename_noisy_sky, obj.step,
                            os.path.join(obj.out_dir, 'ground_system'))

    else:
        ## DHAS file
        # Noisy sky acquisition fits images
        filename_image = os.path.join(obj.out_dir,
                                      'dhas',
                                      '{}_G{}_{}.fits'.format(obj.root,
                                                              obj.guider,
                                                              obj.step))

        utils.write_fits(filename_image, np.uint16(obj.image))
        ## Gound system file
        convert_fits_to_dat(filename_image, obj.step,
                            os.path.join(obj.out_dir, 'ground_system'))

#-------------------------------------------------------------------------------
def convert_fits_to_dat(infile, obsmode, out_dir, root=None):
    '''
    Convert a .fits file to a .dat file for use on the ground system

    If 'infile' is an array, provide 'root'.

    Parameters
    ----------
    infile: str, array-like
        Can be a str (implies that this is a .fits file) or an array/cube
    obsmode: str
        The mode of image (i.e. 'PSF', 'CAL', 'TRK', 'ACQ'/'ACQ1'/'ACQ2', or 'ID')
    outfile: str
        Where to save the file
    root: str
        If infile is array-like, please provide the root image name
    '''

    obsmode = obsmode.upper()

    if isinstance(infile, str):
        data = fits.getdata(infile)

        filename = infile.split('/')[-1]
        root = filename.split('.')[0]
    else:
        root = '{}_{}'.format(root, obsmode)

    outfile = '{}.dat'.format(root)
    #data = swap_if_little_endian(data)
    flat = data.flatten()

    if (obsmode == 'PSF') or (obsmode == 'TRK') or (obsmode == 'LOSTRK'):
        # ascii float format
        fmt = '{:16.7e} '

    elif (obsmode == 'ID') or (obsmode == 'ACQ1') or (obsmode == 'ACQ2') or \
         (obsmode == 'ACQ') or (obsmode == 'CAL'):
        # ascii hex dat format
        fmt = '{:04X} '

    else:
        log.error("FSW File Writing: Observation mode not recognized. Returning.")

    with open(os.path.join(out_dir, outfile), 'w') as file_out:
        for dat in flat.astype(np.uint16):
            file_out.write(fmt.format(dat))

    print("Successfully wrote: {}".format(os.path.join(out_dir, outfile)))
    return

### Write out files
def write_stc(filename, xarr, yarr, countrate, guider):
    """
    Write out stc files using offset, rotated catalog
    """
    xia, yia = coordinate_transforms.Raw2DHAS(xarr, yarr, guider)
    xia = np.asarray(xia)
    yia = np.asarray(yia)
    countrate = np.asarray(countrate)

    inum = np.arange(len(xia))
    data = np.array([inum, xia, yia, countrate]).T

    utils.write_to_file(filename, data, fmt=['%d', '%f', '%f', '%e'])
    print("Successfully wrote: {}".format(filename))


def write_cat(filename, xarr, yarr):
    '''
    Write out star catalog in real pixs
    '''
    coords = np.array([xarr, yarr]).T
    try:
        utils.write_to_file(filename, coords, fmt=['%d', '%d'])
    except AttributeError:
        utils.write_to_file(filename, coords)
    print("Successfully wrote: {}".format(filename))
