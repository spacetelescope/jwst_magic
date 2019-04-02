"""Writes all flight software files for CAL, ID, ACQ, and/or TRK steps

This module takes an FGS simulation object for CAL, ID, ACQ, and/or TRK
steps (created in ``buildfgssteps.py``), which it uses to write the
necessary flight software files for use with the DHAS, the FGSES/
CertLab, or just for user inspection. The files created for each step
are as follows:

    CAL:
        sky.fits (the input file after being time-normalized (converted
            from counts/s to counts)
        bias.fits (bias file used to create noisy FGS image)
        cds.fits (correlated double sample)
        .fits
    ID:
        sky.fits (the input file after being time-normalized (converted
            from counts/s to counts)
        bias.fits (bias file used to create noisy FGS image)
        cds.fits (correlated double sample)
        ff.fits (full frame; image before CDS)
        strips.fits (strips to run in DHAS)
        strips.dat (strips to run in FGSES)
        .gssscat
        .stc
        .prc (to run in DHAS or FGSES)
    ACQ1 or ACQ2:
        sky.fits (the input file after being time-normalized (converted
            from counts/s to counts)
        bias.fits (bias file used to create noisy FGS image)
        cds.fits (correlated double sample)
        .fits (to run in DHAS)
        .dat (to run in FGSES)
        .cat
        .stc
        .prc (to run in DHAS or FGSES)
    LOSTRK:
        .fits
        .dat (to run in FGSES)
    TRK:
        .fits (to run in DHAS)

Authors
-------
    - Keira Brooks
    - Lauren Chambers

Use
---
    This module can be executed in a Python shell as such:
    ::
        from jwst_magic.fsw_file_writer import write_files
        write_files.write_all(obj)

    Required arguments:
        ``obj`` - FGS simulation object for CAL, ID, ACQ, and/or TRK
            stages; created by ``buildfgssteps.py``
"""

# Standard Library Imports
import os
import logging

# Third Party Imports
from astropy.io import fits
import numpy as np

# Local Imports
from .. import utils
from ..fsw_file_writer import mkproc
from jwst_magic._utils import coordinate_transforms

# Paths
__location__ = os.path.realpath(os.path.join(os.getcwd(), os.path.dirname(__file__)))
PACKAGE_PATH = os.path.split(__location__)[0]
OUT_PATH = os.path.split(PACKAGE_PATH)[0]  # Location of out/ and logs/ directory
DATA_PATH = os.path.join(PACKAGE_PATH, 'data')

# Start logger
LOGGER = logging.getLogger(__name__)


def write_all(obj):
    """Create **all** the files and images needed for simulation with
    the flight software and use by analysts for a particular step

    Parameters
    ----------
    obj : obj
        FGS simulation object for CAL, ID, ACQ, and/or TRK stages;
        created by ``buildfgssteps.py``
    """
    # Determine the root filename
    filename_root = '{}_G{}_{}'.format(obj.root, obj.guider, obj.step)
    obj.filename_root = filename_root

    if obj.step == 'CAL':
        # Write files for use by folks at STScI
        write_sky(obj)
        write_bias(obj)
        write_cds(obj)
        write_image(obj)

    elif obj.step == 'ID':
        # Write files for use by folks at STScI
        write_sky(obj)
        write_bias(obj)
        write_stc(obj)
        write_cds(obj)
        write_cat(obj)
        write_image(obj)

        # Write files for use in the DHAS and FGSES
        write_strips(obj)
        write_prc(obj)
        write_dat(obj)

    elif obj.step == 'ACQ1':
        # Write files for use by folks at STScI
        write_sky(obj)
        write_bias(obj)
        write_stc(obj)
        write_cds(obj)
        write_cat(obj)

        # Write files for use in the DHAS and FGSES
        write_image(obj)
        write_prc(obj)
        write_dat(obj)

    elif obj.step == 'ACQ2':
        # Write files for use by folks at STScI
        write_sky(obj)
        write_bias(obj)
        write_stc(obj)
        write_cds(obj)
        write_cat(obj)

        # Write files for use in the DHAS and FGSES
        write_image(obj)
        write_dat(obj)

    elif obj.step == 'TRK':
        # Write files for use by folks at STScI
        write_sky(obj)
        write_bias(obj)
        write_stc(obj)

        # Write files for use in the DHAS
        write_image(obj)

    elif obj.step == 'LOSTRK':
        # Write files for use by folks at STScI
        write_sky(obj)
        write_stc(obj)
        write_image(obj)

        # Write files for use in the FGSES
        write_dat(obj)


def write_sky(obj):
    """Write the time-normed image, or "sky" image

    Units:  counts
    Size:   n_cols x n_rows

    Parameters
    ----------
    obj : obj
        FGS simulation object for CAL, ID, ACQ, and/or TRK stages;
        created by ``buildfgssteps.py``
    """
    filename_sky = os.path.join(obj.out_dir, 'stsci',
                                obj.filename_root + 'sky.fits')
    utils.write_fits(filename_sky, np.uint16(obj.time_normed_im), log=LOGGER)


def write_bias(obj):
    """Write the bias image

    Units:  counts
    Size:   n_cols x n_rows x (n_reads x n_ramps)

    Parameters
    ----------
    obj : obj
        FGS simulation object for CAL, ID, ACQ, and/or TRK stages;
        created by ``buildfgssteps.py``
    """
    if obj.bias is not None:
        filename_bias = os.path.join(obj.out_dir,
                                     'stsci',
                                     obj.filename_root + 'bias.fits')
        utils.write_fits(filename_bias, np.uint16(obj.bias), log=LOGGER)


def write_cds(obj):
    """Write the correlated double sample (CDS) image by subtracting
    the 0th read from the 1st read

    Units:  counts
    Size:   n_cols x n_rows x ((n_reads x n_ramps)/2)

    Parameters
    ----------
    obj : obj
        FGS simulation object for CAL, ID, ACQ, and/or TRK stages;
        created by ``buildfgssteps.py``
    """
    if obj.cds is not None:
        filename_cds = os.path.join(obj.out_dir,
                                    'stsci',
                                    obj.filename_root + 'cds.fits')
        utils.write_fits(filename_cds, np.uint16(obj.cds), log=LOGGER)


def write_image(obj):
    """Write a normal image (i.e. ACQ1, TRK, CAL)

    Units:  counts
    Size:   n_cols x n_rows x (n_reads x n_ramps)

    Parameters
    ----------
    obj : obj
        FGS simulation object for CAL, ID, ACQ, and/or TRK stages;
        created by ``buildfgssteps.py``
    """

    if obj.step == 'ID':
        # Create "full-frame" (rather than strips) image
        location = 'stsci'
        filetype = 'ff.fits'
    elif obj.step == 'LOSTRK':
        # Place the FITS file in stsci/, as it is just for reference
        # to the LOSTRK.dat file
        location = 'stsci'
        filetype = '.fits'
    else:
        location = 'dhas'
        filetype = '.fits'

    # Create image fits file
    filename = os.path.join(obj.out_dir, location,
                            obj.filename_root + filetype)

    if obj.step == 'LOSTRK':
        utils.write_fits(filename, obj.image, log=LOGGER) # Don't make it np.uint16
    else:
        utils.write_fits(filename, np.uint16(obj.image), log=LOGGER)

def write_strips(obj):
    """Write an ID strips image

    Units:  counts
    Size:   2048 x 64 x (36 x n_reads x n_ramps)

    Parameters
    ----------
    obj : obj
        FGS simulation object for CAL, ID, ACQ, and/or TRK stages;
        created by ``buildfgssteps.py``
    """

    # Extract strips from ff img
    filename_id_strips = os.path.join(obj.out_dir,
                                      'dhas',
                                      obj.filename_root + 'strips.fits')
    # Write to strips to fits file
    filename_hdr = os.path.join(DATA_PATH,
                                'newG{}magicHdrImg.fits'.format(obj.guider))
    hdr0 = fits.getheader(filename_hdr, ext=0)
    utils.write_fits(filename_id_strips, np.uint16(obj.strips), header=hdr0,
                     log=LOGGER)


def write_prc(obj):
    """Write a procedure (.prc) file for use with file software

    Parameters
    ----------
    obj : obj
        FGS simulation object for CAL, ID, ACQ, and/or TRK stages;
        created by ``buildfgssteps.py``

    Notes
    -----
        Don't need to include the yoffset in the .prc IF the strips
        offset parameter in DHAS is set to 12
    """
    if obj.step == 'ID':
        step = 'ID'
        acq1_imgsize = None
        acq2_imgsize = None

    elif obj.step == 'ACQ1':
        step = 'ACQ'
        acq1_imgsize = obj.acq1_imgsize
        acq2_imgsize = obj.acq2_imgsize

    else:
        raise ValueError('Unknown step {}; cannot write .prc file.'.format(obj.step))

    mkproc.Mkproc(obj.guider, obj.root, obj.xarr, obj.yarr, obj.countrate,
                  step=step, out_dir=obj.out_dir, acq1_imgsize=acq1_imgsize,
                  acq2_imgsize=acq2_imgsize)


def write_dat(obj):
    """
    Convert a .fits file to a .dat file for use on the ground system

    Parameters
    ----------
    obj : obj
        FGS simulation object for CAL, ID, ACQ, and/or TRK stages;
        created by ``buildfgssteps.py``
    """
    if obj.step == 'ID':
        data_to_write = obj.strips
    else:
        data_to_write = obj.image

    out_dir = os.path.join(obj.out_dir, 'ground_system')

    obsmode = obj.step.upper()

    if isinstance(data_to_write, str):
        data = fits.getdata(data_to_write)
        filename = data.split('/')[-1].split('.')[0]
    else:
        data = data_to_write
        filename = '{}_G{}_{}.dat'.format(obj.root, obj.guider, obsmode)
    filename = os.path.join(out_dir, filename)

    flat = data.flatten()

    # Write out TRK/LOSTRK files in ASCII float format
    if (obsmode == 'PSF') or (obsmode == 'TRK') or (obsmode == 'LOSTRK'):
        # ascii float format
        fmt = '{:16.7e} '
        with open(filename, 'w') as file_out:
            for dat in flat:  # Note: NOT saving out as uint16!!!
                file_out.write(fmt.format(dat))

    # Write out all other files in ASCII hex format (from uint16)
    elif (obsmode == 'ID') or (obsmode == 'ACQ1') or (obsmode == 'ACQ2') or \
         (obsmode == 'ACQ') or (obsmode == 'CAL'):
        # ascii hex dat format
        fmt = '{:04X} '
        with open(filename, 'w') as file_out:
            for dat in flat.astype(np.uint16):
                file_out.write(fmt.format(dat))

    else:
        raise ValueError("FSW File Writing: Observation mode {} not recognized.".format(obsmode))

    LOGGER.info("Successfully wrote: {}".format(filename))
    return


def write_stc(obj):
    """
    Write out stc files using offset, rotated catalog

    Parameters
    ----------
    obj : obj
        FGS simulation object for CAL, ID, ACQ, and/or TRK stages;
        created by ``buildfgssteps.py``
    """
    # STC files using offset, rotated catalog
    filename_stc = os.path.join(obj.out_dir,
                                'stsci',
                                obj.filename_root + '.stc')
    if obj.step == 'ACQ1' or obj.step == 'ACQ2':
        xarr = obj.xarr - obj.imgsize // 2
        yarr = obj.yarr - obj.imgsize // 2
    else:
        xarr = obj.xarr
        yarr = obj.yarr

    xia, yia = coordinate_transforms.Raw2DHAS(xarr, yarr, obj.guider)
    xia = np.asarray(xia)
    yia = np.asarray(yia)
    countrate = np.asarray(obj.countrate)

    inum = np.arange(len(xia))
    data = np.array([inum, xia, yia, countrate]).T

    utils.write_to_file(filename_stc, data, fmt=['%d', '%f', '%f', '%e'])
    LOGGER.info("Successfully wrote: {}".format(filename_stc))


def write_cat(obj):
    """
    Write out star catalog in real pixels

    Parameters
    ----------
    obj : obj
        FGS simulation object for CAL, ID, ACQ, and/or TRK stages;
        created by ``buildfgssteps.py``
    """
    if obj.step == 'ID':
        filetype = '.gssscat'
    else:
        filetype = '.cat'

    filename_starcat = os.path.join(obj.out_dir, 'stsci',
                                    obj.filename_root + filetype)

    coords = np.array([obj.xarr, obj.yarr]).T
    try:
        utils.write_to_file(filename_starcat, coords, fmt=['%d', '%d'])
    except AttributeError:
        utils.write_to_file(filename_starcat, coords)
    LOGGER.info("Successfully wrote: {}".format(filename_starcat))
