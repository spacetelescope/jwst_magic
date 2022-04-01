"""Setup to write all simulated files for CAL, ID, ACQ, and/or TRK steps

This module parses the ``config.ini`` file to create an FGS simulation
object for CAL, ID, ACQ, and/or TRK steps. This object should ultimately
be passed to ``write_files.py`` to create the necessary flight software files
for use with the DHAS, the FGSES/CertLab, or just for user inspection.

Authors
-------
    - Keira Brooks
    - Lauren Chambers
    - Shannon Osborne

Use
---
    This module can be executed in a Python shell as such:
    ::
        from jwst_magic.fsw_file_writer import buildfgssteps
        buildfgssteps.BuildFGSSteps(im, guider, root, step, out_dir)

    Required arguments:
        ``im`` - image array or filepath for the input FGS image
        ``guider`` - number for guider 1 or guider 2
        ``root`` - will be used to create the output directory, ./out/{root}
        ``step`` - name of guiding step for which to create images
            (expecting 'ID', 'ACQ1', 'ACQ2', 'TRK', or 'LOSTRK')
    Optional arguments:
        ``guiding_selections_file`` - file containing X/Y positions and
            countrates for all stars in an image
        ``configfile`` - file defining parameters for each guider step.
            If not defined, defaults to jwst_magic/data/config.ini
        ``out_dir`` - where output files will be saved. If not provided,
            the image(s) will be saved within the repository at
            jwst_magic/
        ``logger_passed`` - denotes if a logger object has already been
            generated.
"""

# Standard Library Imports
import logging
import os

# Third Party Imports
from astropy.io import ascii as asc
from astropy.io import fits
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import numpy as np
from scipy.ndimage import shift

# Local imports
from jwst_magic.utils import utils
from jwst_magic.fsw_file_writer import config, detector_effects
from jwst_magic.star_selector import select_psfs

# Paths
__location__ = os.path.realpath(os.path.join(os.getcwd(), os.path.dirname(__file__)))
PACKAGE_PATH = os.path.split(__location__)[0]
OUT_PATH = os.path.split(PACKAGE_PATH)[0]  # Location of out/ and logs/ directory
DATA_PATH = os.path.join(PACKAGE_PATH, 'data')

# Start logger
LOGGER = logging.getLogger(__name__)

# Global Values
OSS_TRIGGER = 474608.4  # ADU/sec
COUNTRATE_CONVERSION = 0.65
DIM_STAR_THRESHOLD_FACTOR = 0.30
BRIGHT_STAR_THRESHOLD_ADDEND = 332226
GAIN_G1 = 1.74
GAIN_G2 = 1.57

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# MAIN CLASS
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


class BuildFGSSteps(object):
    """Creates an FGS simulation object for ID, ACQ, and/or TRK stages
    to be used with DHAS.
    """
    def __init__(self, im, guider, root, step, guiding_selections_file=None, configfile=None,
                 out_dir=None, thresh_factor=0.6, logger_passed=False, psf_center_file=None,
                 shift_id_attitude=True, use_oss_defaults=False, catalog_countrate=None,
                 override_bright_guiding=False):
        """Initialize the class and call build_fgs_steps().
        """
        # Check path exists
        utils.ensure_dir_exists(out_dir)

        # Set up logger
        if not logger_passed:
            out_dir_root = out_dir.split(root)[0] + root
            utils.create_logger_from_yaml(__name__, out_dir_root=out_dir_root, root=root, level='DEBUG')

        try:
            # Initialize attributes
            self.input_im = im
            self.guider = guider
            self.root = root
            self.step = step
            self.yoffset = 12
            self.out_dir = out_dir
            self.thresh_factor = thresh_factor
            self.use_oss_defaults = use_oss_defaults
            self.catalog_countrate = catalog_countrate
            self.threshold = None
            self.override_bright_guiding = override_bright_guiding
            if 'config' in guiding_selections_file:
                self.config = guiding_selections_file.split('/')[-1].split('.txt')[0].split('config')[-1]
            else:
                self.config = 'unknown'

            # READ IN IMAGE
            if isinstance(im, str):
                data = fits.getdata(im)
                self.input_im = data
            else:
                self.input_im = im

            # Shift image to ID attitude
            self.stsci_dir = 'stsci'
            self.dhas_dir = 'dhas'
            self.ground_system_dir = 'ground_system'
            if shift_id_attitude:
                self.stsci_dir = '{}_shifted'.format(self.stsci_dir)
                self.dhas_dir = '{}_shifted'.format(self.dhas_dir)
                self.ground_system_dir = '{}_shifted'.format(self.ground_system_dir)

            # Build FGS steps
            self.build_fgs_steps(guiding_selections_file, configfile, psf_center_file)

        except Exception as e:
            LOGGER.exception(f'{repr(e)}: {e}')
            raise

    def build_fgs_steps(self, guiding_selections_file, configfile, psf_center_file=None):
        """Creates an FGS simulation object for ID, ACQ, and/or TRK stages
        to be used with DHAS.

        Parameters
        ----------
        guiding_selections_file : str
            File containing X/Y positions and countrates for all stars
            in the provided image. Will be the shifted or unshifted version
            based on the shift_id_attitude keyword
        configfile : str
            File defining parameters for each guider step. If not
            defined, defaults to jwst_magic/data/config.ini
        psf_center_file : str, optional
            Path to psf center file in order to re-center the TRK box
            to not be centered on the guiding elections PSF location,
            but on the actual center of the PSF (found with smoothing;
            tuple of (y,x)). Used when smoothing='low'.
        """
        utils.ensure_dir_exists(os.path.join(self.out_dir, self.dhas_dir))
        utils.ensure_dir_exists(os.path.join(self.out_dir, self.ground_system_dir))
        utils.ensure_dir_exists(os.path.join(self.out_dir, self.stsci_dir))

        self.get_coords_and_counts(guiding_selections_file, psf_center_file)

        section = '{}_dict'.format(self.step.lower())
        config_ini = self.build_step(section, configfile)
        self.image = self.create_img_arrays(section, config_ini)

        # Write the files
        LOGGER.info("FSW File Writing: Creating {} FSW files".format(self.step))

    def get_coords_and_counts(self, guiding_selections_file, psf_center_file=None):
        """Get coordinate information and countrates of guide star and
        reference stars.

        Parameters
        ----------
        guiding_selections_file : str
            File containing X/Y positions and countrates for all stars
            in the provided image. Will be the shifted or unshifted version
            based on the shift_id_attitude keyword
        psf_center_file : str, optional
            Path to psf center file in order to re-center the TRK box
            to not be centered on the guiding elections PSF location,
            but on the actual center of the PSF (found with smoothing;
            tuple of (y,x)). Used when smoothing='low'. Will be the
            shifted or unshifted version based on the shift_id_attitude keyword
        """
        LOGGER.info("FSW File Writing: Using {} as the guiding selections file".format(guiding_selections_file))

        if guiding_selections_file.endswith('reg'):
            self.xarr, self.yarr = np.loadtxt(guiding_selections_file)
            self.countrate = []
            for xa, ya, in zip(self.xarr, self.yarr):
                self.countrate.append(utils.get_countrate_3x3(xa, ya, self.input_im))
        else:
            self.yarr, self.xarr, self.countrate = np.loadtxt(guiding_selections_file,
                                                              delimiter=' ',
                                                              skiprows=1).T

        # Make sure the TRK box is centered on the center of the PSF, not on the brightest point
        if self.step == 'TRK' and psf_center_file is not None:
            self.yarr, self.xarr, _ = np.loadtxt(psf_center_file, delimiter=' ', skiprows=1).T

        # Cover cases where there is only one entry in the reg file
        try:
            len(self.xarr)
        except TypeError:
            self.xarr = np.asarray([self.xarr])
            self.yarr = np.asarray([self.yarr])
            self.countrate = np.asarray([self.countrate])

        # Overwrite threshold with OSS default values for the POF case
        if self.use_oss_defaults:
            if len(self.xarr) != 1:
                raise ValueError('Trying to apply OSS defaults to non-POF case (with multiple star selections)')

            if self.catalog_countrate is None:
                raise ValueError('When creating FSW files with the OSS defaults (use_oss_defaults=True), you'
                                 'must pass in the catalog_countrate of the guide star as well.')

            countrate_3x3 = self.catalog_countrate * COUNTRATE_CONVERSION
            self.countrate = np.asarray([countrate_3x3])

        self.threshold, self.thresh_factor = bright_guiding_check(np.asarray(self.countrate),
                                                                  self.thresh_factor,
                                                                  normal_ops=self.use_oss_defaults,
                                                                  override_bright_guiding=self.override_bright_guiding)

    def build_step(self, section, configfile=None):
        """Read out the parameters in the config.ini, and alter class
        attributes to build the current step accordingly.

        Parameters
        ----------
        section : str
            Name of step within the config file ('{step}_dict')
        configfile : str, optional
            File defining parameters for each guider step. If not
            defined, defaults to jwst_magic/data/config.ini
        Returns
        -------
        config_ini : obj
            Object containing all parameters from the given configfile
        """
        self.nreads = 2

        if configfile is None:
            configfile = os.path.join(DATA_PATH, 'config.ini')
        config_ini = config.load_config_ini(configfile)

        if self.step not in ['ID', 'CAL']:
            # If not ID or CAL (i.e. if making a subarray), only use the guide star
            self.xarr = np.asarray([self.xarr[0]])
            self.yarr = np.asarray([self.yarr[0]])
            self.countrate = np.asarray([self.countrate[0]])
            self.input_im = create_im_subarray(self.input_im, self.xarr,
                                               self.yarr, config_ini.getint(section, 'imgsize'))

            if self.step == 'ACQ1' or self.step == 'ACQ2':
                self.imgsize = config_ini.getint(section, 'imgsize')
                self.acq1_imgsize = config_ini.getint('acq1_dict', 'imgsize')
                self.acq2_imgsize = config_ini.getint('acq2_dict', 'imgsize')

        return config_ini

    def create_img_arrays(self, step, config_ini):
        """Create the needed image arrays for the given step.

        Possible arrays (created as class attributes):
            - time_normed_im : the "sky" image, or the time-normalized
                image (in counts)
            - bias : the FGS bias used to simulate the image; includes
                0th read bias structure and KTC "shot" noise. An array
                of size (nramps x nreads) x n_rows x n_columns
            - cds : subtracts the zeroth read from the first read
            - strips : divides a full-frame array into 36 strips of
                size 64 x 2048

        If creating LOSTRK images, the final "image" array will be
        normalized and resized for use in the FGSES

        Parameters
        ----------
        step : str
            Name of step within the config file ('{step}_dict')
        config_ini : : obj
            Object containing all parameters from the given config file

        Returns
        -------
        image : 2-D numpy array
            The final image including any necessary detector effects,
            either full-frame or the appropriately sized subarray, in counts.
        """

        # Create the time-normalized image (will be in counts, where the
        # input_im is in counts per second)
        self.time_normed_im = self.input_im * config_ini.getfloat(step, 'tframe')

        # Add the bias, and build the array of reads with noisy data
        if config_ini.getboolean(step, 'bias'):
            gain = GAIN_G1 if self.guider == 1 else GAIN_G2

            # Get the bias ramp
            nramps = config_ini.getint(step, 'nramps')
            det_eff = detector_effects.FGSDetectorEffects(
                self.guider, self.xarr, self.yarr, self.nreads, nramps,
                config_ini.getint(step, 'imgsize')
            )
            self.bias = det_eff.add_detector_effects()

            # Take the bias and add a noisy version of the input image
            # (specifically Poisson noise), adding signal over each read
            image = np.zeros_like(self.bias)
            signal_cube = [self.time_normed_im] * nramps
            positive_signal_cube = np.copy(signal_cube)
            positive_signal_cube[positive_signal_cube < 0] = 0

            # Add signal to every first read
            i_read1 = 0
            ndrop1 = int(config_ini.getfloat(step, 'ndrops1'))
            image[i_read1::self.nreads] += self.bias[i_read1::self.nreads] + np.random.poisson(
                (ndrop1 + 1) * positive_signal_cube * gain) / gain

            # Add signal to every second read
            i_read2 = 1
            ndrop2 = int(config_ini.getfloat(step, 'ndrops2'))
            image[i_read2::self.nreads] = image[i_read1::self.nreads] + np.random.poisson(
                (ndrop2 + 1) * positive_signal_cube * gain) / gain

            # Set pixels in the image to 0 when flagged as bad in the DQ array
            image = self.add_fgs_dq(image, self.nreads, nramps, det_eff.array_bounds)

        else:
            # In the case of LOSTRK, just return one frame with one
            # frame readout worth of signal
            self.bias = None
            image = self.time_normed_im

        # # Cut any pixels under zero
        image = utils.correct_image(image)

        # Create the CDS image by subtracting the first read from the second
        # read, for each ramp
        if config_ini.getboolean(step, 'cdsimg'):
            self.cds = create_cds(image)
        else:
            self.cds = None

        # If in ID, split the full-frame image into strips
        if config_ini.getboolean(step, 'stripsimg'):
            self.strips = create_strips(image,
                                        config_ini.getint(step, 'imgsize'),
                                        config_ini.getint(step, 'nstrips'),
                                        config_ini.getint(step, 'nramps'),
                                        self.nreads,
                                        config_ini.getint(step, 'height'),
                                        self.yoffset,
                                        config_ini.getint(step, 'overlap'))

        # Modify further for LOSTRK images (that will be run in FGSES)
        if self.step == 'LOSTRK':
            # Resize image array to oversample by 6 (from 43x43 to 255x255)
            image = image.repeat(6, axis=0)
            image = image.repeat(6, axis=1)
            image = image[1:-2, 1:-2]

            # Normalize to a count sum of 1000
            image = image / np.sum(image) * 1000

        return image

    def add_fgs_dq(self, image, nreads, nramps, array_bounds):
        """Set the bias to 0 where the DQ Array is flagged
        as bad (where it equals 1)

        image :
            The time-normalized image with all noise
            already applied
        nreads : int
            Number of reads in the ramp
        nramps : int
            Number of ramps in the integration
        array_bounds : list, tuple
            The subarray location of the step
        """
        dq_file = os.path.join(DATA_PATH, 'fgs_dq_G{}.fits'.format(self.guider))
        dq_arr = np.copy(fits.getdata(dq_file))

        # Cut and duplicate DQ to match shape of bias file
        xlow, xhigh, ylow, yhigh = array_bounds
        dq_arr = dq_arr[xlow:xhigh, ylow:yhigh]
        dq_arr = np.stack([dq_arr] * nramps * nreads, axis=0)
        image[dq_arr == 1] = 0

        return image

# ------------------------------------------------------------------------------
# ANCILLARY FUNCTIONS
# ------------------------------------------------------------------------------


def create_strips(image, imgsize, nstrips, nramps, nreads, strip_height, yoffset, overlap):
    """Create the ID strips fits image to be passed into DHAS.

    Parameters
    ----------
    image : 2-D numpy array
        Image data
    imgsize : int
        Dimension of full-frame image (pixels)
    nstrips : int
        Number of strips to split the image into
    nramps : int
        Number of ramps in the integration
    nreads : int
        Number of reads in the ramp
    strip_height : int
        Height of each strip (pixels)
    yoffset : int
        The offset at the bottom of the array before the first strip (pixels)
    overlap : int
        The number of pixels of overlap between strips

    Returns
    -------
    strips : 3-D numpy array
        Array of ID strips with dimensions:
        (nramps * nreads * nz, strip_height, imgsize)
    """
    nz = nramps * nreads
    strips = np.zeros((nstrips * nz, strip_height, imgsize))
    nn = 0
    for i in range(nstrips):
        for iz in range(nz):
            ylow = i * (strip_height - overlap) + yoffset
            yhigh = ylow + strip_height
            strips[nn] = image[iz, ylow:yhigh, :]
            nn += 1

    return strips


def create_cds(arr, fix_saturated_pix=True):
    """Create CDS image: Subtract the first read from the second read.
    Option to handle saturated pixels in CDS.

    Parameters
    ----------
    arr : 3-D numpy array
        Image data
    fix_saturated_pix : boolean, optional
        Apply a fix to saturated pixels in the CDS array? Default is True.

    Returns
    -------
    2-D numpy array
        CDS of image
    """
    # Second read minus first read
    first_reads = arr[:-1:2]
    second_reads = arr[1::2]
    cds_arr = second_reads - first_reads

    if fix_saturated_pix:
        # Determine which pixels are saturated in each read
        saturated_read_2 = second_reads >= 55000
        saturated_read_1 = first_reads >= 55000

        # any pixel that saturates in either read should have a CDS value of 25K assigned
        #join saturated_read_1 and saturated_read_2
        saturated = saturated_read_2 | saturated_read_1
        cds_arr[saturated] = 25000

    return cds_arr


def display(image, ind=0, vmin=None, vmax=None, xarr=None, yarr=None):
    """Plot an image array. If the array is a cube, specify the index to
    look at using 'ind'. If you want to add a scatter plot of PSF centers,
    use add_coords=True and pass in x and y (Note: This only works for full
    frame images for the time being).

    Parameters
    ----------
    image : numpy array
        Image data (2- or 3-dimensional)
    ind : int, optional
        The frame of a 3-D array to plot
    vmin : float, optional
        The minimum value of the colorbar
    vmax : float, optional
        The maximum value of the colorbar
    xarr : list, optional
        List of x values to plot as points over image
    yarr : list, optional
        List of y values to plot as points over image
    """
    if image.ndim == 3:
        img = image[ind]
    else:
        img = np.copy(image)

    plt.clf()
    plt.imshow(img, cmap='Greys_r', norm=LogNorm(), vmin=vmin, vmax=vmax)

    if xarr and yarr:
        plt.scatter(xarr, yarr, color='white')
    else:
        plt.colorbar()

    plt.show()


def create_im_subarray(image, xcoord, ycoord, imgsize, show_fig=False):
    """Based on the array size given by nx and ny, created a subarray
    around the guide star.

    Parameters
    ----------
    image : 2-D numpy array
        Image data
    xcoord : int
        X coordinate of guide star (pixels)
    ycoord : int
        Y coordinate of guide star (pixels)
    imgsize : int
        Dimension of subarray (pixels)
    show_fig : bool, optional
        Denotes whether to generate a plot of the created subarray

    Returns
    -------
    img : 2-D numpy array
        Subarray of image data around guide star coordinates
    """
    xlow = int(xcoord - imgsize / 2)
    ylow = int(ycoord - imgsize / 2)

    xhigh = int(xcoord + imgsize / 2)
    yhigh = int(ycoord + imgsize / 2)

    img = image[ylow:yhigh, xlow:xhigh]

    if show_fig:
        plt.figure()
        plt.imshow(img, cmap='Greys_r')
        plt.show()

    return img


def shift_to_id_attitude(image, root, guider, out_dir, guiding_selections_file,
                         all_found_psfs_file, center_pointing_file, psf_center_file=None,
                         logger_passed=False):
    """Shift the FGS image such that the guide star is at the ID
    attitude. Rewrite the FGS FITS file, guiding_selections, and all_found_psfs
    catalog files for this shifted case. (Rename old versions of those files
    to be "unshifted_".)

    Parameters
    ----------
    image : 2D numpy array
        Image data
    root : str
        Name used to create the output directory, {out_dir}/out/{root}
    guider : int
        Guider number (1 or 2)
    out_dir : str
        Where output files will be saved.
    guiding_selections_file : str
        Path to existing unshifted_guiding_selections_{root}_G{guider}.txt
    all_found_psfs_file : str
        Path to existing unshifted_gall_found_psfs_{root}_G{guider}.txt
    center_pointing_file : str
        Path to existing center_pointing_{root}_G{guider}.txt
    psf_center_file : str, optional
        Path to existing unshifted_psf_center_{root}_G{guider}.txt
        center psf file
    logger_passed : bool
        T/F if a logger is passed into the function

    Returns
    -------
    shifted_image : 2D numpy array
        The image data shifted so that the guide star falls on the
        ID attitude
    shifted_guiding_selections : str
        Path to the shifted guiding selections file which should be
        used in future calculations over the unshift data
    psf_center_file : str or None
        Path to shifted psf center file if it exists
    """
    # Set up logger
    if not logger_passed:
        out_dir_root = out_dir.split(root)[0] + root
        utils.create_logger_from_yaml(__name__, out_dir_root=out_dir_root, root=root, level='DEBUG')

    if isinstance(image, str):
        image = fits.getdata(image)

    # 0) Fetch existing information and determine ID attitude
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Define input file info
    if 'guiding_selections_' in guiding_selections_file:
        file_root = guiding_selections_file.split('guiding_selections_')[-1].split('.txt')[0]
    elif 'regfile' in guiding_selections_file:
        file_root = guiding_selections_file.split('/')[-1].split('_regfile.txt')[0]
    else:
        raise ValueError('Guiding selections file {} has a naming structure that cannot be '
                         'parsed.'.format(guiding_selections_file))

    # Add guider number to file root if not present - common for old regfiles
    if 'G{}'.format(guider) not in file_root:
        file_root += '_G{}'.format(guider)

    # Load the catalogs with the unshifted data
    guiding_selections_cat = asc.read(guiding_selections_file)
    all_found_psfs_cat = asc.read(all_found_psfs_file)

    # Set the pixel shift
    xend, yend = (1023.5, 1023.5)  # ID attitude
    hdr_keyword = '{}'.format((xend, yend))

    # 1) Shift the image array
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    xstart, ystart = guiding_selections_cat['x', 'y'][0] # Guide star location
    dx = xend - xstart
    dy = yend - ystart

    if (dx, dy) == (0, 0):
        LOGGER.info('FSW File Writing: No need to shift file; guide star already at ID attitude.')
        return image, guiding_selections_file, psf_center_file

    bkg = np.median(image)

    LOGGER.info("FSW File Writing: Shifting guide star to ID attitude ({}, {})".format(xend, yend))
    shifted_image = shift(image, (dy, dx), mode='constant', cval=bkg, prefilter=True)
    # TODO: ? Should we be shifting the seed image and re-adding the bias on top of it?
    # TODO: Save out seed and bias separately, shift seed here, and re-add bias on top of it

    # 2) Write new shifted guiding_selections*.txt
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Shift the guiding selections catalog
    shifted_guiding_selections_cat = guiding_selections_cat.copy()
    shifted_guiding_selections_cat['x'] += dx
    shifted_guiding_selections_cat['y'] += dy

    shifted_guiding_selections = os.path.join(out_dir, 'shifted_guiding_selections_{}.txt'.format(file_root))

    # Write new guiding_selections*.txts
    utils.write_cols_to_file(shifted_guiding_selections,
                             labels=['y', 'x', 'countrate'],
                             cols=shifted_guiding_selections_cat,
                             log=LOGGER)

    # 3) Write new shifted all_found_psfs*.txt
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Shift the all_found_psfs*.txt catalog
    shifted_all_psfs_cat = all_found_psfs_cat.copy()
    shifted_all_psfs_cat['x'] += dx
    shifted_all_psfs_cat['y'] += dy

    shifted_all_psfs = os.path.join(out_dir, 'shifted_all_found_psfs_{}.txt'.format(file_root))

    all_cols = utils.create_cols_for_coords_counts(shifted_all_psfs_cat['x'],
                                                   shifted_all_psfs_cat['y'],
                                                   shifted_all_psfs_cat['countrate'],
                                                   None,
                                                   labels=shifted_all_psfs_cat['label'],
                                                   inds=range(len(shifted_all_psfs_cat['x'])))

    # Write new all_found_psfs*.txts
    utils.write_cols_to_file(shifted_all_psfs,
                             labels=['label', 'y', 'x', 'countrate'],
                             cols=all_cols,
                             log=LOGGER)

    # 4) Write new shifted center_pointing*.txt
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    center_pointing_cat = asc.read(center_pointing_file, format='commented_header', delimiter=',')

    shifted_center_pointing_cat = center_pointing_cat.copy()
    label = center_pointing_cat.colnames[0]
    if isinstance(shifted_center_pointing_cat[label][0], str):
        y, x = [float(i) for i in shifted_center_pointing_cat[label][0].split(' ')]
        shifted_center_pointing_cat = [y+dy, x+dx]
    else:
        shifted_center_pointing_cat = shifted_center_pointing_cat[label][0]

    shifted_center_pointing = os.path.join(out_dir, 'shifted_center_pointing_{}.txt'.format(file_root))

    # Write new center_pointing*.txt
    utils.write_cols_to_file(shifted_center_pointing,
                             labels=[label],
                             cols=[shifted_center_pointing_cat],
                             log=LOGGER)

    # 5) Write new shifted psf_center*.txt
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if psf_center_file is not None:
        psf_center_cat = asc.read(psf_center_file)

        shifted_psf_center_cat = psf_center_cat.copy()
        shifted_psf_center_cat['x'] += dx
        shifted_psf_center_cat['y'] += dy

        shifted_psf_center = os.path.join(out_dir, 'shifted_psf_center_{}.txt'.format(file_root))
        psf_center_file = shifted_psf_center

        # Write new psf_center*.txts
        utils.write_cols_to_file(shifted_psf_center,
                                 labels=['y', 'x', 'countrate'],
                                 cols=shifted_psf_center_cat,
                                 log=LOGGER)

    # 6) Rewrite the shifted FGS image and save old file as _unshifted.fits
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Load header file (try to use the converted image to get distortion information)
    header_file = os.path.join(out_dir.split(root)[0], root, 'FGS_imgs', f'unshifted_{root}_G{guider}.fits')
    if not os.path.exists(header_file):
        header_file = os.path.join(DATA_PATH, 'header_g{}.fits'.format(guider))
    hdr = fits.getheader(header_file, ext=0)
    hdr['IDATTPIX'] = (hdr_keyword, 'Image shifted to place GS at ID attitude')

    shifted_FGS_img = os.path.join(out_dir, 'FGS_imgs', 'shifted_' + file_root + '.fits')

    # Write new FITS files
    # Get rid of any artifacts from shifting
    saved_shifted_image = np.copy(shifted_image)
    saved_shifted_image[saved_shifted_image < 1e-8] = 0.
    utils.write_fits(shifted_FGS_img, [None, saved_shifted_image], header=[hdr, None], log=LOGGER)

    return shifted_image, shifted_guiding_selections, psf_center_file


def bright_guiding_check(countrate_list, threshold_factor, normal_ops=False, override_bright_guiding=False):
    """
    Check the 3x3 count rate to see if it is above the OSS trigger, if so, adjust
    the threshold accordingly. Returns list of thresholds and a new threshold factor
    based off the guide star, if required

    Parameters
    ----------
    countrate_3x3: array
        Array of all 3x3 count rates of each selected PSF
    threshold_factor: float between 0 and 1
        This is the user defined threshold factor
    normal_ops: boolean
        Whether this is normal operations. This is used during OTE commissioning
        when Use OSS Defaults is selected for making a POF
    override_bright_guiding: boolean
        If the user wants to guarentee that their provided threshold factor will be used,
        regardless of the 3x3 count rate, they will set this parameter to True. When set
        to False, if the 3x3 count rate is above the OSS trigger, the threshold and threshold
        factors will be replaced.
    """
    # If we want to match the normal operations behavior of OSS
    if normal_ops and override_bright_guiding:
        raise ValueError('Cannot use normal ops functionality and and overwrite the bright guiding functionality at the same time. Please use only one, or neither')

    if normal_ops:
        dim_star_threshold_factor = DIM_STAR_THRESHOLD_FACTOR
    else:
        dim_star_threshold_factor = threshold_factor

    # Check if the 3x3 count rate of the guide star is above the OSS trigger value
    if (countrate_list[0] < OSS_TRIGGER):
        threshold = countrate_list * dim_star_threshold_factor
    else:
        bright_threshold = countrate_list - BRIGHT_STAR_THRESHOLD_ADDEND
        user_threshold = countrate_list * threshold_factor

        # If the user wants to use their threshold no matter what
        if override_bright_guiding:
            threshold = user_threshold
            msg = f" but the user has forced a threshold factor of {threshold_factor}"
        # Check that the bright threshold is larger than user supplied threshold
        elif user_threshold[0] > bright_threshold[0] and not normal_ops:
            threshold = user_threshold
            msg = f". The user supplied a threshold that is larger than the suggestion, so the user-supplied threshold factor of {threshold_factor} will be applied."
        else:
            msg = ", which will be written out to the file"
            threshold = bright_threshold

        LOGGER.warning(f"The selected guide star triggers bright guiding. A successful threshold factor was calculated to be greater than or equal to {np.round(bright_threshold[0]/countrate_list[0], 3)}{msg}")

    # Based the count rate factor off of the guide star
    return threshold, threshold[0]/countrate_list[0]
