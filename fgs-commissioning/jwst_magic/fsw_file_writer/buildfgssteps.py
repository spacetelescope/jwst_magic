"""Setup to write all simulated files for CAL, ID, ACQ, and/or TRK steps

This module parses the ``config.ini`` file to create an FGS simulation
object for CAL, ID, ACQ, and/or TRK steps. This object will ultimately
pass to ``write_files.py`` to create the necessary flight software files
for use with the DHAS, the FGSES/CertLab, or just for user inspection.

Authors
-------
    - Keira Brooks
    - Lauren Chambers

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
        ``reg_file`` - file containing X/Y positions and countrates for
            all stars in an image
        ``configfile`` - file defining parameters for each guider step.
            If not defined, defaults to jwst_magic/data/config.ini
        ``out_dir`` - where output files will be saved. If not provided,
            the image(s) will be saved within the repository at
            tools/fgs-commissioning/
        ``logger_passed`` - denotes if a logger object has already been
            generated.
"""

# Standard Library Imports
import os
import logging

# Third Party Imports
from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

# Local imports
from .. import utils
from ..fsw_file_writer import config, detector_effects, write_files

# Paths
__location__ = os.path.realpath(os.path.join(os.getcwd(), os.path.dirname(__file__)))
PACKAGE_PATH = os.path.split(__location__)[0]
OUT_PATH = os.path.split(PACKAGE_PATH)[0]  # Location of out/ and logs/ directory
DATA_PATH = os.path.join(PACKAGE_PATH, 'data')

# Start logger
LOGGER = logging.getLogger(__name__)

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# MAIN CLASS
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


class BuildFGSSteps(object):
    """Creates an FGS simulation object for ID, ACQ, and/or TRK stages
    to be used with DHAS.
    """
    def __init__(self, im, guider, root, step, reg_file=None, configfile=None,
                 out_dir=None, logger_passed=False):
        """Initialize the class and call build_fgs_steps().
        """
        # Set up logger
        if not logger_passed:
            utils.create_logger_from_yaml(__name__, root=root, level='DEBUG')

        try:
            # Initialize attributes
            self.guider = guider
            self.root = root
            self.step = step
            self.yoffset = 12

            # Build FGS steps
            self.build_fgs_steps(im, root, out_dir, reg_file, configfile)

        except Exception as e:
            LOGGER.exception(e)
            raise

    def build_fgs_steps(self, im, root, out_dir, reg_file, configfile):
        """Creates an FGS simulation object for ID, ACQ, and/or TRK stages
        to be used with DHAS.

        Parameters
        ----------
        im : 2-D numpy array or str
            Image data or filepath to image
        root : str
            Name used to create the output directory, {out_dir}/out/{root}
        out_dir : str
            Where output files will be saved. If not provided, the
            image(s) will be saved within the repository at
            tools/fgs-commissioning/
        reg_file : str
            File containing X/Y positions and countrates for all stars
            in the provided image
        configfile : str
            File defining parameters for each guider step. If not
            defined, defaults to jwst_magic/data/config.ini
        """
        # Define paths
        self.out_dir = utils.make_out_dir(out_dir, OUT_PATH, root)
        utils.ensure_dir_exists(os.path.join(self.out_dir, 'dhas'))
        utils.ensure_dir_exists(os.path.join(self.out_dir, 'ground_system'))
        utils.ensure_dir_exists(os.path.join(self.out_dir, 'stsci'))

        # READ IN IMAGE
        if isinstance(im, str):
            data = fits.getdata(im)
            self.input_im = data
        else:
            self.input_im = im

        # *** SHOULD THIS BE HAPPENING HERE?? ***
        # Correct for negative, saturated pixels and other nonsense
        self.input_im = utils.correct_image(self.input_im)

        # # THEN convert to uint16
        # self.input_im = np.uint16(self.input_im)

        # LOGGER.info('FSW File Writing: Max of input image: {}'.format(np.max(self.input_im)))

        self.get_coords_and_counts(reg_file=reg_file)

        section = '{}_dict'.format(self.step.lower())
        config_ini = self.build_step(section, configfile)
        self.image = self.create_img_arrays(section, config_ini)

        # Write the files
        LOGGER.info("FSW File Writing: Creating {} FSW files".format(self.step))
        write_files.write_all(self)

    def get_coords_and_counts(self, reg_file=None):
        """Get coordinate information and countrates of guide star and
        reference stars.

        Parameters
        ----------
        reg_file : str, optional
            File containing X/Y positions and countrates for all stars
            in the provided image
        """
        if reg_file is None:
            reg_file = os.path.join(self.out_dir,
                                    '{0}_G{1}_regfile.txt'.format(self.root,
                                                                  self.guider))
        LOGGER.info("FSW File Writing: Using {} as the reg file".format(reg_file))

        if reg_file.endswith('reg'):
            self.xarr, self.yarr = np.loadtxt(reg_file)
            self.countrate = []
            for xa, ya, in zip(self.xarr, self.yarr):
                self.countrate.append(utils.countrate_3x3(xa, ya, self.input_im))
        else:
            self.yarr, self.xarr, self.countrate = np.loadtxt(reg_file, delimiter=' ',
                                                              skiprows=1).T
        # Add y offset to all coordinates
        # self.yarr = self.yarr - self.yoffset

        # Cover cases where there is only one entry in the reg file
        try:
            len(self.xarr)
        except TypeError:
            self.xarr = np.asarray([self.xarr])
            self.yarr = np.asarray([self.yarr])
            self.countrate = np.asarray([self.countrate])

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

    def create_img_arrays(self, section, config_ini):
        """Create the needed image arrays for the given step.

        Possible arrays (created as class attributes):
            - time_normed_im : the "sky" image, or the time-normalized
                image (in counts per one second, AKA counts)
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
        section : str
            Name of step within the config file ('{step}_dict')
        configfile : str, optional
            File defining parameters for each guider step. If not
            defined, defaults to jwst_magic/data/config.ini

        Returns
        -------
        image : 2-D numpy array
            The final image including any necessary detector effects,
            either full-frame or the appropriately sized subarray, in counts.
        """

        # Create the time-normalized image (will be in counts, where the
        # input_im is in counts per second)
        # ** THIS MIGHT NOT BE CORRECT IF IT OCCURS AFTER THE IMAGE HAS BEEN NORMALIZED TO COUNTS ALREADY **
        self.time_normed_im = self.input_im * config_ini.getfloat(section, 'tframe')

        # Add the bias, and build the array of reads with noisy data
        if config_ini.getboolean(section, 'bias'):
            # Get the bias ramp
            nramps = config_ini.getint(section, 'nramps')
            det_eff = detector_effects.FGSDetectorEffects(
                self.guider, self.xarr, self.yarr, self.nreads, nramps,
                config_ini.getint(section, 'imgsize')
            )
            self.bias = det_eff.add_detector_effects()

            # Take the bias and add a noisy version of the input image
            # (specifically Poisson noise), adding signal over each read
            image = np.copy(self.bias)
            signal_cube = [self.time_normed_im] * nramps

            # First read
            # Calculate how much signal should be in the first read
            i_read = 0
            n_drops_before_first_read = int(config_ini.getfloat(section, 'ndrops1'))
            n_frametimes_in_first_read = n_drops_before_first_read + 1
            # Add signal to every first read
            noisy_signal_cube = np.random.poisson(signal_cube)
            image[i_read::(self.nreads)] += n_frametimes_in_first_read * noisy_signal_cube

            # Second read
            # Calculate how much signal should be in the second read
            i_read = 1
            n_drops_before_second_read = int(config_ini.getfloat(section, 'ndrops2'))
            n_frametimes_in_second_read = n_frametimes_in_first_read + n_drops_before_second_read + 1
            # Add signal to every second read
            noisy_signal_cube = np.random.poisson(signal_cube)
            image[i_read::(self.nreads)] += n_frametimes_in_second_read * noisy_signal_cube

        else:
            # In the case of LOSTRK, just return one frame with one
            # frame readout worth of signal
            self.bias = None
            image = self.time_normed_im

        # Cut any pixels over saturation or under zero
        # image = utils.correct_image(image)

        # Create the CDS image by subtracting the first read from the second
        # read, for each ramp
        if config_ini.getboolean(section, 'cdsimg'):
            self.cds = create_cds(image, section, config_ini)
        else:
            self.cds = None

        # If in ID, split the full-frame image into strips
        if config_ini.getboolean(section, 'stripsimg'):
            self.strips = create_strips(image,
                                        config_ini.getint(section, 'imgsize'),
                                        config_ini.getint(section, 'nstrips'),
                                        config_ini.getint(section, 'nramps'),
                                        self.nreads,
                                        config_ini.getint(section, 'height'),
                                        self.yoffset,
                                        config_ini.getint(section, 'overlap'))

        # Modify further for LOSTRK images (that will be run in FGSES)
        if self.step == 'LOSTRK':
            # Normalize to a count sum of 1000
            image = image / np.sum(image) * 1000
            # Resize image array to oversample by 6 (from 43x43 to 255x255)
            image = image.repeat(6, axis=0)
            image = image.repeat(6, axis=1)
            image = image[1:-2, 1:-2]

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
        Numer of ramps in the integraion
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

    # Make sure the data is between 0 and 65,000 counts and are finite numbers
    # strips = utils.correct_image(strips)

    return strips


def create_cds(arr, section, config_ini, fix_saturated_pix=True):
    """Create CDS image: Subtract the first read from the second read.
    Option to handle saturated pixels in CDS.

    Parameters
    ----------
    arr : 3-D numpy array
        Image data in
    fix_saturated_pix : boolean
        Apply a fix to saturated pixels in the CDS array?

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
        saturated_read_2 = second_reads == 65000
        saturated_read_1 = first_reads == 65000

        # For pixels that are saturated in the second read, calculate their
        # expected CDS value using the count rate from the first read and
        # assuming linearity.
        n_drops_before_first_read = int(config_ini.getfloat(section, 'ndrops1'))
        n_frametimes_in_first_read = n_drops_before_first_read + 1
        time_first_read = config_ini.getfloat(section, 'tframe') * n_frametimes_in_first_read  # seconds
        first_read_countrates = first_reads / time_first_read  # counts / second

        n_drops_before_second_read = int(config_ini.getfloat(section, 'ndrops2'))
        n_frametimes_in_cds_read = n_drops_before_second_read + 1
        time_cds_read = config_ini.getfloat(section, 'tframe') * n_frametimes_in_cds_read  # seconds

        # If the calculated CDS value is less than saturation, set that
        # as the pixel value. Otherwise, set the pixel value to 65000.
        cds_counts = first_read_countrates * time_cds_read  # counts
        cds_counts[cds_counts > 65000] = 65000
        cds_arr[saturated_read_2] = cds_counts[saturated_read_2]

        # n_sat_2 = len([p for p in saturated_read_2[0].flatten() if p])
        n_sat_2 = len(saturated_read_2[0][saturated_read_2[0] == True].flatten())
        print('Adjusting {} pixels that are saturated in read 2.'.format(n_sat_2))

        # For pixels that are saturated in both reads, set their CDS
        # value to the saturated value (65000).
        saturated_both_reads = saturated_read_1 * saturated_read_2
        cds_arr[saturated_both_reads] = 65000
        print(saturated_both_reads)

        # n_sat_both = len([p for p in saturated_both_reads[0].flatten() if p])
        n_sat_both = len(saturated_both_reads[0][saturated_both_reads[0] == True].flatten())
        print('Adjusting {} pixels that are saturated in both reads.'.format(n_sat_both))

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
    if imgsize % 2 == 1:
        xlow = int(xcoord - (imgsize // 2 + 1))
        ylow = int(ycoord - (imgsize // 2 + 1))
    else:
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
