'''Writes all flight software files for ID, ACQ, and/or TRK steps

This module creates an FGS simulation object for ID, ACQ, and/or TRK
steps, which it uses to create the necessary flight software files for
use with the DHAS, the FGSES/CertLab, or just for user inspection. The
files created for each step are as follows:

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
        from jwst_fgs_commissioning_tools.fsw_file_writer import buildfgssteps
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
        ``configfile`` - file definiing parameters for each guider step
        ``out_dir`` - where output FGS image(s) will be saved. If not
            provided, the image(s) will be saved to ../out/{root}.
'''

# STDLIB
import os

# Third Party
from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

# LOCAL
from jwst_fgs_commissioning_tools import utils, log
from jwst_fgs_commissioning_tools.fsw_file_writer import config, getbias, write_files

# DEFINE ALL PATHS
FSW_PATH = os.path.dirname(os.path.realpath(__file__))
PACKAGE_PATH = os.path.split(FSW_PATH)[0]
OUT_PATH = os.path.split(PACKAGE_PATH)[0]  # Location of out/ and logs/ directory
DATA_PATH = os.path.join(PACKAGE_PATH, 'data')

class BuildFGSSteps(object):
    '''
    Creates an FGS simulation object for ID, ACQ, and/or TRK stages to be used
    with DHAS.
    '''
    def __init__(self, im, guider, root, step, reg_file=None, configfile=None,
                 out_dir=None):
        # Practical things
        self.guider = guider
        self.root = root
        self.step = step
        self.yoffset = 12

        ## DEFINE ALL THINGS PATHS
        self.out_dir = utils.make_out_dir(out_dir, OUT_PATH, root)
        utils.ensure_dir_exists(os.path.join(self.out_dir, 'dhas'))
        utils.ensure_dir_exists(os.path.join(self.out_dir, 'ground_system'))
        utils.ensure_dir_exists(os.path.join(self.out_dir, 'stsci'))

        ## READ IN IMAGE
        if isinstance(im, str):
            data = fits.getdata(im) #*_bin_norm from FGS_bin_tool
            self.input_im = data
        else:
            self.input_im = im

        # Correct for negative, saturated pixels and other nonsense
        self.input_im = utils.correct_image(self.input_im)

        # THEN convert to uint16
        self.input_im = np.uint16(self.input_im)

        log.info('FSW File Writing: Max of input image: {}'.format(np.max(self.input_im)))

        self.get_coords_and_counts(reg_file=reg_file)

        section = '{}_dict'.format(self.step.lower())
        config_ini = self.build_step(section, configfile)
        self.image = self.create_img_arrays(section, config_ini)
        self.write()


    ### Guide star and reference star coordinates and countrates
    def get_coords_and_counts(self, reg_file=None):
        '''
        Get coordinate information of guide star and reference stars
        '''
        if reg_file is None:
            reg_file = os.path.join(self.out_dir,
                                    '{0}_G{1}_regfile.txt'.format(self.root,
                                                                  self.guider))
        log.info("FSW File Writing: Using {} as the reg file".format(reg_file))

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
        '''
        Build the step of commissioning based on the parameters in the config.ini
        and return the step section if step is known
        '''
        self.nreads = 2

        if configfile is None:
            configfile = os.path.join(DATA_PATH, 'config.ini')
        config_ini = config.load_config_ini(configfile)

        if self.step != 'ID':
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
        '''
        Create a noisy sky image for ID and ACQ steps.
        There is only an added poissonfactor for ACQ2.
        PoissonNoise should always be set to True for ID.
        '''
        #Create the time-normalized image
        self.time_normed_im = self.input_im * config_ini.getfloat(section, 'tcds')

        ## Grab the expected bias
        if config_ini.getboolean(section, 'bias'):
            self.bias = getbias.getbias(self.guider, self.xarr, self.yarr,
                                        self.nreads, config_ini.getint(section, 'nramps'),
                                        config_ini.getint(section, 'imgsize'))

            ## Take the bias and add a noisy version of the input image, adding signal
            ## over each read
            image = np.copy(self.bias)
            for ireads in range(self.nreads):
                image[ireads::(self.nreads)] += (ireads + 1) * \
                                                np.random.poisson(self.time_normed_im)
        else:
            self.bias = None
            image = self.time_normed_im

        ## Cut any pixels over saturation or under zero
        image = utils.correct_image(image)

        if config_ini.getboolean(section, 'cdsimg'):
            self.cds = create_cds(image)
        else:
            self.cds = None

        if config_ini.getboolean(section, 'stripsimg'):
            self.strips = create_strips(image,
                                        config_ini.getint(section, 'imgsize'),
                                        config_ini.getint(section, 'nstrips'),
                                        config_ini.getint(section, 'nramps'),
                                        self.nreads,
                                        config_ini.getint(section, 'height'),
                                        self.yoffset,
                                        config_ini.getint(section, 'overlap'))
        if self.step == 'LOSTRK':
            # Normalize to a count sum of 1000
            image = image / np.sum(image) * 1000
            # Resize image array to oversample by 6 (from 43x43 to 255x255)
            image = image.repeat(6, axis=0)
            image = image.repeat(6, axis=1)
            image = image[1:-2, 1:-2]

        return image

    def write(self):
        write_files.write_all(self)

#-------------------------------------------------------------------------------
def create_strips(image, imgsize, nstrips, nramps, nreads, h, yoffset, overlap):
    '''
    Create the ID strips fits image to be passed into DHAS
    '''
    nz = nramps * nreads
    strips = np.zeros((nstrips * nz, h, imgsize))
    nn = 0
    for i in range(nstrips):
        for iz in range(nz):
            ylow = i * (h - overlap) + yoffset
            yhigh = ylow + h
            strips[nn] = image[iz, ylow:yhigh, :]
            nn += 1

    # Make sure the data is between 0 and 65,000 counts and are finite numbers
    strips = utils.correct_image(strips)
    strips = np.uint16(strips)

    return strips


def create_cds(arr):
    '''
    Create CDS image: Subtract the first read from the second read.
    '''
    return arr[1::2] - arr[:-1:2]

def display(image, ind=0, vmin=None, vmax=None, xarr=None, yarr=None):
    '''
    Display an image array. If the array is a cube, specify the index to
    look at using 'ind'. If you want to add a scatter plot of PSF centers,
    use add_coords=True and pass in x and y (Note: This only works for full
    frame images for the time being).
    '''
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

### Arrays and array manipulation
def add_jitter(cube, total_shift=3):
    '''
    Add random single pixel jitter

    VERY rudimentary. Uses np.roll so images look kind of funny.
    '''
    cube2 = np.zeros_like(cube)
    for i, img in enumerate(cube):
        # This will generate a random integer up to the total_shift
        shift = np.random.randint(total_shift + 1)
        cube2[i] = np.roll(img, shift)

    return cube2

def create_im_subarray(image, xcoord, ycoord, imgsize, show_fig=False):
    '''
    Based on the array size given by nx and ny, created a subarray around
    the guide star.
    '''
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


def add_background(array, imgsize, nz):
    '''
    Add background to array
    '''
    try:
        array += 500 + 10. * np.random.standard_normal((nz, imgsize, imgsize))
    except NameError:
        array = 500 + 10. * np.random.standard_normal((nz, imgsize, imgsize))

    return array
