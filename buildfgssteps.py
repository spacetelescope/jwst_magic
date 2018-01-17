# STDLIB
import os

# Third Party
from astropy.io import fits
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import numpy as np
import pickle

# LOCAL
from getbias import getbias
import utils
import log
import select_psfs
import write_files

# DEFINE ALL NECESSARY CONSTANTS
LOCAL_PATH = os.path.dirname(os.path.realpath(__file__))
DATA_PATH = os.path.join(LOCAL_PATH, 'data')

class BuildFGSSteps(object):
    '''
    Creates an FGS simulation object for ID, ACQ, and/or TRK stages to be used
    with DHAS.
    '''
    def __init__(self, im, guider, root, step, reg_file=None, nramps=None):
        # Practical things
        self.guider = guider
        self.root = root
        self.step = step
        self.yoffset = 12

        ## DEFINE ALL THINGS PATHS
        self.out_dir = os.path.join(LOCAL_PATH, 'out', root)
        utils.ensure_dir_exists(os.path.join(self.out_dir, 'dhas'))
        utils.ensure_dir_exists(os.path.join(self.out_dir, 'ground_system'))
        utils.ensure_dir_exists(os.path.join(self.out_dir, 'stsci'))

        ## READ IN IMAGE
        if isinstance(im, str):
            data = fits.getdata(im) #*_bin_norm from FGS_bin_tool
            self.input_im = np.uint16(data)
        else:
            self.input_im = np.uint16(im)

        # Correct for negative, saturated pixels and other nonsense
        self.input_im = utils.correct_image(self.input_im)
        log.info('Max of input image: {}'.format(np.max(self.input_im)))

        self.get_coords_and_counts(reg_file=reg_file)

        step_dict = self.build_step(nramps)
        self.image = self.create_img_arrays(step_dict)
        self.write()


    ### Guide star and reference star coordinates and countrates
    def get_coords_and_counts(self, reg_file=None):
        '''
        Get coordinate information of guide star and reference stars
        '''
        if reg_file is None:
            reg_file = os.path.join(os.path.join(LOCAL_PATH, 'out', self.root,
                                                 '{0}_G{1}_regfile.txt'.format(self.root,
                                                                               self.guider)))
        log.info("Using {} as the reg file".format(reg_file))
        if reg_file.endswith('reg'):
            self.xarr, self.yarr = np.loadtxt(reg_file)
            self.countrate = []
            for xa, ya, in zip(self.xarr, self.yarr):
                self.countrate.append(select_psfs.countrate_3x3(xa, ya, self.input_im))
        else:
            self.yarr, self.xarr, self.countrate = np.loadtxt(reg_file, delimiter=' ',
                                                              skiprows=1).T
        # Add y offset to all coordinates
        self.yarr = self.yarr - self.yoffset

        # Cover cases where there is only one entry in the reg file
        try:
            len(self.xarr)
        except TypeError:
            self.xarr = np.asarray([self.xarr])
            self.yarr = np.asarray([self.yarr])
            self.countrate = np.asarray([self.countrate])

    def build_step(self, configfile=None):
        '''
        Build the step of commissioning and return the step dictionary if step
        is known
        '''
        self.nreads = 2

        if not configfile:
            parameters = pickle.load(open(os.path.join(DATA_PATH,
                                                       'master_config.p'), "rb"))
        else:
            parameters = pickle.load(open(configfile), "rb")

        if self.step == 'ID':
            step_dict = parameters['id_dict']
            self.acqnum = step_dict['acqNum']

        else:
            step_dict = parameters['{}_dict'.format(self.step.lower())]

            self.xarr = np.asarray([self.xarr[0]])
            self.yarr = np.asarray([self.yarr[0]])
            self.countrate = np.asarray([self.countrate[0]])
            self.input_im = create_im_subarray(self.input_im, self.xarr,
                                               self.yarr, step_dict['imgsize'])
            self.acqnum = step_dict['acqNum']

        return step_dict


    def create_img_arrays(self, step_dict):
        '''
        Create a noisy sky image for ID and ACQ steps.
        There is only an added poissonfactor for ACQ2.
        PoissonNoise should always be set to True for ID.
        '''
        #Create the time-normalized image
        self.time_normed_im = self.input_im * step_dict['tcds']

        ## Grab the expected bias
        self.bias = getbias(self.guider, self.xarr, self.yarr, self.nreads,
                            step_dict['nramps'], step_dict['imgsize'])

        ## Take the bias and add a noisy version of the input image, adding signal
        ## over each read
        image = np.copy(self.bias)
        for ireads in range(self.nreads):
            image[ireads::(self.nreads)] += (ireads + 1) * \
                                         np.random.poisson(self.time_normed_im)
        ## Cut any pixels over saturation or under zero
        image = utils.correct_image(image)

        if step_dict['cdsimg']:
            self.cds = create_cds(image)
        else:
            self.cds = False
        if step_dict['stripsimg']:
            self.strips = create_strips(image, step_dict['imgsize'],
                                        step_dict['nstrips'],
                                        step_dict['nramps'],
                                        self.nreads,
                                        step_dict['height'],
                                        self.yoffset,
                                        step_dict['overlap'])
        return image

    def write(self):
        write_files.write_all(self)


#-------------------------------------------------------------------------------
def create_strips(image, imgsize, nstrips, nramps, nreads, h, yoffset, overlap):
    '''
    Create the ID strips fits image to be passed into DHAS
    '''
    nz = nramps * nreads
    strips = np.zeros((nstrips*nz, h, imgsize))
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
    return arr[1::2]-arr[:-1:2]

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
        shift = np.random.randint(total_shift+1)
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
