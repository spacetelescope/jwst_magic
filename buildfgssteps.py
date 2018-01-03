# STDLIB
import os

# Third Party
from astropy.io import fits
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import numpy as np

# LOCAL
from get_bias import getbias
import utils
import log
import select_psfs
import write_files

# DEFINE ALL NECESSARY CONSTANTS
TCDSID = 0.338    # neil 18 may 12
TCDSACQ1 = 0.3612 # 20 apr 17 keira # 1cds time = 2*(128*128*10.e-6) = 0.32768s
TCDSACQ2 = 0.0516 # 20 apr 17 keira brooks # 2cds time = 4*(32*32*10e-6) = 0.04096s
TCDSTRK = 0.0256  # neil 18 may 12
TCDSFG = 0.0512   # neil 18 may 12

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

        step_dict = self.build_step(step, nramps)
        self.create_img_arrays(step_dict)
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
        # Cover cases where there is only one entry in the reg file
        self.xarr = np.asarray(self.xarr)
        self.yarr = np.asarray(self.yarr)
        self.countrate = np.asarray(self.countrate)

    def build_step(self, step, nramps):
        '''
        Build the step of commissioning and return the step dictionary if step
        is known
        '''
        self.nreads = 2

        if step == 'ID':
            return self.build_id()
        elif step == 'ACQ1':
            return self.build_acq1()
        elif step == 'ACQ2':
            return self.build_acq2()
        elif step == 'TRK':
            return self.build_trk(nramps=nramps)
        elif step == 'LOSTRK':
            return self.build_lostrk()
        else:
            log.error('No known step requested.')
            raise StandardError('No known step requested. Exiting')

    def build_id(self, yoffset=0, h=64, overlap=8):
        '''
        Build the ID step:
        nx, ny = 2048 (or size of truth image)
        nramps = 2
        nstrips = [32,33,34,35,36,37,39,40,42] for overlap = [0,2,4,6,8,10,12,14,16,18]
        '''
        step_dict = build_id_dict(yoffset=0, h=64, overlap=8)

        return step_dict


    def build_acq1(self):
        '''
        Build the ACQ1 step:
        nx, ny = 128
        nramps = 6
        ACQ1 = repeat 6x(reset, drop, read, drop, read)
        '''
        step_dict = self.build_acq1_dict()

        self.input_im = create_im_subarray(self.input_im, self.xarr,
                                           self.yarr, step_dict['imgsize'])
        return step_dict


    def build_acq2(self):
        '''
        Build the ACQ2 step:
        nx, ny = 32
        nramps = 5
        ACQ2 = repeat 5x(reset, drop, read, drop, drop, drop, read, drop)
        '''
        step_dict = self.build_acq2_dict()

        self.input_im = create_im_subarray(self.input_im, self.xarr,
                                           self.yarr, step_dict['imgsize'])
        return step_dict


    def build_trk(self, nramps):
        '''
        Build the TRK step:
        nx, ny = 32
        nramps = a large number (i.e. 5000+)
        '''
        step_dict = self.build_trk_dict(nramps)

        self.input_im = create_im_subarray(self.input_im, self.xarr,
                                           self.yarr, step_dict['imgsize'])
        return step_dict


    def build_lostrk(self):
        '''
        Build the Line of Sight - TRK step:
        nx, ny = 43
        nramps = a large number (i.e. 5000+)
        '''
        step_dict = self.build_lostrk_dict()

        self.input_im = create_im_subarray(self.input_im, self.xarr,
                                           self.yarr, step_dict['imgsize'])
        return step_dict

    def create_img_arrays(self, step_dict):
        '''
        Create a noisy sky image for ID and ACQ steps.
        There is only an added poissonfactor for ACQ2.
        PoissonNoise should always be set to True for ID.
        '''
        #Create the time-normalized image
        time_normed_im = self.input_im * step_dict['tcds']

        ## Grab the expected bias
        bias = getbias(self.guider, self.xarr, self.yarr, self.nreads,
                       step_dict['nramps'], step_dict['imgsize'])

        ## Take the bias and add a noisy version of the input image, adding signal
        ## over each read
        image = np.copy(bias)
        for ireads in range(self.nreads):
            image[ireads::(self.nreads)] += (ireads + 1) * \
                                         np.random.poisson(time_normed_im)
        ## Cut any pixels over saturation or under zero
        image = utils.correct_image(image)

        if step_dict['cdsimg']:
            self.cds = create_cds(image)
        if step_dict['stripsimg']:
            self.strips = create_strips(image, step_dict['imgsize'],
                                        step_dict['nstrips'],
                                        step_dict['nramps'],
                                        step_dict['nreads'],
                                        step_dict['height'],
                                        step_dict['yoffset'],
                                        step_dict['overlap'])
        return image

    def write(self):
        write_files.write_all(self)


    @staticmethod
    def build_id_dict():
        '''
        Build the ID step:
        nx, ny = 2048 (or size of truth image)
        nramps = 2
        nstrips = [32,33,34,35,36,37,39,40,42] for overlap = [0,2,4,6,8,10,12,14,16,18]
        '''
        id_dict = {}
        id_dict['imgsize'] = 2048
        id_dict['nramps'] = 2

        id_dict['overlap'] = 8
        id_dict['height'] = 64
        id_dict['yoffset'] = 0

        id_dict['tcds'] = TCDSID
        id_dict['step'] = 'ID'

        id_dict['cdsimg'] = False
        id_dict['stripsimg'] = True
        id_dict['acqNum'] = ''

        nstrips_arr = [32, 33, 34, 35, 36, 37, 39, 40, 42]
        id_dict['nstrips'] = nstrips_arr[id_dict['overlap']/2]

        return id_dict


    def build_acq1_dict(self):
        '''
        Build the ACQ1 step:
        nx, ny = 128
        nramps = 6
        ACQ1 = repeat 6x(reset, drop, read, drop, read)
        '''
        acq1_dict = {}
        acq1_dict['imgsize'] = 128
        acq1_dict['nramps'] = 6

        acq1_dict['tcds'] = TCDSACQ1
        acq1_dict['step'] = 'ACQ1'

        acq1_dict['xarr'] = self.xarr[0]
        acq1_dict['yarr'] = self.yarr[0]

        acq1_dict['cdsimg'] = True
        acq1_dict['stripsimg'] = False
        acq1_dict['acqNum'] = 1

        return acq1_dict


    def build_acq2_dict(self):
        '''
        Build the ACQ2 step:
        nx, ny = 32
        nramps = 5
        ACQ2 = repeat 5x(reset, drop, read, drop, drop, drop, read, drop)
        '''
        acq2_dict = {}
        acq2_dict['imgsize'] = 32
        acq2_dict['nramps'] = 5

        acq2_dict['tcds'] = TCDSACQ2
        acq2_dict['step'] = 'ACQ2'

        acq2_dict['xarr'] = self.xarr[0]
        acq2_dict['yarr'] = self.yarr[0]

        acq2_dict['cdsimg'] = True
        acq2_dict['stripsimg'] = False
        acq2_dict['acqNum'] = 2

        return acq2_dict


    def build_trk_dict(self, nramps):
        '''
        Build the TRK step:
        nx, ny = 32
        nramps = a large number (i.e. 5000+)
        '''
        trk_dict = {}

        trk_dict['imgsize'] = 32
        trk_dict['nramps'] = nramps

        trk_dict['tcds'] = TCDSTRK
        trk_dict['step'] = 'TRK'

        trk_dict['xarr'] = self.xarr[0]
        trk_dict['yarr'] = self.yarr[0]

        trk_dict['cdsimg'] = False
        trk_dict['stripsimg'] = False
        trk_dict['acqNum'] = ''

        return trk_dict


    def build_lostrk_dict(self):
        '''
        Build the Line of Sight - TRK step:
        nx, ny = 43
        nramps = a large number (i.e. 5000+)
        '''
        lostrk_dict = {}

        lostrk_dict['imgsize'] = 43
        lostrk_dict['nramps'] = 1

        lostrk_dict['tcds'] = TCDSTRK
        lostrk_dict['step'] = 'LOSTRK'

        lostrk_dict['xarr'] = self.xarr[0]
        lostrk_dict['yarr'] = self.yarr[0]

        lostrk_dict['cdsimg'] = False
        lostrk_dict['stripsimg'] = False
        lostrk_dict['acqNum'] = ''

        return lostrk_dict

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
        xlow = int(xcoord) - (imgsize / 2 + 1)
        ylow = int(ycoord) - (imgsize / 2 + 1)
    else:
        xlow = int(xcoord) - imgsize / 2
        ylow = int(ycoord) - imgsize / 2

    xhigh = int(xcoord) + imgsize / 2
    yhigh = int(ycoord) + imgsize / 2

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
