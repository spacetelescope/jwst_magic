# STDLIB
import os

# Third Party
from astropy.io import fits
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import numpy as np

# LOCAL
import rptoia
import getbias
import utils
from mkproc import Mkproc
import log
import select_psfs

# DEFINE ALL NECESSARY CONSTANTS
TCDSID = 0.338    # neil 18 may 12
TCDSACQ1 = 0.3612 # 20 apr 17 keira # 1cds time = 2*(128*128*10.e-6) = 0.32768s
TCDSACQ2 = 0.0516 # 20 apr 17 keira brooks # 2cds time = 4*(32*32*10e-6) = 0.04096s
TCDSTRK = 0.0256  # neil 18 may 12
TCDSFG = 0.0512   # neil 18 may 12

LOCAL_PATH = os.path.dirname(os.path.realpath(__file__))
# Location of all necessary *.fits files (newmagicHdrImg, bias0, etc)
DATA_PATH = os.path.join(LOCAL_PATH, 'data')

class FGS(object):
    '''
    Creates an FGS simulation object for ID, ACQ, and/or TRK stages to be used
    with DHAS.
    '''
    def __init__(self, im, guider, root, out_dir=None, data_path=None,
                 guide_star_coords=None, reg_file=None, overlap=0, biaszeropt=True,
                 biasktc=True, biasped=True, poissonnoise=True, background=False,
                 scdrift=False):

        self.nreads = 2

        # Noise to add
        self.biaszeropt = biaszeropt
        self.biasktc = biasktc
        self.biasped = biasped
        self.poissonnoise = poissonnoise
        self.background = background

        # Detector motion
        self.scdrift = scdrift ## Do I want this here???

        # Practical things
        self.guider = guider
        self.root = root

        ## DEFINE ALL THINGS PATHS
        # Define output directory
        if out_dir is None:
            self.out_dir = os.path.join(LOCAL_PATH, 'out', self.root)
        else:
            self.out_dir = out_dir

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
        self.input_im = correct_image(self.input_im)
        log.info('Max of input image: {}'.format(np.max(self.input_im)))

        self.get_coords_and_counts(reg_file=reg_file)
        if guide_star_coords is None:
            self.get_guide_star_coords(gs_ind=0)


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

    def get_guide_star_coords(self, gs_ind=0):
        '''
        From the list of coordinates, pull out the guide star coordinates using
        the guide star index (gs_ind)
        '''
        self.xgs = self.xarr[gs_ind]
        self.ygs = self.yarr[gs_ind]
        self.countrategs = self.countrate[gs_ind]
        log.info('Coordinates of Guide Star: x={0}, y={1}'.format(self.xgs, self.ygs))

    def create_noisy_sky(self, bias, time_normed_im, poissonnoise=True, acqNum=None):

        '''
        Create a noisy sky image for ID and ACQ steps.
        There is only an added poissonfactor for ACQ2.
        PoissonNoise should always be set to True for ID.
        '''

        if self.step == 'ACQ':
            time_normed_im = 2. * time_normed_im  # Is this still necessary?

        image = np.copy(bias)

        if poissonnoise:
            # Add (increasing) poisson noise + data to each read in a ramp
            for ireads in range(self.nreads):
                image[ireads::(self.nreads)] += (ireads + 1) * np.random.poisson(time_normed_im)

        else:
            # Add (increasing) data to each read in a ramp
            for ireads in range(self.nreads):
                image[ireads::(self.nreads)] += (ireads + 1) * time_normed_im

        image = correct_image(image)
        return image

    def create_strips(self, arr):
        '''
        Create the ID strips fits image to be passed into DHAS
        '''
        strips = np.zeros((self.nstrips * self.nz, self.h, self.nx))
        nn = 0
        for i in range(self.nstrips):
            for iz in range(self.nz):
                y1 = i * (self.h - self.overlap) + self.yoffset
                y2 = y1 + self.h
                strips[nn] = arr[iz, y1:y2]
                nn += 1

        # Make sure the data is between 0 and 65,000 counts and are finite numbers
        strips = correct_image(strips)
        strips = np.uint16(strips)

        return strips

    # Write out files
    def write_stc(self, filename, xarr, yarr, countrate):
        """
        Write out stc files using offset, rotated catalog
        """
        xi, yi = rptoia.rptoia(xarr, yarr, self.guider)

        try:
            len(xi)
        except TypeError:
            xi = [xi]
            yi = [yi]
            countrate = [countrate]

        inum = np.arange(len(xi))
        data = np.array([inum, xi, yi, countrate]).T
        utils.write_to_file(filename, data, fmt=['%d', '%f', '%f', '%e'])
        print("Successfully wrote: {}".format(filename))

    def write_cat(self, filename, xarr, yarr):
        '''
        Write out star catalog in real pixs
        '''
        coords = np.array([xarr, yarr]).T
        try:
            utils.write_to_file(filename, coords, fmt=['%d', '%d'])
        except AttributeError:
            utils.write_to_file(filename, coords)
        print("Successfully wrote: {}".format(filename))

    # Put it all together
    def setup_step(self, nx, ny, nramps, tcds, step, yoffset=0, h=64, overlap=8,
                   nstrips=None):
        '''
        Set up attributes and basic arrays for ID, acquisitions 1 or 2, and TRK:
        ID:
        nx, ny = 2048 (or size of truth image)
        nramps = 2
        nstrips = [32,33,34,35,36,37,39,40,42] for overlap = [0,2,4,6,8,10,12,14,16,18]

        Acq 1:
        nx, ny = 128
        nramps = 6
        ACQ1 = repeat 6x(reset, drop, read, drop, read)

        Acq 2:
        nx, ny = 32
        nramps = 5
        ACQ2 = repeat 5x(reset, drop, read, drop, drop, drop, read, drop)

        TRK:
        nx, ny = 32
        nramps = a large number (i.e. 5000+)
        '''

        self.nx = nx
        self.ny = ny
        self.nramps = nramps
        self.nz = self.nreads * self.nramps

        self.step = step

        if self.step == 'ID':
            self.yoffset = yoffset
            self.h = h
            self.overlap = overlap
            if nstrips is None:
                nstrips_arr = [32, 33, 34, 35, 36, 37, 39, 40, 42]
                self.nstrips = nstrips_arr[overlap/2]
            else:
                self.nstrips = nstrips

        else:
            self.input_im = create_im_subarray(self.input_im, self.xgs,
                                               self.ygs, self.nx, self.ny)

        self.time_normed_im = self.input_im * tcds

    def create_arrays(self, x, y, acqNum=None, cds=True):
        '''
        Create the necessary arrays:
        bias
        arr (noisy sky image with correct reads and ramps)

        NOTE: Include an acquisition number ("acqNum") for acquisition
        NOTE: For acquisition, x,y should be the coordinates of the guide star
        '''



        self.bias = getbias.getbias(self.guider, self.xgs, self.ygs, self.nreads,
                                    self.nramps,
                                    self.nx, self.ny, bzp=self.biaszeropt,
                                    bktc=self.biasktc, bp=self.biasped)
        self.image = self.create_noisy_sky(self.bias, self.time_normed_im,
                                           self.poissonnoise, acqNum)

        if cds:
            self.cds = create_cds_image(self.image)

        if self.step == 'ID':
            self.strips = self.create_strips(self.image)

    def write_out_files(self, x, y, countrate, acqNum=''):
        '''
        Create **all** the files and images needs for ID & ACQ
        Requires: time_normed_im, x,y,countrate,bias,idarr

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
        <name>_G<guider>_<step><acqnum>.cat :   list of x and y coordinates of the guide
                                            and reference stars
        <name>_G<guider><step><acqnum>.fits :   fits cube of the sky with added bias and noise
                                            to simulate the read and ramp cycle for ID
        '''
        ## STScI only files - mostly just for quick checks of the data
        log.info('Baseline {}'.format(np.max(self.time_normed_im)))

        # Sky imge
        filename_sky = os.path.join(self.out_dir,
                                    'stsci',
                                    '{}_G{}_{}{}sky.fits'.format(self.root,
                                                                 self.guider,
                                                                 self.step,
                                                                 acqNum))
        utils.write_fits(filename_sky, self.time_normed_im)

        # STC files using offset, rotated catalog
        filename_stc = os.path.join(self.out_dir,
                                    'stsci',
                                    '{}_G{}_{}{}.stc'.format(self.root,
                                                             self.guider,
                                                             self.step,
                                                             acqNum))
        self.write_stc(filename_stc, x, y, countrate)

        # Bias image
        filename_bias = os.path.join(self.out_dir,
                                     'stsci',
                                     '{}_G{}_{}{}bias.fits'.format(self.root,
                                                                   self.guider,
                                                                   self.step,
                                                                   acqNum))
        utils.write_fits(filename_bias, self.bias)

        # Create CDS image
        filename_cds = os.path.join(self.out_dir,
                                    'stsci',
                                    '{}_G{}_{}{}cds.fits'.format(self.root,
                                                                 self.guider,
                                                                 self.step,
                                                                 acqNum))
        utils.write_fits(filename_cds, self.cds)

        if self.step == 'ID':
            ## STScI only files - mostly just for quick checks of the data
            # Star catalog in real pixs
            filename_starcat = os.path.join(self.out_dir,
                                            'stsci',
                                            '{}_G{}_{}.gssscat'.format(self.root,
                                                                       self.guider,
                                                                       self.step))
            self.write_cat(filename_starcat, self.xarr, self.yarr)

            # Create ff fits file
            filename_ff = os.path.join(self.out_dir,
                                       'stsci',
                                       '{}_G{}_{}ff.fits'.format(self.root,
                                                                 self.guider,
                                                                 self.step))
            utils.write_fits(filename_ff, np.uint16(self.image))

            ## DHAS file
            # Extract strips from ff img
            filename_id_strips = os.path.join(self.out_dir,
                                              'dhas',
                                              '{}_G{}_{}strips.fits'.format(self.root,
                                                                            self.guider,
                                                                            self.step))
            # Write to strips to fits file
            filename_hdr = os.path.join(DATA_PATH,
                                        'newG{}magicHdrImg.fits'.format(self.guider))
            dummy_data, hdr0 = fits.getdata(filename_hdr, header=True)
            utils.write_fits(filename_id_strips, self.strips, header=hdr0)

            ## Gound system file
            convert_fits_to_dat(filename_id_strips, self.step,
                                os.path.join(self.out_dir, 'ground_system'))

        elif self.step == 'ACQ':
            ## STScI only files - mostly just for quick checks of the data
            # star catalog in real pixs
            filename_starcat = os.path.join(self.out_dir,
                                            'stsci',
                                            '{}_G{}_{}{}.cat'.format(self.root,
                                                                     self.guider,
                                                                     self.step,
                                                                     acqNum))
            self.write_cat(filename_starcat, self.xgs, self.ygs)

            ## DHAS file
            # Noisy sky acquisition fits images
            filename_noisy_sky = os.path.join(self.out_dir,
                                              'dhas',
                                              '{}_G{}_{}{}.fits'.format(self.root,
                                                                        self.guider,
                                                                        self.step,
                                                                        acqNum))
            utils.write_fits(filename_noisy_sky, np.uint16(self.image))

            ## Gound system file
            convert_fits_to_dat(filename_noisy_sky, self.step,
                                os.path.join(self.out_dir, 'ground_system'))

        else:
            pass

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
    fl = data.flatten()

    if (obsmode == 'PSF') or (obsmode == 'TRK'):
        # ascii float format
        f = '{:16.7e} '

    elif (obsmode == 'ID') or (obsmode == 'ACQ1') or (obsmode == 'ACQ2') or \
         (obsmode == 'ACQ') or (obsmode == 'CAL'):
        # ascii hex dat format
        f = '{:04X} '

    else:
        log.error("Observation mode not recognized. Returning.")

    with open(os.path.join(out_dir, outfile), 'w') as file_out:
        for d in fl.astype(np.uint16):
            file_out.write(f.format(d))

    print("Successfully wrote: {}".format(os.path.join(out_dir, outfile)))
    return
#-------------------------------------------------------------------------------
def create_cds_image(arr):
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

def create_im_subarray(image, x, y, nx, ny, show_fig=False):
    '''
    Based on the array size given by nx and ny, created a subarray around
    the guide star.
    '''

    img = np.copy(image)
    if (nx % 2 == 1):
        x1 = int(x) - (nx / 2 + 1)
        y1 = int(y) - (ny / 2 + 1)
    else:
        x1 = int(x) - nx / 2
        y1 = int(y) - ny / 2

    x2 = int(x) + nx / 2
    y2 = int(y) + ny / 2

    img = img[y1:y2, x1:x2]

    if show_fig:
        plt.figure()
        plt.imshow(img, cmap='Greys_r')
        plt.show()

    return img

def correct_image(image, upper_threshold=65000, upper_limit=65000):
    '''
    Correct image for negative and saturated pixels
    '''
    img = np.copy(image)
    img[img < 0] = 0            # neg pixs -> 0
    img[img >= upper_threshold] = upper_limit    # sat'd pixs -> 65K
    img[np.isfinite(img) == 0] = 0

    return img


def add_background(array, nx, ny, nz):
    '''
    Add background to array
    '''
    try:
        array += 500 + 10. * np.random.standard_normal((nz, ny, nx))
    except NameError:
        array = 500 + 10. * np.random.standard_normal((nz, ny, nx))

    return array

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
def run_id(image, guider, root, out_dir=None, interactive=False):

    '''
    Create an ID object and create all necessary files to fun the ID simulation
    in DHAS. Also creates CECIL proc file. Returns ID object.
    '''
    # ID
    id0 = FGS(image, guider, root, out_dir)
    id0.setup_step(2048, 2048, 2, tcds=TCDSID, step='ID')
    log.info("Step: {}".format(id0.step))
    id0.create_arrays(id0.xarr, id0.yarr, acqNum=None)
    id0.write_out_files(id0.xarr, id0.yarr, id0.countrate)

    # Make CECIL proc file
    Mkproc(guider, root, id0.xarr, id0.yarr, id0.countrate, step='ID',
           out_dir=out_dir)

    if interactive:
        return id0


def run_acq(image, guider, root, out_dir=None, interactive=False):
    '''
    Creates two ACQ objects (for acquisitions #1 and #2) and create all necessary
    files to fun the ACQ simulation in DHAS. Also creates CECIL proc file.
    Returns ACQ1 and ACQ2 objects.
    '''
    # Acquisition #1
    acq1 = FGS(image, guider, root, out_dir)
    acq1.setup_step(nx=128, ny=128, nramps=6, tcds=TCDSACQ1, step='ACQ')
    log.info("Step: {}".format(acq1.step))
    acq1.create_arrays(acq1.xgs, acq1.ygs, acqNum=1)
    acq1.write_out_files(acq1.xgs - (acq1.nx / 2.), acq1.ygs - (acq1.ny / 2.),
                         acq1.countrategs, acqNum=1)

    # Acquistion #2
    acq2 = FGS(image, guider, root, out_dir)
    acq2.setup_step(nx=32, ny=32, nramps=5, tcds=TCDSACQ2, step='ACQ')
    acq2.create_arrays(acq2.xgs, acq2.ygs, acqNum=2)
    acq2.write_out_files(acq2.xgs - (acq2.nx / 2.), acq2.ygs - (acq2.ny / 2.),

                         acq2.countrategs, acqNum=2)

    # Make CECIL proc file
    Mkproc(guider, root, acq1.xgs, acq1.ygs, acq1.countrategs, step='ACQ',
           out_dir=out_dir)

    if interactive:
        return acq1, acq2

def run_trk(image, guider, root, num_frames, out_dir=None, jitter=True, interactive=False):
    '''
    Create a TRK object and create all necessary files to fun the TRK simulation
    in DHAS. Returns TRK object.
    '''
    if num_frames is None:
        num_frames = 5000
    trk = FGS(image, guider, root, out_dir)
    trk.setup_step(nx=32, ny=32, nramps=num_frames, tcds=TCDSTRK, step='TRK')
    log.info("Step: {}".format(trk.step))
    trk.create_arrays(trk.xgs, trk.ygs, cds=False)

    if jitter:
        trk.image = add_jitter(trk.image, total_shift=1)

    filename_noisy_sky = os.path.join(trk.out_dir,
                                      'dhas', '{}_G{}_{}.fits'.format(trk.root,
                                                                      trk.guider,
                                                                      trk.step))
    utils.write_fits(filename_noisy_sky, np.uint16(trk.image))

    if interactive:
        return trk


def create_lostrk(image, guider, root, nx, ny, out_dir=None):
    trk = FGS(image, guider, root, out_dir)
    trk.nreads = 1  #want only single frame image
    trk.setup_step(nx, ny, nramps=1, tcds=TCDSTRK, step='TRK')
    trk.create_arrays(trk.xgs, trk.ygs, cds=False)

    # This is a ground system file, but needs to be in .dat format
    filename_noisy_sky = os.path.join(trk.out_dir,
                                      'dhas',
                                      '{}_G{}_LOS{}.fits'.format(trk.root,
                                                                 trk.guider,
                                                                 trk.step))
    utils.write_fits(filename_noisy_sky, trk.image)

    ## Gound system file
    convert_fits_to_dat(filename_noisy_sky, trk.step, os.path.join(trk.out_dir,
                                                                   'ground_system'))

def run_all(im, guider, root, num_frames=None, out_dir=None,
            jitter=True):
    # Not interatctive!!
    run_id(im, guider, root, out_dir, nref)
    run_acq(im, guider, root, out_dir)
    run_trk(im, guider, root, num_frames, out_dir, jitter)
