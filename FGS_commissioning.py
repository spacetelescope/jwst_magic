import argparse
from astropy.io import fits
import numpy as np
import csv
import os
from skimage.filters import threshold_otsu
from scipy import ndimage
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

#local
import grptoia
import getbias
import utils
from mkproc import mkproc


class FGS(object):
    '''
    Creates an FGS simulation object for ID, ACQ, and/or TRK stages to be used
    with DHAS.
    '''
    def __init__(self, im, guider, root, out_dir=None, data_path=None,
                guide_star_coords=None, reg_file=None, overlap=0, biasZeroPt=False,
                biasKTC=False, biasPed=False, poissonNoise=False, background=False,
                SCdrift=False):


        # DEFINE ALL NECESSARY CONSTANTS
        self.tcdsID = 0.338    # neil 18 may 12
        self.tcdsACQ1 = 0.3612 # 20 apr 17 keira # 1cds time = 2*(128*128*10.e-6) = 0.32768s
        self.tcdsACQ2 = 0.0516 # 20 apr 17 keira brooks # 2cds time = 4*(32*32*10e-6) = 0.04096s
        self.tcdsTRK = 0.0256  # neil 18 may 12
        self.tcdsFG = 0.0512   # neil 18 may 12

        self.nreads = 2

        # Noise to add
        self.biasZeroPt = biasZeroPt
        self.biasKTC = biasKTC
        self.biasPed = biasPed
        self.poissonNoise = poissonNoise
        self.background = background

        # Detector motion
        self.SCdrift = SCdrift ## Do I want this here???

        # Practical things
        self.guider = guider
        self.root = root


        ## DEFINE ALL THINGS PATHS
        self.local_path = os.path.dirname(os.path.realpath(__file__))
        # Define location of all necessary *.fits files (newmagicHdrImg, bias0, etc)
        if data_path is None:
            self.data_path = os.path.join(self.local_path,'data')
        else:
            self.data_path = data_path
        # Define output directory
        if out_dir is None:
            self.out_dir = os.path.join(self.local_path,'out', self.root,'dhas_files')
        else:
            self.out_dir = out_dir

        if not os.path.exists(self.out_dir):
            os.makedirs(self.out_dir)


        ## READ IN IMAGE
        if isinstance(im,str):
            hdu=fits.open(im) #*_bin_norm from FGS_bin_tool
            self.input_im = np.uint16(hdu[0].data)
        else:
            self.input_im = np.uint16(im)

        # Correct for negative, saturated pixels and other nonesense
        self.input_im = correct_image(self.input_im)
        print('Max of input image: {}'.format(np.max(self.input_im)))

        self.get_coords_and_counts(reg_file=reg_file)
        if guide_star_coords is None:
            self.get_guide_star_coords(gs_ind=0)


    ### Guide star and reference star coordinates and countrates
    def get_coords_and_counts(self,reg_file=None):
        '''
        Get coordinate information of guide star and reference stars
        '''

        if reg_file is None:
            reg_file = os.path.join(os.path.join(self.local_path,'out', self.root,
                                    '{0}_G{1}_regfile.txt'.format(self.root,self.guider)))#'{}.reg'.format(root)
        print("Using {} as the reg file".format(reg_file))
        if reg_file.endswith('reg'):
            self.x,self.y = np.loadtxt(reg_file)
            self.countrate = utils.get_countrate(self.x,self.y,self.input_im)
        else:
            self.y,self.x,self.countrate = np.loadtxt(reg_file,delimiter=' ',skiprows=1).T #includes y, x, coords and countrate


    def get_guide_star_coords(self, gs_ind=0):
        '''
        From the list of coordinates, pull out the guide star coordinates using
        the guide star index (gs_ind)
        '''
        self.xgs=self.x[gs_ind]
        self.ygs=self.y[gs_ind]
        self.countrategs = self.countrate[gs_ind]
        print('Coordinates of Guide Star: x={0}, y={1}'.format(self.xgs,self.ygs))


    def create_noisy_sky(self,bias,sky,poissonNoise=True,acqNum=None):
        '''
        Create a noisy sky image for ID and ACQ steps.
        There is only an added poissonfactor for ACQ2.
        PoissonNoise should always be set to True for ID.
        '''
        ### 4/24 NEED TO CHECK THIS SECTION
        if acqNum == 2:
            poissonfactor = 2.
        else:
            poissonfactor = 1.

        if self.step == 'ACQ':
            sky = 2.*sky

        im = np.copy(bias)
        if poissonNoise:
            for ireads in range(self.nreads):
                im[ireads::(self.nreads)]+=(ireads+1) * 0.25 * np.random.poisson(poissonfactor*sky)
        else:
            for ireads in range(self.nreads):
                im[ireads::(self.nreads)]+=(ireads+1) * poissonfactor * sky

        im = correct_image(im)

        return im

    def create_cds_image(self,arr):
        '''
        Create CDS image
        '''
        return arr[1::2]-arr[:-1:2]


    def create_strips(self, arr):
        '''
        Create the ID strips fits image to be passed into DHAS
        '''
        strips=np.zeros((self.nstrips*self.nz,self.h,self.nx))
        nn=0
        for i in range(self.nstrips):
           for iz in range(self.nz):
               y1 = i*(self.h-self.overlap)+self.yoffset
               y2 = y1+self.h
               strips[nn,:,:]=arr[iz,y1:y2,:]
               nn+=1

        # Make sure the data is between 0 and 65,000 counts and are finite numbers
        strips= correct_image(strips)
        strips=np.uint16(strips)

        return strips

    ### Write out files
    def write_stc(self,filename,x,y,countrate):
        """
        Write out stc files using offset, rotated catalog
        """
        if self.guider == 1:
            xi, yi = grptoia.g1RPtoIA(x,y)
        elif self.guider == 2:
            xi, yi = grptoia.g2RPtoIA(x,y)

        try:
            len(xi)
        except TypeError:
            xi = [xi]
            yi = [yi]
            countrate = [countrate]

        inum = np.arange(len(xi))
        write_file = open(filename,'w+')
        data = np.array([inum,xi,yi,countrate]).T
        np.savetxt(write_file,data,fmt=['%d','%f','%f','%e'])
        write_file.close()
        print("Successfully wrote: {}".format(filename))


    def write_cat(self,filename,x,y):
        '''
        Write out star catalog in real pixs
        '''
        write_file = open(filename,'w+')
        coords = np.array([self.x,self.y]).T
        try:
            np.savetxt(write_file,coords,fmt=['%d','%d'])
        except AttributeError:
            np.savetxt(write_file,coords)
        write_file.close()
        print("Successfully wrote: {}".format(filename))




    ### Put it all together
    def setup_step(self, nx, ny, nramps, tcds, step, yoffset=0, h=64, overlap=8,
                   nstrips=None):
        '''
        Set up attributes and basic arrays for ID, or acquisitions 1 or 2:
        ID:
        nx, ny = 2048 (or size of truth image)
        nramps = 2
        nstrips = [32,33,34,35,36,37,39,40,42] for overlap =  [0,2,4,6,8,10,12,14,16,18]

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
                nstrips_arr = [32,33,34,35,36,37,39,40,42]
                self.nstrips = nstrips_arr[overlap/2]
            else:
                self.nstrips = nstrips

        else:
            self.input_im = create_im_subarray(self.input_im,self.xgs,self.ygs,self.nx,self.ny)


        self.sky = self.input_im * tcds


    def create_arrays(self, x, y, acqNum=None, cds=True):
        '''
        Create the necessary arrays:
        bias
        arr (noisy sky image with correct reads and ramps)

        NOTE: Include an acquisition number ("acqNum") for acquisition
        NOTE: For acquisition, x,y should be the coordinates of the guide star
        '''
        self.bias = getbias.getbias(self.guider, x, y, self.nreads,
                                    self.nramps, self.nx, self.ny,
                                    self.biasZeroPt, self.biasKTC, self.biasPed,
                                    data_path = self.data_path)
        self.image = self.create_noisy_sky(self.bias, self.sky, self.poissonNoise, acqNum)

        if cds:
            self.cds = self.create_cds_image(self.image)

        if self.step == 'ID':
            self.strips = self.create_strips(self.image)


    def write_out_files(self,x,y,countrate,acqNum=''):
        '''
        Create **all** the files and images needs for ID & ACQ
        Requires: sky, x,y,countrate,bias,idarr

        ID & ACQ:
        <name>_G<guider>_<step>sky.fits :   fits file of they sky array

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

        print('baseline {}'.format(np.max(self.sky)))
        filename_sky = os.path.join(self.out_dir,'{}_G{}_{}{}sky.fits'.format(self.root,self.guider,self.step,acqNum))
        utils.write_fits(filename_sky, self.sky)


        # Write stc files using offset, rotated catalog
        filename_stc = os.path.join(self.out_dir,'{}_G{}_{}{}.stc'.format(self.root,self.guider,self.step,acqNum))
        self.write_stc(filename_stc,x,y,countrate)

        # Bias image
        filename_bias = os.path.join(self.out_dir,'{}_G{}_{}{}bias.fits'.format(self.root,self.guider,self.step,acqNum))
        utils.write_fits(filename_bias,self.bias)

        # Create CDS image
        filename_cds = os.path.join(self.out_dir,'{}_G{}_{}{}cds.fits'.format(self.root,self.guider,self.step,acqNum))
        utils.write_fits(filename_cds, self.cds)

        if self.step == 'ID':
            # Star catalog in real pixs
            filename_starcat = os.path.join(self.out_dir,'{}_G{}_{}.gssscat'.format(self.root,self.guider,self.step))
            self.write_cat(filename_starcat,self.x,self.y)

            # Create ff fits file
            filename_ff = os.path.join(self.out_dir,'{}_G{}_{}ff.fits'.format(self.root,self.guider,self.step))
            utils.write_fits(filename_ff,np.uint16(self.image))

            # Extract strips from ff img
            filename_id_strips = os.path.join(self.out_dir,'{}_G{}_{}strips.fits'.format(self.root,self.guider,self.step))
            # Write to strips to fits file
            hdr0 = fits.open(os.path.join(self.data_path,'newG{}magicHdrImg.fits'.format(self.guider)))[0].header
            utils.write_fits(filename_id_strips,self.strips,header=hdr0)

        elif self.step == 'ACQ':
            # star catalog in real pixs
            filename_starcat = os.path.join(self.out_dir,'{}_G{}_{}{}.cat'.format(self.root,self.guider,self.step,acqNum))
            self.write_cat(filename_starcat,self.xgs,self.ygs)

            # Noisy sky acquisition fits images
            filename_noisy_sky = os.path.join(self.out_dir,'{}_G{}_{}{}.fits'.format(self.root,self.guider,self.step,acqNum))
            utils.write_fits(filename_noisy_sky,np.uint16(self.image))

        else:
            pass

#-------------------------------------------------------------------------------
def display(image, ind=0, vmin=None, vmax=None, cmap='Greys_r', x=None, y=None,
            add_coords=False):
    '''
    Display an image array. If the array is a cube, specify the index to
    look at using 'ind'. If you want to add a scatter plot of PSF centers,
    use add_coords=True and pass in x and y (Note: This only works for full
    frame images for the time being).
    '''
    if image.ndim == 3:
        im = image[ind]
    else:
        im = np.copy(image)

    plt.clf()
    plt.imshow(im,cmap=cmap,norm=LogNorm(),vmin=vmin,vmax=vmax)

    if add_coords:
        plt.scatter(x,y,color='white')
    else:
        plt.colorbar()

    plt.show()

### Arrays and array manipulation
def add_jitter(cube,x,y,nx,ny,total_shift=3):
    '''
    Add random single pixel jitter

    VERY rudimentary. Uses np.roll so images look kind of funny.
    '''
    cube2 = np.zeros_like(cube)
    for i in range(len(cube)):
        shift = np.random.randint(total_shift+1) # This will generate a random integer up to the total_shift
        cube2[i] = np.roll(cube[i],shift)

    return cube2

def create_im_subarray(image,x,y,nx,ny,show_fig=False):
    '''
    Based on the array size given by nx and ny, created a subarray around
    the guide star.
    '''
    im = np.copy(image)
    if (nx % 2 == 1):
        x1 = int(x)-(nx/2+1)
        y1 = int(y)-(ny/2+1)
    else:
        x1 = int(x)-nx/2
        y1 = int(y)-ny/2

    x2 = int(x)+nx/2
    y2 = int(y)+ny/2

    im = im[y1:y2,x1:x2]

    if show_fig:
        plt.figure()
        plt.imshow(im,cmap='Greys_r')
        plt.show()

    return im

def correct_image(image, upper_threshold=65000, upper_limit=65000):
    '''
    Correct image for negative and saturated pixels
    '''
    im = np.copy(image)
    im[im < 0] = 0            # neg pixs -> 0
    im[im >= upper_threshold] = upper_limit    # sat'd pixs -> 65K
    im[np.isfinite(im) == 0] = 0

    return im


def add_background(array,nx,ny,nz):
    '''
    Add background to array
    '''
    try:
        array += 500+10.*np.random.standard_normal((nz,ny,nx))
    except NameError:
        array = 500+10.*np.random.standard_normal((nz,ny,nx))

    return array

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
def run_ID(im, guider, root, out_dir=None, template_path=None, nref=10, interactive=False):
    '''
    Create an ID object and create all necessary files to fun the ID simulation
    in DHAS. Also creates CECIL proc file. Returns ID object.
    '''
    # ID
    id0 = FGS(im, guider, root, out_dir)
    id0.setup_step(2048, 2048, 2, tcds=id0.tcdsID, step='ID')
    print(id0.step)
    id0.create_arrays(id0.x,id0.y,acqNum=None)
    id0.write_out_files(id0.x,id0.y,id0.countrate)

    # Make CECIL proc file
    mkproc(guider, root, id0.x, id0.y, id0.countrate, step='ID', nref=nref, out_dir=out_dir,
           template_path=template_path)

    if interactive:
        return id0


def run_ACQ(im, guider, root, out_dir=None, template_path=None, interactive=False):
    '''
    Creates two ACQ objects (for acquisitions #1 and #2) and create all necessary
    files to fun the ACQ simulation in DHAS. Also creates CECIL proc file.
    Returns ACQ1 and ACQ2 objects.
    '''
    # Acquisition #1
    acq1 = FGS(im, guider, root, out_dir)
    acq1.setup_step(nx=128, ny=128, nramps=6, tcds=acq1.tcdsACQ1,step='ACQ')
    print(acq1.step)
    acq1.create_arrays(acq1.xgs,acq1.ygs,acqNum=1)
    acq1.write_out_files(acq1.xgs-(acq1.nx/2.),acq1.ygs-(acq1.ny/2.),acq1.countrategs,acqNum=1)

    # Acquistion #2
    acq2 = FGS(im, guider, root, out_dir)
    acq2.setup_step(nx=32, ny=32, nramps=5, tcds=acq2.tcdsACQ2,step='ACQ')
    acq2.create_arrays(acq2.xgs,acq2.ygs,acqNum=2)
    acq2.write_out_files(acq2.xgs-(acq2.nx/2.),acq2.ygs-(acq2.ny/2.),acq2.countrategs,acqNum=2)

    # Make CECIL proc file
    mkproc(guider, root, acq1.xgs, acq1.ygs, acq1.countrategs, step='ACQ', nref=None,
           out_dir=out_dir, template_path=template_path)

    if interactive:
        return acq1, acq2

def run_TRK(im, guider, root, num_frames, out_dir=None, jitter=True, interactive=False):
    '''
    Create a TRK object and create all necessary files to fun the TRK simulation
    in DHAS. Returns TRK object.
    '''
    if num_frames is None:
        num_frames = 5000
    trk = FGS(im, guider, root, out_dir)
    trk.setup_step(nx=32, ny=32, nramps=num_frames, tcds=trk.tcdsTRK, step='TRK')
    print(trk.step)
    trk.create_arrays(trk.xgs,trk.ygs,cds=False)

    if jitter:
        trk.image = add_jitter(trk.image,trk.xgs,trk.ygs,trk.nx,trk.ny,total_shift=1)

    filename_noisy_sky = os.path.join(trk.out_dir,'{}_G{}_{}.fits'.format(trk.root,trk.guider,trk.step))
    utils.write_fits(filename_noisy_sky,np.uint16(trk.image))

    if interactive:
        return trk

def create_LOSTRK(im, guider, root, nx, ny, out_dir=None,interactive=False):
    trk = FGS(im, guider, root, out_dir)
    trk.nreads = 1 #want only single frame image
    trk.setup_step(nx, ny, nramps=1, tcds=trk.tcdsTRK, step='TRK')
    trk.create_arrays(trk.xgs,trk.ygs,cds=False)

    filename_noisy_sky = os.path.join(trk.out_dir,'{}_G{}_LOS{}.fits'.format(trk.root,trk.guider,trk.step))
    utils.write_fits(filename_noisy_sky,trk.image)

def run_all(im, guider, root, num_frames=None, out_dir=None, template_path=None, jitter=True,nref=10):
    # Not interatctive!!
    run_ID(im, guider, root, out_dir, template_path,nref)
    run_ACQ(im, guider, root, out_dir, template_path)
    run_TRK(im, guider, root, num_frames, out_dir, jitter)
