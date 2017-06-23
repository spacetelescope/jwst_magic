import numpy as np
np.seterr(all='raise')
import os
import csv
from astropy.io import fits
import scipy.signal
from glob import glob

import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

from skimage.filters import threshold_otsu
from scipy import ndimage

#Local
from mkproc import mkproc
import utils
import counts_to_jmag

'''
Python Img Re-Binning Tool for FGS comissioning data

Input NIRCam images, rebin to FGS platescale and renormalize to desired FGS

NIRCam FOV: 1.09' by 1.09'
NIRCam pixel scale: 0.032"/pixel

FGS FOV: 2.4' by 2.4'
FGS pixel scale: 0.070"/pixel


For development of this code, used mock-up NIRCam files from Ball Aerospace's ITM tool,
simulating Global Alignment using the short wavelength channel

(Sherie was told they include detector noise and bias offsets between the readout channels)
'''


#-------------------------------------------------------------------------------
def bad_pixel_correction(data,BP_thresh):
    '''Finds and smooths out bad pixels with a median filter'''
    #apply median filter
    smooth = scipy.signal.medfilt(data,3)

    #set negative values to zero
    j = smooth.copy()
    j[j < 0] = 0

    #create array of zeros, for the bad pixel mask
    bpmask = np.zeros(data.shape)

    #difference between image and smoothed image; leaves the background behind so we can filter out the bad pixels
    delta = data-smooth

    #Locating the bad pixels. If there are still bpix in masked image, fiddle with delta threshold
    j = np.where(delta > BP_thresh)

    #using location of the bad pixels, replace the bpix value with median value of the smoothed image
    #also get rid of any negative numbers
    data[j] = np.median(smooth)
    data[data<0]=0

    #recast as unsigned integers
    data = np.int_(data)

    #clip any over saturated/hot pixels left, replace with integer form of median value of smoothed image
    data[data > 50000] = np.int_(np.median(smooth))

    return data


def rotate_nircam_image(im,fgs_guider,header,nircam_mod):
    '''
    Given NIRCAM module A or B (given by the header in your original NIRCAM image),
    rotate/flip to put in correct orientation for FGS 1 and 2.
    '''
    # The Dectector keyword retruns 'NRCA*' or 'NRCB*' so to simplify matters
    # I just pull out the 4th character in the string
    if nircam_mod is not None:
        module = nircam_mod
    else:
        module = header['DETECTOR'][3]

    if module == 'A':
        ## NIRCAM Module A
        if fgs_guider == 1:
            ## FGS guider = 2; Perform a Left-Right flip
            im =  np.fliplr(im)# equivalent to im[:,::-1]
        else:
            ## FGS guider = 2; Perform a 180 degree rotation
            im = np.rot90(im,k=20)

    elif module == 'B':
        ## NIRCAM Module B
        if fgs_guider == 1:
            ## FGS guider = 2; Perform a Up-Down flip
            im =  np.flipud(im)# equivalent to im[::-1,...]
        else:
            ## FGS guider = 2; No change necessary!
            pass
    else:
        print('Check the header keyword "DETECTOR" for the NIRCAM module, \
              then re-run using the "nircam_mod" keyword to bypass the header query.')
    return im


def pad_data(data, padding):
    """
    Pad data with median of data with Poisson noise
    """
    size = np.shape(data)[0]

    # Create an array of size binned data + 2xpadding
    background = np.zeros((size+2*(padding), size+2*(padding)))

    # Have a different median value based on pedestal
    # For each specified region in 'background', fill with median value of region in data image
    ped_size = size/4
    medians = []
    for i in range(4):
        medians.append(np.median(data[:,i*ped_size:(i+1)*ped_size-1]))
        background[:,padding+(i*ped_size):padding+((i+1)*ped_size)] = medians[i]
    background[:,:padding]=medians[0]  # Add median value of last pedestal to first padded section
    background[:,-padding:]=medians[3] # Add median value of last pedestal to last padded section
    # Add Poisson noise
    padded_data = np.random.poisson(background)
    # Add in real data
    padded_data[padding:padding+size,padding:padding+size] = data

    return padded_data

def resize_nircam_image(data, NIRCam_scale,FGS_pix,FGS_plate_size):
    cropped = data[4:-4,4:-4] # crop 4pixel zero-padding
    binned_pix = int(round((data.shape[0]*NIRCam_scale*FGS_pix)/(FGS_plate_size*60)))
    data_resized = utils.resize_array(cropped,binned_pix,binned_pix)

    padding = (cropped.shape[0] - binned_pix)/2
    data_pad = pad_data(data_resized, padding)
    fgs_data = np.pad(data_pad, 4, 'constant')

    return fgs_data


def normalize_data(data, fgs_counts, objects=None):
    '''
    Passing in the objects array will ensure that the total counts in the image
    are only distributed to the psfs
    '''
    if objects is not None:
        objs = objects!=0
        psfs_only = data * objs
        data_norm = (float(fgs_counts)/psfs_only.sum()) * psfs_only.astype(np.float64)

    else:
        data_norm = (float(fgs_counts)/data.sum()) * data.astype(np.float64)

    return data_norm



def add_bias_to_data(bias_data_path, FGS_data, root, guider='', output_path='',
                     save_to_fits=True):
    """
    OUT OF DATE - 6/21/17
    Adds in the bias from the guider (two seperate files) to FGS
    (or simulated FGS image that has been padded and normalized)

    This assumes that the format of the guider bias filename is
       job<ID>_g<guider_number>bias.fits.
    If you have a different filename format, pass in the guider_name.

    Default guider bias files are found in the following location:
       guider1: "g1bias.fits"
       guider2: "g2bias.fits"
    """
    header, guider = utils.read_fits(bias_data_path,index=0)

    binned_pad_norm_bias = FGS_data + guider

    if save_to_fits:
        if not os.path.exists(os.path.join(output_path,'bin_norm_bias_imgs')):
            os.makedirs(os.path.join(output_path,'bin_norm_bias_imgs'))

        if guider_name is None:
            guider_name = bias_data_path.split('/')[-1].split('.')[0][-6:]
        out_path =  os.path.join(output_path,'bin_norm_bias_imgs',
                      '{}_G{}_binned_pad_norm.fits'.format(root,guider))
        utils.write_fits(out_path,binned_pad_norm_bias)

    return binned_pad_norm_bias

def find_objects(smoothed_data):
    '''
    Smooth image and pull out individual psfs

    Returns 'objects' - image array with all pixels belonging to one object
    given the same value - and 'num_objects' - the total number of objects
    '''
    # Use Otsu thresholding to pull out only bright regions
    global_otsu = utils.threshold_im(smoothed_data)
    # Use labels to pull out the 18 psfs
    objects, num_objects = ndimage.measurements.label(global_otsu)

    return objects, num_objects


def find_centroids(data, objects, num_objects, root, guider, output_path='',
                   write_out_file=False):
    """
    Given an image with multiple PSFs, smooth first with a Gaussian filter before
    using Otsu thresholding to pull out of the PSFs. Then, label each of these
    PSFs so that their centroids can be found and these coordinates returned.
    """
    ## Pull out coords of each psf
    ## **Coords are Y,X**
    coords = []
    for i in range(1,num_objects+1):
        im = np.copy(objects)
        im[objects!=i]=0
        coords.append(ndimage.measurements.center_of_mass(im*data))

    if write_out_file:
        filename= os.path.join(output_path,'{}_G{}_psf_centers.csv'.format(root,guider))
        utils.write_to_file(filename,coords,labels=['y','x'])

    return coords

def plot_centroids(data,coords,root,guider,output_path):
    plt.clf()
    plt.imshow(data,cmap='Greys_r',norm=LogNorm())
    for i in range(len(coords)):
        plt.scatter(coords[i][1],coords[i][0],color='blue')
    plt.title('Centroids found for {}'.format(root))
    plt.xlim(0,np.shape(data)[0])
    plt.ylim(0,np.shape(data)[0])
    plt.savefig(os.path.join(output_path,'{}_G{}_centers.png'.format(root,guider)))
    plt.close()


def count_rate_total(data, smoothed_data, gs_points, threshold = None,
                     counts_3x3=True):
    """
    Get the x,y, and counts for each psf in the image

    The threshold default of 150, assumes that your data type is uint16. If not,
    a good start is 40.
    """
    if threshold is None:
        threshold = smoothed_data.max() * 0.05
    objects, num_objects = utils.isolate_psfs(smoothed_data,threshold)

    counts = []
    coords = []
    val = []
    for i in range(1,num_objects+1):
        im = np.copy(objects)
        im[objects!=i]=False
        im[objects==i]=True
        coord = ndimage.measurements.center_of_mass(im*data)
        coords.append(coord)
        if counts_3x3:
            counts.append(utils.countrate_3x3(coord,data))
        else:
            counts.append(np.sum(im*data))
        val.append(np.sum(im*1.))

    coords_master = []
    #Pad the range over which the coordinates can match over a sufficiently large
    # range, but not large enough that another PSF center could be within this area
    pad = 12
    #Make a master coordinate list that matches with the coords list, but gives
    # the gs_points coordinates
    for i in range(len(coords)):
        for j in range(len(gs_points)):
            if (coords[i][1]-pad < gs_points[j][1] < coords[i][1]+pad) and (coords[i][0]-pad < gs_points[j][0] < coords[i][0]+pad):
                coords_master.append(gs_points[j])

    y,x = map(list,zip(*coords_master))

    return y, x, counts, val

def create_cols_for_coords_counts(y,x,counts,val):
    """
    Create an array of columns of y, x, and counts of each PSF to be written out.
    Put the PSF with the most compact PSF first in the list to make it the
    Guide Star.
    """

    cols = []
    for i, (yy,xx,co) in enumerate(zip(y,x,counts)):
        cols.append([yy,xx,co]) ## these coordinates are y,x

    #Find most compact PSF
    min_ind = np.where(val == np.min(val))[0][0]
    #Move most compact PSF to top of the list
    cols.insert(0, cols.pop(min_ind))

    return cols

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
def convert_im(input_im, guider, fgs_counts=None, jmag=None, nircam_mod=None):
    '''
    Takes NIRCam image and turns it into an FGS-like image, gets count rate and location of
    each guide star in each image

    Parameters
    ==========
    input_im: str,list of strings
        This can be the path to a file, the path to a directory where all .fits
        files are images you want to convert, or a list of paths (i.e. use glob to
        get this list)
        input_path = '/user/kbrooks/itar/FGS/ga/subset_060617'
    guider: int
        1 or 2
    fgs_counts: int
        The FGS counts in the star. If set to 'None', pass in the jmag. SET TO 65353 FOR NOW
    jmag: float
        The J magnitude of the star. If set to 'None' and fgs_counts set to 'None',
        will defaul to 11.
    nircam_mod: str
        The NIRCAM module, otherwise the header will be parsed
    '''

    ## Establish paths and necessary files
    # 'local_path' is the path where this script exists
    local_path = os.path.dirname(os.path.realpath(__file__))
    # Data path is the directory that includes *.fits files (ie newmagicHdrImg,bias0, etc)
    data_path = os.path.join(local_path,'data')
    # Guider-dependent files
    header_file = os.path.join(data_path,'newG{}magicHdrImg.fits'.format(guider))
    bias_data_path = os.path.join(data_path,'g{}bias0.fits'.format(guider))

    # ---------------------------------------------------------------------
    ## Constants
    NIRCam_scale = 0.032 #NIRCam pixel scale
    FGS_pix = 2048 #FGS image size in pixels
    FGS_plate_size = 2.4 #FGS image size in arcseconds
    ## Constants to change
    BP_thresh = 2000 #Bad pixel threshold

    # ---------------------------------------------------------------------
    ## Find FGS counts to be used for normalization
    if fgs_counts is None:
        if jmag is None:
            print('No counts or J magnitude given, setting to default')
            jmag=11
        fgs_counts = counts_to_jmag.jmag_to_fgs_counts(jmag,guider)
    else:
        jmag = counts_to_jmag.fgs_counts_to_jmag(fgs_counts,guider)

    print('J magnitude = {:.1f}'.format(jmag))

    # ---------------------------------------------------------------------
    ## Get list of images from input path: can take file,list,dir
    if isinstance(input_im,list):
        im_list = input_im
    elif os.path.isfile(input_im):
        im_list = [input_im]
    elif os.path.isdir(input_im):
        im_list = (glob(os.path.join(input_im, '*.fits')))
    else:
        print("Input format not recognized. Exiting.")
        return

    # ---------------------------------------------------------------------
    ## For the images requested, convert to FGS images
    for im in im_list:
        basename = os.path.basename(im)
        root = basename.split('.')[0]
        print('Beginning to create FGS image from {}'.format(root))

        output_path = os.path.join(local_path,'out',root)
        if not os.path.exists(output_path):
            os.makedirs(output_path)

        header, data = utils.read_fits(im)

        # ---------------------------------------------------------------------
        ## Create FGS image
        # Mask out bad pixels
        data_masked = bad_pixel_correction(data,BP_thresh)
        # Rotate the NIRCAM image into FGS frame
        data_rot = rotate_nircam_image(data_masked,guider,header,nircam_mod)
        # Pad image
        data_pad = resize_nircam_image(data, NIRCam_scale,FGS_pix,FGS_plate_size)
        # Find individual psfs
        smoothed_data = ndimage.gaussian_filter(data_pad,sigma=25)
        objects, num_objects = find_objects(smoothed_data)
        # Normalize image
        data_norm = normalize_data(data_pad, fgs_counts)

        out_path = os.path.join(output_path,
                                'FGS_imgs','{}_G{}_binned_pad_norm.fits'.format(root,guider))
        data_norm[data_norm >= 65535] = 65535 #any value about 65535 will wrap when converted to uint16
        hdr = utils.read_fits(header_file,0)[0]
        utils.write_fits(out_path,np.uint16(data_norm),header=hdr)

        # ---------------------------------------------------------------------
        ## Get coordinates and count rate
        coords = find_centroids(data_norm, objects, num_objects, root, guider,
                                output_path=output_path)
        plot_centroids(data_norm,coords,root,guider,output_path)
        y, x, counts, val=count_rate_total(data_norm, smoothed_data,coords,
                                             counts_3x3=True)
        cols=create_cols_for_coords_counts(y,x,counts,val)
        utils.write_cols_to_file(output_path,
                           filename='{0}_G{1}_psf_count_rates.txt'.format(root,guider),
                           labels=['y','x','count rate'],
                           cols=cols)

        # ---------------------------------------------------------------------

        print ("Finished for {}, Guider = {}".format(root,guider))
