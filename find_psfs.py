import argparse
from astropy.io import fits
import numpy as np
import csv
import os

from skimage.filters import threshold_otsu
from scipy import ndimage

import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle
from matplotlib.text import Text
from matplotlib.image import AxesImage

#local
from utils import *

def isolate_psfs(smoothed_data,threshold):
    psfs_only = smoothed_data>threshold
    objects, num_objects = ndimage.measurements.label(psfs_only)

    #Check to make sure you have 18 PSFs.
    #The problem with this is that if we have two PSFs really close to each other,
    # the threshold will be very high which will ruin the fun for all the other PSFs.
    # **This is NOT a permanent fix.**

    while num_objects != 18: # assume you want 18 PSFs
        threshold +=5
        psfs_only = smoothed_data > threshold
        objects, num_objects = ndimage.measurements.label(psfs_only)

    return objects, num_objects


def threshold_im(im):
    """
    Uses Otsu thresholding to pull out the smoothed PSFs in the image
    """
    threshold_global_otsu = threshold_otsu(im)
    global_otsu = im >= threshold_global_otsu
    return global_otsu


def find_objects(smoothed_data):
    '''
    Smooth image and pull out individual psfs

    Returns 'objects' - image array with all pixels belonging to one object
    given the same value - and 'num_objects' - the total number of objects
    '''
    # Use Otsu thresholding to pull out only bright regions
    global_otsu = threshold_im(smoothed_data)
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
        write_to_file(filename,coords,labels=['y','x'])

    return coords

def get_psfs_from_plot(data):
    fig, ax = plt.subplots()
    ax.set_title('click on points', picker=True)
    ax.set_ylabel('ylabel', picker=True, bbox=dict(facecolor='red'))
    line, = ax.plot(np.random.rand(100), 'o', picker=5)
    fig.canvas.mpl_connect('pick_event', onpick1)


def onpick1(event):
    coords=[]
    if isinstance(event.artist, Line2D):
        thisline = event.artist
        xdata = thisline.get_xdata()
        ydata = thisline.get_ydata()
        ind = event.ind
        coords.append((np.take(xdata, ind)[0],np.take(ydata, ind)[0]))
        print('X= {}'.format(np.take(xdata, ind)[0])) # Print X point
        print('Y= {}'.format(np.take(ydata, ind)[0])) # Print Y point



def distance_calc((x1,y1),(x2,y2)):
    return np.sqrt( (x2 - x1)**2 + (y2 - y1)**2 )
