import argparse
from astropy.io import fits
import numpy as np
import csv
import os

from skimage.filters import threshold_otsu
from scipy import ndimage,signal

import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle
from matplotlib.text import Text
from matplotlib.image import AxesImage
import matplotlib
matplotlib.rcParams['image.origin'] = 'upper'

#local
from utils import *



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
        coords.append(ndimage.measurements.center_of_mass(im*data)) #y,x

    if write_out_file:
        filename= os.path.join(output_path,'{}_G{}_psf_centers.csv'.format(root,guider))
        write_to_file(filename,coords,labels=['y','x'])

    return coords


def distance_calc((x1,y1),(x2,y2)):
    return np.sqrt( (x2 - x1)**2 + (y2 - y1)**2 )


def plot_centroids(data,coords,root,guider,output_path,compact=False):
    #labels = ['gs','rs1','rs2','rs3','rs4','rs5','rs6','rs7','rs8','rs9','rs10',
    #          'rs11','rs12','rs13','rs14','rs15','rs16','rs17']
    if compact:
        pad = 150
    else:
        pad = 500

    plt.figure(figsize=(17,17))
    plt.imshow(data,cmap='Greys',norm=LogNorm())
    for i in range(len(coords)):
        plt.scatter(coords[i][1],coords[i][0])
        plt.annotate('({},{})'.format(int(coords[i][0]),int(coords[i][1])),(coords[i][1]-(pad*0.05),coords[i][0]+(pad*0.05)))
    plt.title('Centroids found for {}'.format(root))
    plt.xlim(np.shape(data)[1]/2-pad,np.shape(data)[1]/2+pad)
    plt.ylim(np.shape(data)[0]/2+pad,np.shape(data)[0]/2-pad)
    plt.savefig(os.path.join(output_path,'{}_G{}_centers.png'.format(root,guider)))
    plt.close()


def count_rate_total(data, smoothed_data, gs_points, dist=None, threshold = None,
                     counts_3x3=True,num_psfs=18):
    """
    Get the x,y, and counts for each psf in the image

    The threshold default of 150, assumes that your data type is uint16. If not,
    a good start is 40.
    """
    if threshold is None:
       threshold = smoothed_data.max() * 0.05
    objects, num_objects = isolate_psfs(smoothed_data,threshold,num_psfs)

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
            counts.append(countrate_3x3(coord,data))
        else:
            counts.append(np.sum(im*data))
        val.append(np.sum(im*1.))

    coords_master = []
    #Pad the range over which the coordinates can match over a sufficiently large
    # range, but not large enough that another PSF center could be within this area
    if dist is None:
        pad = 12
    else:
        pad = dist-1
    #Make a master coordinate list that matches with the coords list, but gives
    # the gs_points coordinates
    for i in range(len(coords)):
        for j in range(len(gs_points)):
            if (coords[i][1]-pad < gs_points[j][1] < coords[i][1]+pad) and (coords[i][0]-pad < gs_points[j][0] < coords[i][0]+pad):
                coords_master.append(gs_points[j])

    y,x = map(list,zip(*coords_master))
    return y, x, counts, val

def create_cols_for_coords_counts(y,x,counts,val,inds=None):
    """
    Create an array of columns of y, x, and counts of each PSF to be written out.
    Use the inds returned from pick_stars based on user input.

    If no inds are given, put the PSF with the most compact PSF first in the list to make it the
    Guide Star. **This method is not fool-proof, use at own risk***
    """

    cols = []
    for i, (yy,xx,co) in enumerate(zip(y,x,counts)):
        cols.append([yy,xx,co]) ## these coordinates are y,x

    if inds is None:
        min_ind = np.where(val == np.min(val))[0][0] #Find most compact PSF
        cols.insert(0, cols.pop(min_ind)) #Move most compact PSF to top of the list
    else:
        cols = [cols[i] for i in inds]

    return cols


class SelectStars(object):
    def __init__(self,fig,ax,data,xarray,yarray,dist=15):
        self.data = data
        self.xt = xarray
        self.yt = yarray

        self._ind = None
        self.coords=[]

        self.epsilon = dist

        self.fig = fig
        self.ax = ax
        self.canvas = fig.canvas
        self.inds=[]

        self.cid = []
        self.cid.extend((
            self.canvas.mpl_connect('button_press_event', self.button_press_callback),
            ))

        self.canvas.mpl_disconnect(self.canvas.manager.key_press_handler_id)
        plt.show()

        return


    def button_press_callback(self, event):
        'whenever a mouse button is pressed'
        if event.inaxes==None: return
        if event.button != 1: return
        self._ind = self.get_ind_under_point(event)


    def get_ind_under_point(self, event):
        'get the index of the vertex under point if within epsilon tolerance'
        # xt and yt are arrays of the points to compare with
        d = np.sqrt((self.xt-event.xdata)**2 + (self.yt-event.ydata)**2)
        i = np.unravel_index(np.nanargmin(d), d.shape)
        indseq = np.nonzero(np.equal(d, np.amin(d)))[0]
        ind = indseq[0]

        if d[ind]>=self.epsilon:
            print('No star within {} pixels. No star selected.'.format(self.epsilon))
            return
        elif ind in self.inds:
            print('Star already selected, please choose another star')
        else:
            print('Star selected: x={:.1f}, y={:.1f}'.format(self.xt[ind],self.yt[ind]))
            self.inds.append(ind)

        return ind


def pick_stars(data,xarray,yarray,dist,root='',compact=False):
        '''
        Parameters:
            data: str, ndarray
                Either the full path to the data or an array
            xarray: array, list
                Array of x coordinates of each star in the image
            yarray: array, list
                Array of y coordinates of each star in the image
            root: str
                Root name of the images
        '''
        if isinstance(data, str):
            f = fits.open(data)
            data = f[0].data

        if compact:
            pad = 150
        else:
            pad = 300

        print("Click as near to the center of the star as possible.\n \
               The first star that is choosen will be the guide star.\n \
               All additional stars that are clicked on, are the\n \
               reference stars.")

        fig = plt.figure(figsize=(8,8))
        ax = fig.add_subplot(111)
        ax.imshow(data,cmap='coolwarm',interpolation='nearest',norm=LogNorm())
        ax.set_ylim(np.shape(data)[0]/2+pad,np.shape(data)[0]/2-pad)
        ax.set_xlim(np.shape(data)[1]/2-pad,np.shape(data)[1]/2+pad)

        fig.show()

        obj = SelectStars(fig,ax,data,xarray,yarray,dist=dist)

        return obj.inds


def create_reg_file(data, root, guider, output_path, return_nref=False,
                    num_psfs=18, compact=False):
    if isinstance(data,str):
        data = read_fits(filename)[1]

    if compact:
        smoothed_data = signal.medfilt(data,5)
        objects, num_objects = isolate_psfs(smoothed_data,threshold=None,num_psfs=num_psfs)
    else:
        smoothed_data = ndimage.gaussian_filter(data,sigma=25)
        objects, num_objects = find_objects(smoothed_data)


    coords = find_centroids(data, objects, num_objects, root, guider,
                            output_path=output_path)

    dist = np.floor(np.min(find_dist_between_points(coords))) - 1. #find the minimum distance between PSFs

    plot_centroids(data,coords,root,guider,output_path,compact=compact)
    y, x, counts, val = count_rate_total(data, smoothed_data,coords,dist=dist,
                                         counts_3x3=True,num_psfs=num_psfs)
    inds = pick_stars(data,x,y,dist,root=root)

    cols = create_cols_for_coords_counts(y,x,counts,val,inds=inds)
    write_cols_to_file(output_path,
                       filename='{0}_G{1}_regfile.txt'.format(root,guider),
                       labels=['y','x','count rate'],
                       cols=cols)
    if return_nref:
        nref = len(inds)-1
        return nref
