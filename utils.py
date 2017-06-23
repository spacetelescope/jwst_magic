#!/usr/bin/python

from astropy.io import fits
import numpy as np
import csv
import os
from skimage.filters import threshold_otsu
from scipy import ndimage

def write_to_fits(outfile,data,header=None):
    '''
    Write data to a simple fits. Assumes one extension and no header.
    '''
    hdul = fits.PrimaryHDU(data=data)
    if header is not None:
        hdul.header = header

    hdul.writeto(outfile,clobber=True)
    print("Successfully wrote: {}".format(outfile))


def get_fits_data(filename, index = 1):
    '''
    Get the header and data (for the specified index) from a fits file
    '''
    f = fits.open(filename)
    data = f[index].data
    header = f[0].header
    f.close()

    return header, data


def write_to_file(filename,rows,labels=None,mode='w'):
    """ Write out results to a csv, dat, txt, etc file.

    Parameters
    -----------
    filename : str
        Filename of output file

    rows : list
        List of results, one entry per row in file

    labels : list, optional
        List of labels for each column
    """
    mode=mode
    if filename.endswith('csv'):
        with open(filename, mode) as fout:
            csvwriter = csv.writer(fout)
            if labels is not None:
                # Write out labels for each column
                csvwriter.writerow(labels)
            for i in rows:
                csvwriter.writerow(i)
    else:
        np.savetxt(filename,rows,fmt='%.4f',header=' '.join(labels))


def write_cols_to_file(output_path, filename, labels, cols):
    filename= os.path.join(output_path,filename)
    write_to_file(filename,cols,labels=labels)



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

## All the ways to get the countrate:
def count_rate_3x3(root, data, gs_points, radius):
    """
    Using the coordinates of each PSF,place a 3x3 box around center pixel, and sum
    the counts of the pixels in this box.
    """

    for point in gs_points:

        xx = int(point[1])
        yy = int(point[0])

        tmp2 = []
        tmp = []
        for ii in range(-1,2):
            for jj in range(-1,2):
                tmp.append(p_data[yy+jj,xx+ii])

        tmp2.append(sum(tmp))
    return tmp2

def get_countrate(x,y,arr):
    """
    If the countrate isn't given by a file, use the coords given by the
    reg file to find the countrates for each PSF.
    """
    countrate=[] #create countrate array of len nstars
    for i in range(len(x)):
        x1 = x[i]-50
        x2 = x[i]+50
        y1 = y[i]-50
        y2 = y[i]+50
        if x1 < 50:
            x1=0
        if x2 > 2047:
            x2=2047
        if y1 < 50:
            y1=0
        if y2> 2047:
            y2=2047
        countrate.append(np.sum(arr[y1:y2,x1:x2])) # counts/sec NOT ADU
        print('Countrate: {}'.format(np.sum(arr[y1:y2,x1:x2])))

    countrate = np.asarray(float(countrate))
    return countrate
