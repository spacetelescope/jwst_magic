#!/usr/bin/python

from astropy.io import fits
import numpy as np
import csv
import os
from skimage.filters import threshold_otsu
from scipy import ndimage
import itertools

import matplotlib.pyplot as plt

def write_fits(outfile,data,header=None):
    '''
    Write data to a simple fits. Assumes one extension and no header.
    '''
    out_dir = os.path.dirname(outfile)
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    hdul = fits.PrimaryHDU(data=data)
    if header is not None:
        hdul.header = header

    hdul.writeto(outfile,overwrite=True)
    print("Successfully wrote: {}".format(outfile))


def read_fits(filename):
    '''
    Get the header and data (for the specified index) from a fits file
    '''
    f = fits.open(filename)
    try:
        if f[1].name == 'SCI':
            data = f[1].data #Try first for the SCI extension
    except IndexError:
        data = f[0].data #If SCI extension doesn't exist, use PrimaryHDU

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
    out_dir = os.path.dirname(filename)
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

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
        np.savetxt(filename,rows,fmt='%.4f',header=' '.join(labels))#,fmt='%.4f'


def write_cols_to_file(output_path, filename, labels, cols):
    filename= os.path.join(output_path,filename)
    write_to_file(filename,cols,labels=labels)


def swap_if_little_endian(data):
    if sys.byteorder == 'little':
        data = data.byteswap()
    return data

def swap_if_big_endian(data):
    if sys.byteorder == 'big':
        data = data.byteswap()
    return data

def convert_fits_to_dat(infile,obsmode,out_dir,root=None):
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
        header,data = read_fits(infile)

        filename = infile.split('/')[-1]
        root = filename.split('.')[0]
    else:
        root = '{}_{}'.format(root,obsmode)

    print("Converting {}.fits to .dat".format(root))

    outfile = '{}.dat'.format(root)
    #data = swap_if_little_endian(data)
    fl = data.flatten()

    if (obsmode == 'PSF') or (obsmode == 'TRK'):
        # ascii float format
        f = '{:16.7e} '

    elif (obsmode == 'ID') or (obsmode == 'ACQ1') or (obsmode == 'ACQ2') or (obsmode == 'ACQ') or (obsmode == 'CAL'):
        # ascii hex dat format
        f = '{:04X} '

    else:
        print("Observation mode not recognized. Returning.")

    with open(os.path.join(out_dir,outfile), 'w') as file_out:
        for i,d in enumerate(fl.astype(np.uint16)):
            file_out.write(f.format(d))

    print("Successfully wrote: {}".format(os.path.join(out_dir,outfile)))
    return


# http://stackoverflow.com/questions/8090229/resize-with-averaging-or-rebin-a-numpy-2d-array
def resize_array(a, new_rows, new_cols):
    '''
    This function takes an 2D numpy array and produces a smaller array
    of size new_rows, new_cols. new_rows and new_cols must be less than
    or equal to the number of rows and columns in a. new_rows and new_columns
    do not have to be integer factors of the original array rows and columns.
    '''
    rows = len(a)
    cols = len(a[0])
    yscale = float(rows) / new_rows
    xscale = float(cols) / new_cols

    # first average across the cols to shorten rows
    new_a = np.zeros((rows, new_cols))
    for j in range(new_cols):
        # get the indices of the original array we are going to average across
        the_x_range = (j*xscale, (j+1)*xscale)
        firstx = int(the_x_range[0])
        lastx = int(the_x_range[1])
        # figure out the portion of the first and last index that overlap
        # with the new index, and thus the portion of those cells that
        # we need to include in our average
        x0_scale = 1 - (the_x_range[0]-int(the_x_range[0]))
        xEnd_scale =  (the_x_range[1]-int(the_x_range[1]))
        # scale_line is a 1d array that corresponds to the portion of each old
        # index in the_x_range that should be included in the new average
        scale_line = np.ones((lastx-firstx+1))
        scale_line[0] = x0_scale
        scale_line[-1] = xEnd_scale
        # Make sure you don't screw up and include an index that is too large
        # for the array. This isn't great, as there could be some floating
        # point errors that mess up this comparison.
        if scale_line[-1] == 0:
            scale_line = scale_line[:-1]
            lastx = lastx - 1
        # Now it's linear algebra time. Take the dot product of a slice of
        # the original array and the scale_line
        new_a[:,j] = np.dot(a[:,firstx:lastx+1], scale_line)/scale_line.sum()
    # Then average across the rows to shorten the cols. Same method as above.
    # It is probably possible to simplify this code, as this is more or less
    # the same procedure as the block of code above, but transposed.
    # Here I'm reusing the variable a. Sorry if that's confusing.
    a = np.zeros((new_rows, new_cols))
    for i in range(new_rows):
        the_y_range = (i*yscale, (i+1)*yscale)
        firsty = int(the_y_range[0])
        lasty = int(the_y_range[1])
        y0_scale = 1 - (the_y_range[0]-int(the_y_range[0]))
        yEnd_scale =  (the_y_range[1]-int(the_y_range[1]))
        scale_line = np.ones((lasty-firsty+1))
        scale_line[0] = y0_scale
        scale_line[-1] = yEnd_scale
        if scale_line[-1] == 0:
            scale_line = scale_line[:-1]
            lasty = lasty - 1
        a[i:,] = np.dot(scale_line, new_a[firsty:lasty+1,])/scale_line.sum()

    return a

def find_xy_between_two_points(coords1,coords2):
    diff1 = np.abs((coords1[0] - coords2[0]))
    diff2 = np.abs((coords1[1] - coords2[1]))

    return diff1, diff2


def find_resultant(coords1, coords2):
    diff1, diff2 = find_xy_between_two_points(coords1,coords2)

    z = np.sqrt(diff1**2 + diff2**2)

    return z

def find_dist_between_points(coords):
    dists = []
    for c1, c2 in itertools.combinations(coords, 2):
        dists.append(find_resultant(c1,c2))

    return dists


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


def isolate_psfs(smoothed_data,threshold,num_psfs=18):
    if threshold is None:
        threshold = smoothed_data.max() * 0.05
    psfs_only = smoothed_data>threshold
    objects, num_objects = ndimage.measurements.label(psfs_only)

    #Check to make sure you have num_psfs PSFs.
    #The problem with this is that if we have two PSFs really close to each other,
    # the threshold will be very high which will ruin the fun for all the other PSFs.
    # **This is NOT a permanent fix.**

    while num_objects > num_psfs: # assume you want 18 PSFs
        #print("Num of objects detected: {}".format(num_objects))
        threshold = .95*threshold
        psfs_only = smoothed_data > threshold
        objects, num_objects = ndimage.measurements.label(psfs_only)

    while num_objects < num_psfs: # assume you want 18 PSFs
        #print("Num of objects detected: {}".format(num_objects))
        threshold = 1.05*threshold
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
def countrate_3x3(coords,data):
    """
    Using the coordinates of each PSF,place a 3x3 box around center pixel, and sum
    the counts of the pixels in this box.
    """
    xx = int(coords[1])
    yy = int(coords[0])

    counts = np.sum(data[yy-1:yy+2,xx-1:xx+2])
    return counts

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
