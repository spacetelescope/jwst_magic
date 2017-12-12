""" Common utilities for the FGS commissioning tool"""
# STDLIB
import csv
import itertools
import os
import sys
import time

# Third Party
from astropy.io import fits
import numpy as np


def ensure_dir_exists(fullpath):
    """Creates dirs from ``fullpath`` if they do not already exist.
    """
    if not os.path.exists(fullpath):
        os.makedirs(fullpath)


def get_logname(logdir, taskname):
    """Generate log filename based on time stamp and task name.

    Parameters
    ----------
    logdir : str
        Path where log file is stored.

    taskname : str
        Name of the task associated with log file.

    Returns
    -------
    logname : str
        Log filename.

    """
    timestamp = time.strftime('%Y_%m_%d_%a_%H%M%S')
    logname = '{0}_{1}.log'.format(timestamp, taskname)
    return os.path.join(logdir, logname)

def write_fits(outfile, data, header=None):
    '''
    Write data to a simple fits. Assumes one extension and no header.
    '''
    out_dir = os.path.dirname(outfile)
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    hdul = fits.PrimaryHDU(data=data)
    if header is not None:
        hdul.header = header

    hdul.writeto(outfile, overwrite=True)
    print("Successfully wrote: {}".format(outfile))


def write_to_file(filename, rows, labels='', mode='w', fmt='%.4f'):
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

    mode = mode
    if filename.endswith('csv'):
        with open(filename, mode) as fout:
            csvwriter = csv.writer(fout)
            if labels is not None:
                # Write out labels for each column
                csvwriter.writerow(labels)
            for i in rows:
                csvwriter.writerow(i)
    else:
        try:
            np.savetxt(filename, rows, fmt=fmt, header=' '.join(labels))
        except TypeError:
            f = open(filename, 'w')
            f.write('# ' + ' '.join(labels) + '\n')
            for row in rows:
                f.write(' '.join(row) + '\n')


def write_cols_to_file(output_path, filename, labels, cols):
    '''
    Write columns of data to a file
    '''
    filename = os.path.join(output_path, filename)
    write_to_file(filename, cols, labels=labels)


def swap_if_little_endian(data):
    '''
    Swap byte order if little endian
    '''
    if sys.byteorder == 'little':
        data = data.byteswap()
    return data

def swap_if_big_endian(data):
    '''
    Swap byte order if big endian
    '''
    if sys.byteorder == 'big':
        data = data.byteswap()
    return data


# http://stackoverflow.com/questions/8090229/resize-with-averaging-or-rebin-a-numpy-2d-array
def resize_array(arr, new_rows, new_cols):
    '''
    This function takes an 2D numpy array and produces a smaller array
    of size new_rows, new_cols. new_rows and new_cols must be less than
    or equal to the number of rows and columns in a. new_rows and new_columns
    do not have to be integer factors of the original array rows and columns.
    '''
    rows, cols = np.shape(arr)
    yscale = float(rows) / new_rows
    xscale = float(cols) / new_cols

    # first average across the cols to shorten rows
    new_a = np.zeros((rows, new_cols))
    for j in range(new_cols):
        firstx, lastx = j*xscale, (j+1)*xscale
        scale_line = rescale_array(firstx, lastx)
        new_a[:, j] = np.dot(arr[:, int(firstx):int(lastx)+1], scale_line)/scale_line.sum()

    new_arr = np.zeros((new_rows, new_cols))
    for i in range(new_rows):
        firsty, lasty = i*yscale, (i+1)*yscale
        scale_line = rescale_array(firsty, lasty)
        new_arr[i:,] = np.dot(scale_line, new_a[int(firsty):int(lasty)+1,])/scale_line.sum()

    return new_arr
def rescale_array(first, last):
    '''
    Rows can be rows or columns. To be used with resize_array.
    '''
    scale_line = np.ones((int(last)-int(first)+1))
    scale_line[0] = 1 - (first-int(first))
    scale_line[-1] = (last-int(last))
    if last == int(last):
        scale_line = scale_line[:-1]
        last = int(last) - 1
    return scale_line


def find_xy_between_two_points(coords1, coords2):
    '''
    Find the x and y differences between two points
    '''
    diff1 = np.abs((coords1[0] - coords2[0]))
    diff2 = np.abs((coords1[1] - coords2[1]))

    return diff1, diff2


def find_resultant(coords1, coords2):
    '''
    Find the magnitude of the resultant vector described by two points
    '''
    diff1, diff2 = find_xy_between_two_points(coords1, coords2)
    resultant = np.sqrt(diff1**2 + diff2**2)

    return resultant

def find_dist_between_points(coords):
    '''
    Find distances between all combinations of points
    '''
    dists = []
    for c1, c2 in itertools.combinations(coords, 2):
        dists.append(find_resultant(c1, c2))

    return dists
