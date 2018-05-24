""" Common utilities for the JWST MaGIC package"""
# STDLIB
import csv
import itertools
import os
import sys
import time
import logging.config
import yaml

# Third Party
from astropy.io import fits
import numpy as np

PACKAGE_PATH = os.path.dirname(os.path.realpath(__file__))
LOG_CONFIG_FILE = os.path.join(PACKAGE_PATH, 'data', 'logging.yaml')

def create_logger_from_yaml(module_name, path=LOG_CONFIG_FILE, root='',
                            level=logging.DEBUG):
    '''Set up logger using YAML file

    References
    ----------
    https://fangpenlin.com/posts/2012/08/26/good-logging-practice-in-python/
    https://www.blog.pythonlibrary.org/2016/06/09/python-how-to-create-an-exception-logging-decorator/
    https://docs.python.org/2/howto/logging.html
    '''
    # Parse logging level input
    if type(level) != int:
        if level.upper() == 'DEBUG':
            level = logging.DEBUG
        elif level.upper() == 'INFO':
            level = logging.INFO
        elif level.upper() == 'WARNING':
            level = logging.WARNING
        elif level.upper() == 'ERROR':
            level = logging.ERROR
        elif level.upper() == 'CRITICAL':
            level = logging.CRITICAL
        else:
            raise ValueError('Unknown logging level {}'.format(level))

    # Try to open the config file
    if os.path.exists(path):
        with open(path, 'rt') as f:
            config = yaml.safe_load(f.read())

        # Update the stdout handler level
        config['handlers']['console']['level'] = level

        # Add user's log filename to file handler
        task_name = module_name.split('.')[-1].split('_')
        task_name = ''.join(task_name)
        log_label = '_'.join([task_name, root])
        log_path = os.path.join(os.path.dirname(PACKAGE_PATH), 'logs')
        logfile = get_logname(log_path, log_label)
        config['handlers']['file_handler']['filename'] = logfile

        # # Replace filler "my_module" logger with module-specific logger
        # config['loggers'][name] = config['loggers']['my_module']
        # del config['loggers']['my_module']

        # Create logger from modified dictionary
        logging.config.dictConfig(config)

    # If the config file doesn't exist, just make a super basic logger
    else:
        raise FileNotFoundError('No log config file at {}'.format(path))
        logging.basicConfig(module_name, level=level)

    logger = logging.getLogger(module_name)
    logger.info('Started logging to file {}'.format(logfile))

    return logger

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
    print("Successfully wrote: {}".format(filename))


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

    # First average across the cols to shorten rows
    new_a = np.zeros((rows, new_cols))
    for j in range(new_cols):
        # Calculate the (fractional) starting and ending columns that will be
        # averaged into one column (e.g. averaging from column 4.3 to 7.8)
        firstx, lastx = j * xscale, (j + 1) * xscale
        # Calculate list of scaling factors of each column that is being averaged
        # e.g. if avging from 4.3 to 7.8: [0.7, 1, 1, 0.8]
        #                        column:    4   5  6   7
        scale_line = rescale_array(firstx, lastx)

        # Fill new array with averaged columns
        try:
            new_a[:, j] = np.dot(arr[:, int(firstx):int(lastx) + 1], scale_line) / scale_line.sum()
        except ValueError:
            # If needed, crop the scaling list to match the number of columns
            scale_line = scale_line[:-1]
            new_a[:, j] = np.dot(arr[:, int(firstx):int(lastx) + 1], scale_line) / scale_line.sum()

    # Then average across the rows to produce the final array
    new_arr = np.zeros((new_rows, new_cols))
    for i in range(new_rows):
        # Calculate the (fractional) starting and ending rows that will be
        # averaged into one row
        firsty, lasty = i * yscale, (i + 1) * yscale
        # Calculate scaling factors of each row that is being averaged
        scale_line = rescale_array(firsty, lasty)

        # Fill new array with averaged rows
        try:
            new_arr[i:,] = np.dot(scale_line, new_a[int(firsty):int(lasty) + 1,]) / scale_line.sum()
        except ValueError:
            # If needed, crop the scaling list to match the number of rows
            scale_line = scale_line[:-1]
            new_arr[i:,] = np.dot(scale_line, new_a[int(firsty):int(lasty) + 1,]) / scale_line.sum()

    return new_arr

def rescale_array(first, last):
    '''
    Rows can be rows or columns. To be used with resize_array.
    '''
    scale_line = np.ones((int(last) - int(first) + 1))
    scale_line[0] = 1 - (first - int(first))
    scale_line[-1] = (last - int(last))
    if last == int(last):
        # scale_line = scale_line[:-1]
        scale_line[-1] = 0.  # Changed 1/16/18 to fix bug with truncating scale_line
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


def correct_image(image, upper_threshold=65000, upper_limit=65000):
    '''
    Correct image for negative and saturated pixels
    '''
    img = np.copy(image)
    img[img < 0] = 0.            # neg pixs -> 0
    img[img >= upper_threshold] = upper_limit    # sat'd pixs -> 65K
    img[np.isfinite(img) == 0] = 0.

    return img

def countrate_3x3(x, y, data):
    """
    Using the coordinates of each PSF, place a 3x3 box around center pixel and sum
    the counts of the pixels in this box.
    """
    x = int(x)
    y = int(y)

    counts = np.sum(data[y - 1:y + 2, x - 1:x + 2])
    return counts

def get_guider(header, guider=None, log=None):
    '''Read the guider information from the FITS header'''
    try:
        hdr_guider = int(header['DETECTOR'][-1])
    except KeyError:
        return

    # Handle NIRCam images being passed if the guider hasn't been commanded
    if not guider and not header['DETECTOR'].startswith('GUIDER'):
        try:
            log.warning("The header indicates that this is a NIRCam image; the " +
                        "guider number cannot be parsed from the header. Setting " +
                        " to default value of guider = 1.")
        except AttributeError:
            print("The header indicates that this is a NIRCam image; the " +
                  "guider number cannot be parsed from the header. Setting " +
                   " to default value of guider = 1.")
        hdr_guider = 1

    # Can compare with human-commanded guider number
    if guider:
        # Override just in case the human gets it wrong for an FGS image
        if hdr_guider != guider and header['DETECTOR'].startswith('GUIDER'):
            try:
                log.warning("Image Conversion: The header indicates that this is a " +
                            "guider {0} image. Processing as a guider {0} image.".format(hdr_guider))
            except:
                print("Image Conversion: The header indicates that this is a " +
                      "guider {0} image. Processing as a guider {0} image.".format(hdr_guider))
        # Otherwise, if it is a NIRCam image, just take what the user commands
        else:
            hdr_guider = guider

    # Make sure it's a real guider!
    if hdr_guider not in [1, 2]:
        try:
            log.warning("Invalid guider number provided/found: {}".format(hdr_guider))
        except AttributeError:
            print("Invalid guider number provided/found: {}".format(hdr_guider))

    return hdr_guider

def get_data_and_header(filename):
    '''Open a FITS file, get the data from either the SCI or PRIMARY
    extension, and get the header from the PRIMARY extension.
    '''
    with fits.open(filename) as hdulist:
        if len(hdulist) > 1:
            data = hdulist[1].data
        else:
            data = hdulist[0].data
        header = hdulist[0].header

    return data, header

def make_root(root, filename):
    '''If no root has been provided, extract from base filename'''
    if root is None:
        if filename[-11:] == 'ALLpsfs.txt':
            root = os.path.basename(filename).split('.')[0][:-11]
        else:
            root = os.path.basename(filename).split('.')[0]

    return root

def make_out_dir(out_dir, default_out_path, root):
    # Determine output directory
    if out_dir is None:
        out_dir = os.path.join(default_out_path, 'out', root)
    else:
        out_dir = os.path.join(out_dir, 'out', root)

    return out_dir
