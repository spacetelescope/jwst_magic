""" Common utilities for the JWST MaGIC package.

Authors
-------
    - Keira Brooks
    - Lauren Chambers
    - Shannon Osborne

Use
---
    ::
        from jwst_magic import utils
        utils.ensure_dir_exists(dir)
"""
# STDLIB
from collections import OrderedDict
import csv
import itertools
import logging
import logging.config
import os
import re
import requests
import string
import socket
import sys
import time
import yaml

# Third Party
from astropy.io import fits
import numpy as np
import pandas as pd

PACKAGE_PATH = os.path.dirname(os.path.realpath(__file__)).split('utils')[0]
LOG_CONFIG_FILE = os.path.join(PACKAGE_PATH, 'data', 'logging.yaml')

# Start logger
LOGGER = logging.getLogger(__name__)


class CustomFormatter(logging.Formatter):
    """Logging Formatter to add colors and count warning / errors"""

    grey = "\x1b[38;21m"
    yellow = "\x1b[33;21m"
    red = "\x1b[31;21m"
    bold_red = "\x1b[31;1m"
    reset = "\x1b[0m"
    format = "%(asctime)s - %(levelname)s - %(message)s"

    FORMATS = {
        logging.DEBUG: grey + format + reset,
        logging.INFO: grey + format + reset,
        logging.WARNING: yellow + format + reset,
        logging.ERROR: red + format + reset,
        logging.CRITICAL: bold_red + format + reset
    }

    def format(self, record):
        log_fmt = self.FORMATS.get(record.levelno)
        formatter = logging.Formatter(log_fmt)
        return formatter.format(record)


def create_logger_from_yaml(module_name, path=LOG_CONFIG_FILE, root='',
                            level=logging.DEBUG):
    """Set up logger using YAML file

    References
    ----------
    https://fangpenlin.com/posts/2012/08/26/good-logging-practice-in-python/
    https://www.blog.pythonlibrary.org/2016/06/09/python-how-to-create-an-exception-logging-decorator/
    https://docs.python.org/2/howto/logging.html
    """
    # Ensure the logs directory exists
    log_path = determine_log_path()

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
        logfile = get_logname(log_path, log_label)
        config['handlers']['file_handler']['filename'] = logfile

        # # Replace filler "my_module" logger with module-specific logger
        # config['loggers'][name] = config['loggers']['my_module']
        # del config['loggers']['my_module']

        # Create logger from modified dictionary
        logging.config.dictConfig(config)

    # If the config file doesn't exist, raise an error.
    else:
        raise FileNotFoundError('No log config file at {}'.format(path))
        # Possible write-around if need be: use a simple logger
        # logging.basicConfig(module_name, level=level)

    logger = logging.getLogger(module_name)
    logger.info('Started logging to file {}'.format(logfile))

    return logger, logfile


def determine_log_path():
    """Determine whether to save log files in a shared log directory on
    SOGS, or in the default ``logs`` directory in the package directory.
    Ensure the chosen log directory exists.
    """
    if on_sogs_network():
        log_path = "***REMOVED***/guiding/MAGIC_logs/"
    else:
        log_path = os.path.join(os.path.dirname(PACKAGE_PATH), 'logs')

    ensure_dir_exists(log_path)

    return log_path


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


def on_sogs_network():
    """Determine whether MAGIC is being run on the SOGS network.
    """
    return "sogs" in socket.gethostname()


def write_fits(outfile, data, header=None, log=None):
    """
    Write data to a simple fits. Assumes one extension and no header.
    """
    out_dir = os.path.dirname(outfile)
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    if not any([isinstance(header, fits.header.Header), header is None]):
        raise TypeError(f'Header to be written out in {outfile} is not either "None" or of type fits.header.Header')

    hdul = fits.PrimaryHDU(data=data, header=header)

    hdul.writeto(outfile, overwrite=True)

    if log is not None:
        log.info(f"Successfully wrote: {outfile}")
    else:
        print(f"Successfully wrote: {outfile}")


def write_multiext_fits(outfile, data_list, header_list, log=None):
    """
    Write data to a multi-extension fits file. Expands on write_fits()
    """
    out_dir = os.path.dirname(outfile)
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    for header in header_list:
        if not any([isinstance(header, fits.header.Header), header is None]):
            raise TypeError('Header to be written out in {} is not either "None" or of type fits.header.Header'.format(
                outfile))

    hdu_list = []
    for i, (data, header) in enumerate(zip(data_list, header_list)):
        if i == 0:
            hdu = fits.PrimaryHDU(data=data, header=header)
        else:
            hdu = fits.ImageHDU(data=data, header=header)
        hdu_list.append(hdu)
    hdul = fits.HDUList(hdu_list)

    hdul.writeto(outfile, overwrite=True)

    if log is not None:
        log.info("Successfully wrote: {}".format(outfile))
    else:
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

    mode : str, optional
        Mode with which to write the file (e.g. 'w', 'a', etc.)

    fmt : str or list of strings, optional
        Format string for the text being written.
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


def write_cols_to_file(file_path, labels, cols, log=None):
    """
    Write columns of data to a file
    """
    write_to_file(file_path, cols, labels=labels)

    if log is not None:
        log.info("Successfully wrote: {}".format(file_path))
    else:
        print("Successfully wrote: {}".format(file_path))


def swap_if_little_endian(data):
    """
    Swap byte order if little endian
    """
    if sys.byteorder == 'little':
        data = data.byteswap()
    return data


def swap_if_big_endian(data):
    """
    Swap byte order if big endian
    """
    if sys.byteorder == 'big':
        data = data.byteswap()
    return data


# http://stackoverflow.com/questions/8090229/resize-with-averaging-or-rebin-a-numpy-2d-array
def resize_array(arr, new_rows, new_cols):
    """
    This function takes an 2D numpy array and produces a smaller array
    of size new_rows, new_cols. new_rows and new_cols must be less than
    or equal to the number of rows and columns in a. new_rows and new_columns
    do not have to be integer factors of the original array rows and columns.
    """
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
            new_arr[i:, ] = np.dot(scale_line, new_a[int(firsty):int(lasty) + 1, ]) / scale_line.sum()
        except ValueError:
            # If needed, crop the scaling list to match the number of rows
            scale_line = scale_line[:-1]
            new_arr[i:, ] = np.dot(scale_line, new_a[int(firsty):int(lasty) + 1, ]) / scale_line.sum()

    return new_arr


def rescale_array(first, last):
    """
    Rows can be rows or columns. To be used with resize_array.
    """
    scale_line = np.ones((int(last) - int(first) + 1))
    scale_line[0] = 1 - (first - int(first))
    scale_line[-1] = (last - int(last))
    if last == int(last):
        # scale_line = scale_line[:-1]
        scale_line[-1] = 0.  # Changed 1/16/18 to fix bug with truncating scale_line
    return scale_line


def find_xy_between_two_points(coords1, coords2):
    """
    Find the x and y differences between two points
    """
    diff1 = np.abs((coords1[0] - coords2[0]))
    diff2 = np.abs((coords1[1] - coords2[1]))

    return diff1, diff2


def find_resultant(coords1, coords2):
    """
    Find the magnitude of the resultant vector described by two points
    """
    diff1, diff2 = find_xy_between_two_points(coords1, coords2)
    resultant = np.sqrt(diff1**2 + diff2**2)

    return resultant


def find_dist_between_points(coords):
    """
    Find distances between all combinations of points
    """
    dists = []
    for c1, c2 in itertools.combinations(coords, 2):
        dists.append(find_resultant(c1, c2))

    return dists


def correct_image(image, upper_threshold=65000, upper_limit=65000):
    """
    Correct image for negative and saturated pixels
    """
    img = np.copy(image)
    img[img < 0] = 0.            # neg pixs -> 0
    img[img >= upper_threshold] = upper_limit    # sat'd pixs -> 65K
    img[np.isfinite(img) == 0] = 0.

    return img


def count_rate_total(data, objects, num_objects, x, y, countrate_3x3=True):
    """Get the count rates within each object from a segmentation image.

    Parameters
    ----------
    data : 2-D numpy array
        Image data
    objects : 2-D numpy array
        Segmentation of image data
    num_objects : int
        Number of individual objects in the segmentation data
    x : list
        List of x-coordinates of identified PSFs
    y : list
        List of y-coordinates of identified PSFs
    countrate_3x3 : bool, optional
        Calculate the value of the 3x3 square (True), or of the entire
        object (False)

    Returns
    -------
    countrate : list
        List of count rates of each segmentation object
    val : list
        List of number of pixels within each segmentation object
    """

    countrate = []
    val = []
    for i in range(1, num_objects + 1):
        im = np.copy(objects)
        im[objects != i] = False
        im[objects == i] = True

        if countrate_3x3:
            countrate.append(get_countrate_3x3(x[i - 1], y[i - 1], np.array(data)))
        else:
            countrate.append(np.sum(im * data))
        val.append(np.sum(im * 1.))  # Number of pixels in object

    return countrate, val


def get_countrate_3x3(x, y, data):
    """
    Using the coordinates of each PSF, place a 3x3 box around center pixel and sum
    the countrate of the pixels in this box.
    """
    data = correct_image(data)
    data = np.uint16(data)
    x = int(x)
    y = int(y)

    countrate = np.sum(data[y - 1:y + 2, x - 1:x + 2])
    return countrate


def get_guider(header, guider=None, log=None):
    """Read the guider information from the FITS header"""
    try:
        hdr_guider = int(header['DETECTOR'][-1])
    except KeyError:
        return

    # Handle NIRCam images being passed if the guider hasn't been commanded
    if not guider and not header['DETECTOR'].startswith('GUIDER'):
        if log is not None:
            log.warning("The header indicates that this is a NIRCam image; the " +
                        "guider number cannot be parsed from the header. Setting " +
                        " to default value of guider = 1.")
        else:
            print("The header indicates that this is a NIRCam image; the " +
                  "guider number cannot be parsed from the header. Setting " +
                  " to default value of guider = 1.")
        hdr_guider = 1

    # Can compare with human-commanded guider number
    if guider:
        # Override just in case the human gets it wrong for an FGS image
        if hdr_guider != guider and header['DETECTOR'].startswith('GUIDER'):
            if log is not None:
                log.warning("Image Conversion: The header indicates that this is a " +
                            "guider {0} image. Processing as a guider {0} image.".format(hdr_guider))
            else:
                print("Image Conversion: The header indicates that this is a " +
                      "guider {0} image. Processing as a guider {0} image.".format(hdr_guider))
        # Otherwise, if it is a NIRCam image, just take what the user commands
        else:
            hdr_guider = guider

    # Make sure it's a real guider!
    if hdr_guider not in [1, 2]:
        if log is not None:
            log.warning("Invalid guider number provided/found: {}".format(hdr_guider))
        else:
            print("Invalid guider number provided/found: {}".format(hdr_guider))

    return hdr_guider


def get_data_and_header(filename):
    """Open a FITS file, get the data from either the SCI or PRIMARY
    extension, and get the header from the PRIMARY extension.
    """
    with fits.open(filename) as hdulist:
        if len(hdulist) > 1:
            data = hdulist[1].data
        else:
            data = hdulist[0].data
        header = hdulist[0].header

    return data, header


def make_root(root, filename):
    """If no root has been provided, extract from base filename"""
    if root is None:
        if filename[-11:] == 'ALLpsfs.txt':
            root = os.path.basename(filename).split('.')[0][:-11]
        else:
            root = os.path.basename(filename).split('.')[0]

    return root


def make_out_dir(out_dir, default_out_path, root):
    # Determine output directory
    if root is None:
        raise TypeError('Please define the root directory.')

    if out_dir is None:
        out_dir = os.path.join(default_out_path, 'out', root)
    else:
        out_dir = os.path.join(out_dir, 'out', root)

    return out_dir


def create_cols_for_coords_counts(x, y, countrate, val=None, labels=None, inds=None):
    """Format position and count rate data to be written to file.

    Create an array of columns of y, x, and countrate of each PSF to be
    written out. Use the inds returned from pick_stars based on user
    input. If no inds are given, put the PSF with the most compact PSF
    first in the list to make it the guide star.


    Parameters
    ----------
    x : list
        List of x-coordinates of identified PSFs
    y : list
        List of y-coordinates of identified PSFs
    countrate : list
        List of count rates of identified PSFs
    val : list, optional
        List of the number of pixels in each PSF's segmentation object.
        Used to re-order PSFs by compact-ness
    labels : list, optional
        Denotes whether the PSF alphabetic labels should be included as
        a column to write out
    inds : list, optional
        List of the indices of the guide and reference stars

    Returns
    -------
    cols : list
        List of segment positions, count rates, and maybe labels for
        each selected segments
    """
    if labels is not None:
        # NOTE: these coordinates are y, x
        cols = [[ll, '{:.4f}'.format(yy),
                 '{:.4f}'.format(xx),
                 '{:.4f}'.format(co)] for ll, yy, xx, co in zip(labels, y, x, countrate)]
    else:
        # NOTE: these coordinates are y, x
        cols = [[yy, xx, co] for yy, xx, co in zip(y, x, countrate)]

    if inds is None and val is not None:
        min_ind = np.where(val == np.min(val))[0][0]  # Find most compact PSF
        cols.insert(0, cols.pop(min_ind))  # Move most compact PSF to top of the list
    else:
        cols = [cols[i] for i in inds]

    return cols


def match_psfs_to_segments(x, y, smoothing):
    """Match PSFs found in the image to their alphabetic label (between A and R)

    Parameters
    ----------
    x : list
        List of x-coordinates of identified PSFs
    y : list
        List of y-coordinates of identified PSFs
    smoothing: str
        Options are "low" for minimal smoothing (e.g. MIMF), "high" for large
        smoothing (e.g. GA), or "default" for medium smoothing for other cases

    Returns
    -------
    matched_labels : list
        List of alphabetic labels for each identified PSF
    """
    labels = string.ascii_uppercase[:18]

    # Determine boundaries of array
    x_min = min(x)
    x_max = max(x)
    y_min = min(y)
    y_max = max(y)

    if (x_max - x_min) > (y_max - y_min):
        # Horizontal orientation
        x_list = [1, 2, 2, 3, 3, 3, 4, 4, 5, 5, 6, 6, 7, 7, 7, 8, 8, 9]
        y_list = [3, 2, 4, 1, 3, 5, 2, 4, 1, 5, 2, 4, 1, 3, 5, 2, 4, 3]

        x_coords = np.linspace(x_min, x_max, 9)
        y_coords = np.linspace(y_min, y_max, 5)[::-1]
    else:
        # Vertical orientation
        x_list = [3, 2, 4, 1, 3, 5, 2, 4, 1, 5, 2, 4, 1, 3, 5, 2, 4, 3]
        y_list = [9, 8, 8, 7, 7, 7, 6, 6, 5, 5, 4, 4, 3, 3, 3, 2, 2, 1]

        x_coords = np.linspace(x_min, x_max, 5)
        y_coords = np.linspace(y_min, y_max, 9)[::-1]

    seg_coords = np.array([[x_coords[i_x - 1],
                            y_coords[i_y - 1]] for i_x, i_y in zip(x_list, y_list)])

    # Match actual blob coordinates to segment name
    matched_labels = []
    for x_pos, y_pos in zip(x, y):
        seg_distance = 2048
        for i_sc, sc in enumerate(seg_coords):
            x_distance = x_pos - sc[0]
            y_distance = y_pos - sc[1]
            distance = (x_distance**2 + y_distance**2)**0.5
            if distance < seg_distance:
                seg_distance = distance
                i_seg = i_sc
        matched_labels.append(labels[i_seg])

    if len(set(matched_labels)) != len(matched_labels) and smoothing == 'high':
        LOGGER.warning('Could not accurately map labels to segments. It will not '
                       'be possible to run fsw_file_writer.rewrite_prc using the '
                       'all_found_psfs*.txt file generated here.')

    return matched_labels


def get_car_data():
    """
    Generate and return a dictionary in the form of {car: {apt: #, observations: #}} for
    all the OTE and LOS CARs in commissioning. CARs and APT#s scraped from
    http://www.stsci.edu/ftp/presto/ops/public/jwst-pit-status.html and # of observations
    pulled from jwst_magic/jwst_magic/data/commissioning_activities.yaml
    """
    url = 'http://www.stsci.edu/ftp/presto/ops/public/jwst-pit-status.html'
    html_page = requests.get(url).text
    l = pd.read_html(html_page)
    df = l[0]
    df_set = df[df['Activity ID'].str.contains("OTE|LOS|NIRCam-005|NIRCam-035", case=False)]  # only include guiding CARs

    commissioning_dict = {car.lower(): str(int(apt))
                          for car, apt in zip(df_set['Activity ID'].values, df_set['Program'].values)}

    return commissioning_dict


def natural_keys(text):
    """
    alist.sort(key=natural_keys) sorts in human order
    http://nedbatchelder.com/blog/200712/human_sorting.html
    (from stack overflow:
    https://stackoverflow.com/questions/5967500/how-to-correctly-sort-a-string-with-a-number-inside)
    """

    # For sorting strings by numbers
    def atoi(text):
        return int(text) if text.isdigit() else text

    return [atoi(c) for c in re.split(r'(\d+)', text)]

def setup_yaml():
    """
    Set up a yaml file that will have an OrderedDict written to it
    """
    represent_dict_order = lambda self, data: self.represent_mapping('tag:yaml.org,2002:map', data.items())
    yaml.add_representer(OrderedDict, represent_dict_order)