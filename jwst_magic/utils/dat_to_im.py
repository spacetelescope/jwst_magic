# View a DAT file, either with matplotlib or writing out to FITS

# Standard Library Imports
import re
import sys

# Third Party Imports
from matplotlib.colors import LogNorm
import matplotlib.pyplot as plt
import numpy as np

# Local Imports
from jwst_magic.utils import utils


def dat_to_array(dat_file):
    with open(dat_file) as f:
        data = f.readlines()

    if data[0][4] == ' ':
        # ASCII Hex
        cds = True
        data = re.split(' ', data[0])[:-1]
        data = [float.fromhex(num) for num in data]

        if 'ACQ1' in dat_file:
            data = np.reshape(data, (12, 128, 128))
        elif 'ACQ2' in dat_file:
            data = np.reshape(data, (10, 32, 32))
        elif 'IDstrips' in dat_file:
            data = np.reshape(data, (144, 64, 2048))
        else:
            print('Unrecognized file type; cannot reshape flattened array.')

    else:
        # ASCII Float
        # Only TRK files are ASCII floats, generally
        cds = False
        if '\n' in data[0]:
            # If carriage returns
            data = np.array([re.split('  |\n', row)[1:-1] for row in data]).astype(float)
        else:
            data = np.array(data[0].split()).astype(float)

            if 'LOSTRK' in dat_file:
                data = np.reshape(data, (255, 255))
            else:
                print('Unrecognized file type; cannot reshape flattened array.')

    return data, cds


def dat_to_fits(file):
    data, cds = dat_to_array(file)
    fits_filename = file[:-3] + 'fits'
    utils.write_fits(fits_filename, data)
    print('Saved .fits file to: ', fits_filename)


if __name__ == '__main__':
    if len(sys.argv) < 2:
        print("USAGE:")
        print("./dat_to_im.py FILENAME\n")
        sys.exit(1)

    file = sys.argv[1]

    data, cds = dat_to_array(file)
    if cds:
        print("Displaying CDS image")
        plt.imshow(data[1] - data[0], norm=LogNorm())
    else:
        plt.imshow(data, norm=LogNorm())
    plt.show()
