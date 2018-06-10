"""Apply different types of bias to FGS images.

Authors
-------
    - Keira Brooks
    - Lauren Chambers

Use
---
    This module can be used as such:
    ::
        from jwst_magic.fsw_file_writer import getbias
        bias = getbias.getbias(guider, xcoord, ycoord, nreads,
            nramps, imgsize)

    Required arguments:
        ```guider``` - guider number (1 or 2)
        ```xcoord``` - X coordinate of guide star (pixels)
        ```ycoord``` - Y coordinate of guide star (pixels)
        ```nreads``` - number of reads in the ramp
        ```nramps``` - number of ramps in the FGS step
        ```imgsize``` - dimension of the array (pixels)

Notes
-----
    Adapted from code originally written by Sherie Holfeltz in IDL.
"""

# Standard Library Imports
import os

# Third Party Imports
import numpy as np
from astropy.io import fits

# Local
from .. import utils

# Paths
__location__ = os.path.realpath(os.path.join(os.getcwd(), os.path.dirname(__file__)))
PACKAGE_PATH = os.path.split(__location__)[0]
DATA_PATH = os.path.join(PACKAGE_PATH, 'data')
BIASZERO_G1 = fits.getdata(os.path.join(DATA_PATH, 'g1bias0.fits'))
BIASZERO_G2 = fits.getdata(os.path.join(DATA_PATH, 'g2bias0.fits'))


def getbias(guider, xcoord, ycoord, nreads, nramps, imgsize):
    """Fetch the bias file for the specified guider, crop to the
    appropriate array size, multiply to have the appropriate number of
    ramps and reads, and add kTc noise.

    For ID: nreads = 2, nramps = 2, nx = 2048, ny = 2048
    For ACQ1: nreads = 2, nramps = 6, nx, ny = 128, 128
    For ACQ2: nreads = 2, nramps = 5, nx, ny = 32, 32

    Parameters
    ----------
    guider : int
        Guider number (1 or 2)
    xcoord : int
        X coordinate of guide star (pixels)
    ycoord : int
        Y coordinate of guide star (pixels)
    nramps : int
        Numer of ramps in the integraion
    nreads : int
        Number of reads in the ramp
    imgsize : int
        Dimension of the image array (pixels)

    Returns
    -------
    bias : numpy array
        Array with bias data with dimension (nramps * nreads, imgsize, imgsize)

    Raises
    ------
    ValueError
        Invalid guider number provided
    """
    # Check that the guider is valid
    if int(guider) not in [1, 2]:
        raise ValueError("Guider {} not recognized.".format(guider))

    # Open the zeroth read bias structure file
    bias_file = os.path.join(DATA_PATH, 'g{}bias0.fits'.format(guider))
    bias0 = np.copy(fits.getdata(bias_file))

    # Determine the size of the needed array
    if imgsize == 2048:
        xlow, ylow = 0, 0
        xhigh, yhigh = imgsize, imgsize

    else:
        xlow = int(np.fix(xcoord) - imgsize / 2)
        xhigh = int(np.fix(xcoord) + imgsize / 2)
        ylow = int(np.fix(ycoord) - imgsize / 2)
        yhigh = int(np.fix(ycoord) + imgsize / 2)

    bias = np.zeros((nramps * nreads, imgsize, imgsize))

    # Get zeroth read bias structure from FITS file
    bias += bias0[xlow:xhigh, ylow:yhigh]
    bias = utils.correct_image(bias, upper_threshold=40000., upper_limit=10000.)

    # Add kTc noise (imprints at reset)
    ktc = 10. * np.random.random_sample((nramps, imgsize, imgsize))
    ktc_full = np.repeat(ktc, nreads, axis=0) #repeat the KTC for all reads
    bias += ktc_full

    # # Pedestal imprints at reset (uniform within each ramp) - - - - - - - - - -
    # if imgsize == 2048: #for full frame images - ID
    #     pedestal = np.zeros((nramps * nreads, 2048, 2048)) #Create the full frame pedestal image
    #     for i in range(nramps * nreads):
    #         pedestal[i, :, 0:511] = np.fix(np.round(25 * np.random.random_sample()))
    #         pedestal[i, :, 512:1023] = np.fix(np.round(25 * np.random.random_sample()))
    #         pedestal[i, :, 1024:1535] = np.fix(np.round(25 * np.random.random_sample()))
    #         pedestal[i, :, 1536:2047] = np.fix(np.round(25 * np.random.random_sample()))

    # # For subarrays
    # else:
    #     # Determine if subarray spans multiple pedestals:
    #     spanning = False
    #     for border in [512, 1024, 1536]:
    #         if (border - imgsize) < xlow < border:
    #             spanning = border

    #     # If subarray spans multiple pedestals:
    #     if type(spanning) == int:
    #         pedestal = np.zeros((nramps * nreads, imgsize, imgsize))

    #         # Determine transition x value
    #         xborder = spanning - 1 - xlow

    #         ped_noise = np.fix(25 * np.random.standard_normal(size=(nramps, 2)))
    #         for iramp, ped in zip(range(nramps), ped_noise):
    #             pedestal[iramp * nreads:(iramp + 1) * nreads, :, :xborder] = ped[0]
    #             pedestal[iramp * nreads:(iramp + 1) * nreads, :, xborder:] = ped[1]

    #     # Else if subarray is completely within just one pedestal:
    #     elif not spanning:
    #         ped_noise = np.zeros((nramps, 1, 1))
    #         ped_noise[:, 0, 0] = np.fix(25 * np.random.standard_normal(size=nramps))

    #         # Resize array to match pedestal array
    #         pedestal = np.repeat(ped_noise, imgsize, axis=1)
    #         pedestal = np.repeat(pedestal, imgsize, axis=2)
    #         pedestal = np.repeat(pedestal, nreads, axis=0)

    # bias += pedestal

    # Correct bias image for saturated or negative pixels
    bias = utils.correct_image(bias, upper_limit=0)
    bias.astype(np.uint16)

    return bias
