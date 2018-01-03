'''
Apply different types of bias
'''
# STDLIB
import os

# THIRD PARTY
import numpy as np
from astropy.io import fits

#LOCAL
import utils
import log

LOCAL_PATH = os.path.dirname(os.path.realpath(__file__))
DATA_PATH = os.path.join(LOCAL_PATH, 'data')
BIASZERO_G1 = fits.getdata(os.path.join(DATA_PATH,
                                        'g1bias0.fits'))
BIASZERO_G2 = fits.getdata(os.path.join(DATA_PATH,
                                        'g2bias0.fits'))

def getbias(guider, xcoord, ycoord, nreads, nramps, imgsize):
    """
    For ID: nreads = 2, nramps = 2, nx = 2048, ny = 2048
    For acquisition 1: nreads = 2, nramps = 6, nx, ny = 128, 128
    For acquisition 2: nreads = 2, nramps = 5, nx, ny = 32, 32
    """
    if guider == 1:
        bias0 = np.copy(BIASZERO_G1)
    elif guider == 2:
        bias0 = np.copy(BIASZERO_G2)
    else:
        log.error("Do not recognize guider {}. Exiting.".format(guider))
        raise StandardError("Guider not recognized.")

    if imgsize == 2048:
        xlow, ylow = 0, 0
        xhigh, yhigh = imgsize, imgsize

    else:
        xlow = np.fix(xcoord) - imgsize / 2
        xhigh = np.fix(xcoord) + imgsize / 2
        ylow = np.fix(ycoord) - imgsize / 2
        yhigh = np.fix(ycoord) + imgsize / 2

    bias = np.zeros((nramps*nreads, imgsize, imgsize))

    ## 0th read bias structure
    bias += bias0[xlow:xhigh, ylow:yhigh]
    bias = utils.correct_image(bias, upper_threshold=40000., upper_limit=10000.)

    ## KTC imprints at reset
    ktc = 10. * np.random.random_sample((nramps, imgsize, imgsize))
    ktc_full = np.repeat(ktc, nreads, axis=0) #repeat the KTC for all reads
    bias += ktc_full

    ## Pedestal imprints at read
    pedestal = np.zeros((nramps*nreads, 2048, 2048)) #Create the full frame pedestal image
    for i in range(nramps*nreads):
        pedestal[i, :, 0:511] = np.fix(np.round(25 * np.random.random_sample()))
        pedestal[i, :, 512:1023] = np.fix(np.round(25 * np.random.random_sample()))
        pedestal[i, :, 1024:1535] = np.fix(np.round(25 * np.random.random_sample()))
        pedestal[i, :, 1536:2047] = np.fix(np.round(25 * np.random.random_sample()))
    bias += pedestal[:, xlow:xhigh, ylow:yhigh]

    ## rectify bias img
    bias = utils.correct_image(bias, upper_limit=0) #set any saturated pixels to zero
    bias.astype(np.uint16)

    return bias
