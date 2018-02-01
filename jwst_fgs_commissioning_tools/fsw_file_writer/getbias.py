'''
Apply different types of bias
'''
# STDLIB
import os

# THIRD PARTY
import numpy as np
from astropy.io import fits

#LOCAL
from jwst_fgs_commissioning_tools import utils
from jwst_fgs_commissioning_tools import log

FSW_PATH = os.path.dirname(os.path.realpath(__file__))
PACKAGE_PATH = os.path.split(FSW_PATH)[0]
DATA_PATH = os.path.join(PACKAGE_PATH, 'data')
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
        xlow = int(np.fix(xcoord) - imgsize / 2)
        xhigh = int(np.fix(xcoord) + imgsize / 2)
        ylow = int(np.fix(ycoord) - imgsize / 2)
        yhigh = int(np.fix(ycoord) + imgsize / 2)

    bias = np.zeros((nramps * nreads, imgsize, imgsize))

    ## 0th read bias structure
    bias += bias0[xlow:xhigh, ylow:yhigh]
    bias = utils.correct_image(bias, upper_threshold=40000., upper_limit=10000.)

    ## KTC imprints at reset
    ktc = 10. * np.random.random_sample((nramps, imgsize, imgsize))
    ktc_full = np.repeat(ktc, nreads, axis=0) #repeat the KTC for all reads
    bias += ktc_full


    # Pedestal imprints at reset (uniform within each ramp) - - - - - - - - - -
    if imgsize == 2048: #for full frame images - ID
        pedestal = np.zeros((nramps * nreads, 2048, 2048)) #Create the full frame pedestal image
        for i in range(nramps * nreads):
            pedestal[i, :, 0:511] = np.fix(np.round(25 * np.random.random_sample()))
            pedestal[i, :, 512:1023] = np.fix(np.round(25 * np.random.random_sample()))
            pedestal[i, :, 1024:1535] = np.fix(np.round(25 * np.random.random_sample()))
            pedestal[i, :, 1536:2047] = np.fix(np.round(25 * np.random.random_sample()))

    # For subarrays
    else:
        # Determine if subarray spans multiple pedestals:
        spanning = False
        for border in [512, 1024, 1536]:
            if (border - imgsize) < xlow < border:
                spanning = border

        # If subarray spans multiple pedestals:
        if type(spanning) == int:
            pedestal = np.zeros((nramps * nreads, imgsize, imgsize))

            # Determine transition x value
            xborder = spanning - 1 - xlow

            ped_noise = np.fix(25 * np.random.standard_normal(size=(nramps, 2)))
            for iramp, ped in zip(range(nramps), ped_noise):
                pedestal[iramp * nreads:(iramp + 1) * nreads, :, :xborder] = ped[0]
                pedestal[iramp * nreads:(iramp + 1) * nreads, :, xborder:] = ped[1]

        # Else if subarray is completely within just one pedestal:
        elif not spanning:
            ped_noise = np.zeros((nramps, 1, 1))
            ped_noise[:, 0, 0] = np.fix(25 * np.random.standard_normal(size=nramps))

            # Resize array to match pedestal array
            pedestal = np.repeat(ped_noise, imgsize, axis=1)
            pedestal = np.repeat(pedestal, imgsize, axis=2)
            pedestal = np.repeat(pedestal, nreads, axis=0)

    bias += pedestal[:, xlow:xhigh, ylow:yhigh]

    ## rectify bias img
    bias = utils.correct_image(bias, upper_limit=0) #set any saturated pixels to zero
    bias.astype(np.uint16)

    return bias
