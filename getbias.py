'''
Apply different types of bias
'''
# STDLIB
import os

# Third Party
from astropy.io import fits
import numpy as np

LOCAL_PATH = os.path.dirname(os.path.realpath(__file__))

def getbias(guider, xcoord, ycoord, nreads, nramps, nx, ny, bzp=True,
            bktc=True, bp=True):
    """
    For ID: nreads = 2, nramps = 2, nx = 2048, ny = 2048
    For acquisition 1: nreads = 2, nramps = 6, nx, ny = 128, 128
    For acquisition 2: nreads = 2, nramps = 5, nx, ny = 32, 32
    """
    # Establish location of 'data' directory. Includes *magicHdrImg.fits,
    #*bias.fits, etc. If data path not given, script will assume that the
    #'data' directory lives in the same directory as this script
    data_path = os.path.join(LOCAL_PATH, 'data')

    nz = nramps * nreads
    if (nx == 2048) and (ny == 2048):
        xlow, ylow = 0, 0
        xhigh, yhigh = nx, ny

    else:
        xlow = np.fix(xcoord) - nx / 2
        xhigh = np.fix(xcoord) + nx / 2
        ylow = np.fix(ycoord) - ny / 2
        yhigh = np.fix(ycoord) + ny / 2

    bias = np.zeros((nz, ny, nx))

    # 0th read bias structure (present in every read)
    if bzp:
        bias0file = os.path.join(data_path, 'g{}bias0.fits'.format(guider))
        read0 = fits.getdata(bias0file)#.astype(np.uint16)

        if bias.shape[1] % 2 == 0:
            bias += read0[xlow:xhigh, ylow:yhigh]
        elif bias.shape[1] % 2 == 1:  # If odd, add another pixel (?!?)
            bias += read0[xlow:xhigh + 1, ylow:yhigh + 1]

        bias[bias < 0] = 0.
        bias[bias > 40000] = 10000.

        # bias *= 0.5  # 22 may 13 sth damp 0th read structure for ernie's testing

    # KTC imprints at reset (uniform within each ramp)
    if bktc:
        ktc = 10. * np.random.random_sample((nramps, ny, nx))
        ktc_full = np.repeat(ktc, nreads, axis=0) #repeat the KTC for all reads
        bias += ktc_full

    # Pedestal imprints at reset (uniform within each ramp)
    if bp:
        pedestal = np.zeros((nz, 2048, 2048))
        for iramp in range(nramps):
            ped_noise = np.fix(25 * np.random.standard_normal(size=4))
            pedestal[iramp * nramps:(iramp + 1) * nramps, :, 0:511] = ped_noise[0]
            pedestal[iramp * nramps:(iramp + 1) * nramps, :, 512:1023] = ped_noise[1]
            pedestal[iramp * nramps:(iramp + 1) * nramps, :, 1024:1535] = ped_noise[2]
            pedestal[iramp * nramps:(iramp + 1) * nramps, :, 1536:2047] = ped_noise[3]

        if bias.shape[1] % 2 == 0:
            bias += pedestal[:, xlow:xhigh, ylow:yhigh]
        elif bias.shape[1] % 2 == 1:  # If odd, add another pixel (?!?)
            bias += pedestal[:, xlow:xhigh + 1, ylow:yhigh + 1]

    # rectify bias img
    bias[bias < 0] = 0.
    bias[bias > 65000] = 0.

    bias.astype(np.uint16)

    return bias
