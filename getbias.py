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

    # 0th read bias structure
    if bzp:
        bias0file = os.path.join(data_path, 'g{}bias0.fits'.format(guider))
        read0 = fits.getdata(bias0file)#.astype(np.uint16)
        bias += read0[xlow:xhigh, ylow:yhigh]

        bias[bias < 0] = 0.
        bias[bias > 40000] = 10000.

        bias *= 0.5 # 22 may 13 sth damp 0th read structure for ernie's testing


    # KTC imprints at reset
    if bktc:
        ktc = 10. * np.random.standard_normal((nramps, ny, nx))
        ktc_full = np.repeat(ktc, nreads, axis=0) #repeat the KTC for all reads
        bias += ktc_full

    # pedestal imprints at read
    if bp:
        nnx = 2048
        nny = 2048
        pedestal = np.zeros((nz, nny, nnx))
        for i in range(nz):
            pedestal[i, :, 0:511] = np.fix(np.round(25 * np.random.standard_normal()))
            pedestal[i, :, 512:1023] = np.fix(np.round(25 * np.random.standard_normal()))
            pedestal[i, :, 1024:1535] = np.fix(np.round(25 * np.random.standard_normal()))
            pedestal[i, :, 1536:2047] = np.fix(np.round(25 * np.random.standard_normal()))
        #if pedestal doesnt change throughout cube:
        #pedestal[:, :, 0:511] = np.fix(np.round(25 * np.random.standard_normal())
        #pedestal[:, :, 512:1023] = np.fix(np.round(25 * np.random.standard_normal())
        #pedestal[:, :, 1024:1535] = np.fix(np.round(25 * np.random.standard_normal())
        #pedestal[:], :, 1536:2047] = np.fix(np.round(25 * np.random.standard_normal())
        bias += pedestal[:, xlow:xhigh, ylow:yhigh]

    # rectify bias img
    bias[bias < 0] = 0.
    bias[bias > 65000] = 0.

    bias.astype(np.uint16)


    return bias
