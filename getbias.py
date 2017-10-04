# STDLIB
import os

# Third Party
from astropy.io import fits
import numpy as np

def getbias(guider, x, y, nreads, nramps, nx, ny, bzp=True,
            bktc=True, bp=True, data_path=None):
    """
    For ID: nreads = 2, nramps = 2, nx = 2048, ny = 2048
    For acquisition 1: nreads = 2, nramps = 6, nx, ny = 128, 128
    For acquisition 2: nreads = 2, nramps = 5, nx, ny = 32, 32
    """
    # Establish location of 'data' directory. Includes *magicHdrImg.fits,
    #*bias.fits, etc. If data path not given, script will assume that the
    #'data' directory lives in the same directory as this script
    local_path = os.path.dirname(os.path.realpath(__file__))

    if data_path is None:
        data_path = os.path.join(local_path, 'data')
    else:
        data_path = data_path

    nz = nramps*nreads
    if (nx == 2048) and (ny == 2048):
        x1, y1 = 0, 0
        x2 = nx
        y2 = ny

    else:
        x1 = np.fix(x) - nx / 2
        x2 = np.fix(x) + nx / 2
        y1 = np.fix(y) - ny / 2
        y2 = np.fix(y) + ny / 2

    bias = np.zeros((nz, ny, nx))

    # 0th read bias structure
    if bzp:
        bias0file = os.path.join(data_path, 'g{}bias0.fits'.format(guider))
        read0 = fits.open(bias0file)[0].data#.astype(np.uint16)
        bias += read0[x1:x2, y1:y2]

        bias[bias < 0] = 0.
        bias[bias > 40000] = 10000.

        bias *= 0.5 # 22 may 13 sth damp 0th read structure for ernie's testing


    # KTC imprints at reset
    if bktc:
        ktc = 10. * np.random.standard_normal((nramps, ny, nx))
        iz = 0
        for iramp in range(nramps):
            for iread in range(nreads):
                bias[iz] += ktc[iramp]
                iz += 1

    # pedestal imprints at read
    if bp:
        nnx = 2048
        nny = 2048
        pedestal = np.zeros((nz, nny, nnx))
        tmp = np.fix(np.round(25*np.random.standard_normal((4, nz))))
        for iz in range(nz):
            pedestal[iz, :, 0:511] = tmp[0, iz]
            pedestal[iz, :, 512:1023] = tmp[1, iz]
            pedestal[iz, :, 1024:1535] = tmp[2, iz]
            pedestal[iz, :, 1536:2047] = tmp[3, iz]
            bias[iz] += pedestal[iz, x1:x2, y1:y2]

    # rectify bias img
    bias[bias < 0] = 0.
    bias[bias > 65000] = 0.

    bias.astype(np.uint16)


    return bias
