'''
Apply different types of bias
'''
# STDLIB
import os

# Third Party
from astropy.io import fits
import numpy as np

LOCAL_PATH = os.path.dirname(os.path.realpath(__file__))
PACKAGE_PATH = os.path.split(LOCAL_PATH)[0]

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
    data_path = os.path.join(PACKAGE_PATH, 'data')

    nz = nramps * nreads
    if (nx == 2048) and (ny == 2048):
        xlow, ylow = 0, 0
        xhigh, yhigh = nx, ny

    else:
        xlow = int(np.fix(xcoord) - nx / 2)
        xhigh = int(np.fix(xcoord) + nx / 2)
        ylow = int(np.fix(ycoord) - ny / 2)
        yhigh = int(np.fix(ycoord) + ny / 2)

    bias = np.zeros((nz, ny, nx))

    # 0th read bias structure (present in every read) - - - - - - - - - - - - -
    if bzp:
        bias0file = os.path.join(data_path, 'g{}bias0.fits'.format(guider))
        read0 = fits.getdata(bias0file)#.astype(np.uint16)

        bias += read0[xlow:xhigh, ylow:yhigh]

        bias[bias < 0] = 0.
        bias[bias > 40000] = 10000.

        # bias *= 0.5  # 22 may 13 sth damp 0th read structure for ernie's testing

    # KTC imprints at reset (uniform within each ramp) - - - - - - - - - - - - -
    if bktc:
        ktc = 10. * np.random.random_sample((nramps, ny, nx))
        ktc_full = np.repeat(ktc, nreads, axis=0) #repeat the KTC for all reads
        bias += ktc_full

    # Pedestal imprints at reset (uniform within each ramp) - - - - - - - - - -
    if bp:
        # For full frame images (ID)
        if nx == 2048:
            pedestal = np.zeros((nz, ny, nx))
            ped_noise = np.fix(25 * np.random.standard_normal(size=(nramps, 4)))
            for iramp, ped in zip(range(nramps), ped_noise):
                pedestal[iramp * nreads:(iramp + 1) * nreads, :, 0:511] = ped[0]
                pedestal[iramp * nreads:(iramp + 1) * nreads, :, 512:1023] = ped[1]
                pedestal[iramp * nreads:(iramp + 1) * nreads, :, 1024:1535] = ped[2]
                pedestal[iramp * nreads:(iramp + 1) * nreads, :, 1536:2047] = ped[3]

            pedestal = pedestal[:, xlow:xhigh, ylow:yhigh]

        # For subarrays
        else:
            # Determine if subarray spans multiple pedestals:
            spanning = False
            for border in [512, 1024, 1536]:
                if (border - nx) < xlow < border:
                    spanning = border

            # If subarray spans multiple pedestals:
            if type(spanning) == int:
                pedestal = np.zeros((nz, ny, nx))

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
                pedestal = np.repeat(ped_noise, ny, axis=1)
                pedestal = np.repeat(pedestal, nx, axis=2)
                pedestal = np.repeat(pedestal, nreads, axis=0)

        bias += pedestal

    # rectify bias img
    bias[bias < 0] = 0.
    bias[bias > 65000] = 0.

    bias.astype(np.uint16)

    return bias
