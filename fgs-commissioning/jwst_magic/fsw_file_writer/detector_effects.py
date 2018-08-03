"""Generate various types of bias and noise for FGS images.

Generate an array of appropriate dimensions that includes the following
bias and noise features: zeroth-read bias structure, kTc noise, read
noise, and amp-to-amp pedestal (not currently included).

Authors
-------
    - Keira Brooks
    - Lauren Chambers

Use
---
    This module can be used as such:
    ::
        from jwst_magic.fsw_file_writer import detector_effects
        det_eff = detector_effects.FGSDetectorEffects(guider, xcoord,
            ycoord, nreads, nramps, imgsize)
        bias = det_eff.add_detector_effects()

    Required arguments:
        ``guider`` - guider number (1 or 2)
        ``xcoord`` - X coordinate of guide star (pixels)
        ``ycoord`` - Y coordinate of guide star (pixels)
        ``nreads`` - number of reads in the ramp
        ``nramps`` - number of ramps in the FGS step
        ``imgsize`` - dimension of the array (pixels)

Notes
-----
    Adapted from code originally written by Sherie Holfeltz in IDL.
"""

# Standard Library Imports
import os

# Third Party Imports
import numpy as np
import yaml
from astropy.io import fits

# Local
from .. import utils

# Paths
__location__ = os.path.realpath(os.path.join(os.getcwd(), os.path.dirname(__file__)))
PACKAGE_PATH = os.path.split(__location__)[0]
DATA_PATH = os.path.join(PACKAGE_PATH, 'data')
BIASZERO_G1 = fits.getdata(os.path.join(DATA_PATH, 'g1bias0.fits'))
BIASZERO_G2 = fits.getdata(os.path.join(DATA_PATH, 'g2bias0.fits'))
READ_NOISE = os.path.join(DATA_PATH, 'readnoise.yaml')


class FGSDetectorEffects():
    """Fetch the bias file for the specified guider, crop to the
    appropriate array size, multiply to have the appropriate number of
    ramps and reads, and add kTc  and read noise.

    Parameters
    ----------
    guider : int
        Guider number (1 or 2)
    xcoord : int
        X coordinate of guide star (pixels)
    ycoord : int
        Y coordinate of guide star (pixels)
    nreads : int
        Number of reads in the ramp
    nramps : int
        Numer of ramps in the integraion
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
    def __init__(self, guider, xcoord, ycoord, nreads, nramps, imgsize):
        # Check that the guider is valid
        if int(guider) not in [1, 2]:
            raise ValueError("Guider {} not recognized.".format(guider))

        # Initialize attributes
        self.guider = guider
        self.gs_coords = (xcoord, ycoord)
        self.nreads = nreads
        self.nramps = nramps
        self.imgsize = imgsize

        # Create an empty array of the appropriate size
        self.bias = np.zeros((nramps * nreads, imgsize, imgsize))

        # Determine the bounds of the needed array
        self.array_bounds = self.get_subarray_location()

    def add_detector_effects(self):
        """Create an array with zeroth read bias, read noise, and kTc
        noise effects included

        Returns
        -------
        bias : numpy array
            Array with bias data with dimension (nramps * nreads,
            imgsize, imgsize) that contains all detector effects
        """
        # Add bias and noise
        self.add_zeroth_read_bias()
        self.add_read_noise()
        self.add_ktc_noise()
        # bias = add_pedestal()

        # Correct bias image for saturated or negative pixels
        self.bias = utils.correct_image(self.bias, upper_limit=0)

        # *** Does the bias need to be in the uint16 type?? ***
        # self.bias = self.bias.astype(np.uint16)

        return self.bias

    def add_ktc_noise(self):
        """Add kTc noise, which imprints at reset, to the bias.
        """
        ktc = 10. * np.random.random_sample((self.nramps, self.imgsize, self.imgsize))
        # Repeat the KTC for all reads
        ktc_full = np.repeat(ktc, self.nreads, axis=0)
        self.bias += ktc_full

    def add_pedestal(self):
        """Add pedestal, which imprints at reset, to the bias.
        """
        # For full frame images (ID, CAL)
        if self.imgsize == 2048:
            # Create the full frame pedestal image
            pedestal = np.zeros((self.nramps * self.nreads, 2048, 2048))
            for i in range(self.nramps * self.nreads):
                pedestal[i, :, 0:511] = np.fix(np.round(25 * np.random.random_sample()))
                pedestal[i, :, 512:1023] = np.fix(np.round(25 * np.random.random_sample()))
                pedestal[i, :, 1024:1535] = np.fix(np.round(25 * np.random.random_sample()))
                pedestal[i, :, 1536:2047] = np.fix(np.round(25 * np.random.random_sample()))

        # For subarrays
        else:
            xlow, xhigh, ylow, yhigh = self.array_bounds

            # Determine if subarray spans multiple pedestals:
            spanning = False
            for border in [512, 1024, 1536]:
                if (border - self.imgsize) < xlow < border:
                    spanning = border

            # If subarray spans multiple pedestals:
            if type(spanning) == int:
                pedestal = np.zeros((self.nramps * self.nreads, self.imgsize, self.imgsize))

                # Determine transition x value
                xborder = spanning - 1 - xlow

                ped_noise = np.fix(25 * np.random.standard_normal(size=(self.nramps, 2)))
                for iramp, ped in zip(range(self.nramps), ped_noise):
                    pedestal[iramp * self.nreads:(iramp + 1) * self.nreads, :, :xborder] = ped[0]
                    pedestal[iramp * self.nreads:(iramp + 1) * self.nreads, :, xborder:] = ped[1]

            # Else if subarray is completely within just one pedestal:
            elif not spanning:
                ped_noise = np.zeros((self.nramps, 1, 1))
                ped_noise[:, 0, 0] = np.fix(25 * np.random.standard_normal(size=self.nramps))

                # Resize array to match pedestal array
                pedestal = np.repeat(ped_noise, self.imgsize, axis=1)
                pedestal = np.repeat(pedestal, self.imgsize, axis=2)
                pedestal = np.repeat(pedestal, self.nreads, axis=0)

        self.bias += pedestal

    def add_read_noise(self):
        """Add read noise to every frame.
        """
        # Load all the read noise values from the yaml
        with open(READ_NOISE) as f:
            read_noise_dict = yaml.load(f.read())

        # Get the read noise value for the current step
        array = 'full' if self.imgsize == 2048 else 'subarray'
        read_noise = read_noise_dict['guider{}'.format(self.guider)][array]

        # Add normally distributed read noise to the bias
        # *** What should the standard deviation be??? ***
        std_dev = read_noise * 0.05
        self.bias += np.random.normal(loc=read_noise, scale=std_dev,
                                      size=np.shape(self.bias)).astype(int)

    def add_zeroth_read_bias(self):
        """Add zeroth read bias structure to every frame.
        """
        # Open the zeroth read bias structure file
        bias_file = os.path.join(DATA_PATH, 'g{}bias0.fits'.format(self.guider))
        bias0 = np.copy(fits.getdata(bias_file))

        xlow, xhigh, ylow, yhigh = self.array_bounds

        # Get zeroth read bias structure from FITS file
        self.bias += bias0[xlow:xhigh, ylow:yhigh]
        self.bias = utils.correct_image(self.bias, upper_threshold=40000., upper_limit=10000.)

    def get_subarray_location(self):
        """Get the bounds of the subarray for the given step.

        Returns
        -------
        array_bounds : tuple
            Tuple containing the xlow, xhigh, ylow, and yhigh values
            defining the bounds of the subarray on the full-frame array
        """
        if self.imgsize == 2048:
            xlow, ylow = 0, 0
            xhigh, yhigh = self.imgsize, self.imgsize

        else:
            xcoord, ycoord = self.gs_coords
            xlow = int(np.fix(xcoord) - self.imgsize / 2)
            xhigh = int(np.fix(xcoord) + self.imgsize / 2)
            ylow = int(np.fix(ycoord) - self.imgsize / 2)
            yhigh = int(np.fix(ycoord) + self.imgsize / 2)

        array_bounds = (xlow, xhigh, ylow, yhigh)
        return array_bounds