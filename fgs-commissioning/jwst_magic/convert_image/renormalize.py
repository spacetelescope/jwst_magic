"""For given FGS countrate, J mag or Jab mag, calculate the others.

Convert between FGS countrate, J magnitude, and J_ab magnitude. This
module is very basic and is only for a single bandpass; look at
count_rate.f in Sherie Holfeltz's code for the procedure for more
bandpasses.

Authors
-------
    - Keira Brooks
    - Lauren Chambers
    - Kathryn St.Laurent

Use
---
    This module can be imported in a Python shell as such:
    ::
        from jwst_magic.convert_image import renormalize

Notes
-----
    There is also some functionality for this in Pysynphot so that
    should be looked into at some point for consistency with other
    systems
"""
# Standard Library Imports
import logging

# Third Party Imports
import numpy as np

# Local Imports
from .. import utils

# Constants
FGS_ZERO_POINT = 29.057
J_ZERO_POINT = 0.90
CONVERSION_FACTOR = 3.1418185

# Start logger
LOGGER = logging.getLogger(__name__)

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# MAIN CLASS
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


class NormalizeToCountrate(object):
    """Input the user-defined value and unit (FGS Countrate, FGS Magnitude
    or J Magnitude) and convert into FGS Countrate.
    """
    def __init__(self, value, unit, guider):
        """Initialize the class.
        """
        self.value = value
        self.unit = unit
        self.guider = guider

    def to_countrate(self):
        """Convert self.value to FGS countrate

        Returns
        -------
        float
            Converted value in FGS countrate
        """
        if self.unit.lower() == 'fgs countrate':
            return self.value
        elif self.unit.lower() == 'fgs magnitude':
            return fgs_mag_to_countrate(self.value, self.guider)
        else:
            raise ValueError("Unknown unit:" + self.unit)

    def to_fgs_mag(self):
        """Convert self.value to FGS magnitude

        Returns
        -------
        float
            Converted value in FGS magnitude
        """
        if self.unit.lower() == 'fgs magnitude':
            return self.value
        elif self.unit.lower() == 'fgs countrate':
            return countrate_to_fgs_mag(self.value, self.guider)
        else:
            raise ValueError("Unknown unit:" + self.unit)

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# CONVERSION FUNCTIONS
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


def countrate_to_electrons(countrate, guider):
    """Convert count rate to electrons per second

    Parameters
    ----------
    countrate : float
        Input value in FGS countrate
    guider : int
        Guider number (1 or 2)

    Returns
    -------
    float
        Value converted to electrons per second
    """
    return countrate * find_gain(guider)


def electrons_to_countrate(electrons, guider):
    """Convert electrons per second to count rate

    Parameters
    ----------
    electrons : float
        Input value in electrons per second
    guider : int
        Guider number (1 or 2)

    Returns
    -------
    float
        Value converted to FGS countrate
    """
    return electrons / find_gain(guider)


def countrate_to_fgs_mag(countrate, guider):
    """Convert FGS countrate to FGS magnitude

    Parameters
    ----------
    countrate : float
        Input value in FGS countrate
    guider : int
        Guider number (1 or 2)

    Returns
    -------
    fgs_mag : float
        Value converted to FGS magnitude
    """
    electrons = countrate_to_electrons(countrate, guider)
    fgs_mag = -2.5 * np.log10(electrons/CONVERSION_FACTOR) + FGS_ZERO_POINT
    return fgs_mag


def fgs_mag_to_countrate(fgs_mag, guider):
    """Convert FGS magnitude to FGS countrate

    Parameters
    ----------
    fgs_mag : float
        Input value in FGS magnitude
    guider : int
        Guider number (1 or 2)

    Returns
    -------
    countrate : float
        Value converted to FGS countrate
    """
    electrons = 10**((fgs_mag - FGS_ZERO_POINT)/-2.5) * CONVERSION_FACTOR
    countrate = electrons_to_countrate(electrons, guider)
    return countrate


def fgs_mag_to_j_mag(fgs_mag):
    """Convert FGS magnitude to J magnitude

    Parameters
    ----------
    fgs_mag : float
        Input value in FGS magnitude

    Returns
    -------
    j_mag : float
        Value converted to J magnitude
    """
    j_mag = fgs_mag - (FGS_ZERO_POINT + J_ZERO_POINT)
    return j_mag


def fgs_countrate_to_j_mag(countrate, guider):
    """Convert FGS countrate to J magnitude

    Parameters
    ----------
    countrate : float
        Input value in FGS countrate
    guider : int
        Guider number (1 or 2)

    Returns
    -------
    j_mag : float
        Value converted to J magnitude
    """
    fgs_mag = countrate_to_fgs_mag(countrate, guider)
    j_mag = fgs_mag_to_j_mag(fgs_mag)
    return j_mag


def j_mag_to_fgs_countrate(j_mag, guider):
    """Convert J magnitude to FGS countrate

    Parameters
    ----------
    j_mag : float
        Input value in FGS magnitude
    guider : int
        Guider number (1 or 2)

    Returns
    -------
    countrate : float
        Value converted to FGS countrate
    """
    fgs_mag = j_mag + (FGS_ZERO_POINT + J_ZERO_POINT)
    countrate = fgs_mag_to_countrate(fgs_mag, guider)
    return countrate


def find_gain(guider):
    """Find the gain for each guider to convert from ADU/sec to e-/sec

    Parameters
    ----------
    guider : int
        Guider number (1 or 2)

    Returns
    -------
    conversion : float
        Appropriate conversion factor from ADU/sec to e-/sec
    """
    if guider == 1:
        gain = 1.74
    else:
        gain = 1.57
    return gain
