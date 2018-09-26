"""For given FGS counts, J mag or Jab mag, calculate the others.

Convert between FGS counts, J magnitude, and J_ab magnitude. This
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

# Constants
FGS_ZERO_POINT = 29.057
J_ZERO_POINT = 0.90
CONVERSION_FACTOR = 3.1418185

# Start logger
LOGGER = logging.getLogger(__name__)

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# MAIN CLASS
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


class NormalizeToCounts(object):
    """Input the user-defined value and unit (FGS Counts, FGS Magnitude
    or J Magnitude) and convert into FGS Counts.
    """
    def __init__(self, value, unit, guider, logger_passed=False):
        """Initialize the class.
        """
        # Set up logger
        if not logger_passed:
            utils.create_logger_from_yaml(__name__, root=root, level='DEBUG')

        self.value = value
        self.unit = unit
        self.guider = guider
        LOGGER.warning('* * * AN IMAGE THAT IS NORMALIZED TO COUNTS MIGHT HAVE INCORRECT UNITS DUE TO IT BEING MULTIPLIED BY THE FRAME TIME IN create_img_arrays() * * *')

    def to_counts(self):
        """Convert self.value to FGS counts

        Returns
        -------
        float
            Converted value in FGS counts
        """
        if self.unit == 'FGS Counts':
            return self.value
        elif self.unit == 'FGS Magnitude':
            return fgs_mag_to_counts(self.value, self.guider)

    def to_fgs_mag(self):
        """Convert self.value to FGS magnitude

        Returns
        -------
        float
            Converted value in FGS magnitude
        """
        if self.unit == 'FGS Magnitude':
            return self.value
        elif self.unit == 'FGS Counts':
            return counts_to_fgs_mag(self.value, self.guider)

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# CONVERSION FUNCTIONS
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


def counts_to_electrons(counts, guider):
    """Convert count rate to electrons per second

    Parameters
    ----------
    counts : float
        Input value in FGS counts
    guider : int
        Guider number (1 or 2)

    Returns
    -------
    float
        Value converted to electrons per second
    """
    return counts * find_gain(guider)


def electrons_to_counts(electrons, guider):
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
        Value converted to FGS counts
    """
    return electrons / find_gain(guider)


def counts_to_fgs_mag(counts, guider):
    """Convert FGS counts to FGS magnitude

    Parameters
    ----------
    counts : float
        Input value in FGS counts
    guider : int
        Guider number (1 or 2)

    Returns
    -------
    fgs_mag : float
        Value converted to FGS magnitude
    """
    electrons = counts_to_electrons(counts, guider)
    fgs_mag = -2.5 * np.log10(electrons/CONVERSION_FACTOR) + FGS_ZERO_POINT
    return fgs_mag


def fgs_mag_to_counts(fgs_mag, guider):
    """Convert FGS magnitude to FGS counts

    Parameters
    ----------
    fgs_mag : float
        Input value in FGS magnitude
    guider : int
        Guider number (1 or 2)

    Returns
    -------
    counts : float
        Value converted to FGS counts
    """
    electrons = 10**((fgs_mag - FGS_ZERO_POINT)/-2.5) * CONVERSION_FACTOR
    counts = electrons_to_counts(electrons, guider)
    return counts


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


def fgs_counts_to_j_mag(counts, guider):
    """Convert FGS counts to J magnitude

    Parameters
    ----------
    counts : float
        Input value in FGS counts
    guider : int
        Guider number (1 or 2)

    Returns
    -------
    j_mag : float
        Value converted to J magnitude
    """
    fgs_mag = counts_to_fgs_mag(counts, guider)
    j_mag = fgs_mag_to_j_mag(fgs_mag)
    return j_mag


def j_mag_to_fgs_counts(j_mag, guider):
    """Convert J magnitude to FGS counts

    Parameters
    ----------
    j_mag : float
        Input value in FGS magnitude
    guider : int
        Guider number (1 or 2)

    Returns
    -------
    counts : float
        Value converted to FGS counts
    """
    fgs_mag = j_mag + (FGS_ZERO_POINT + J_ZERO_POINT)
    counts = fgs_mag_to_counts(fgs_mag, guider)
    return counts


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
