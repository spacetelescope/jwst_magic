'''For given FGS counts, J mag or Jab mag, get the others.

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
        from jwst_fgs_commissioning_tools.convert_image import counts_to_jmag

Notes
-----
    There is also some functionality for this in Pysynphot so that
    should be looked into at some point for consistency with other
    systems
'''
# Imports
# Third PARTY
import numpy as np

# GLOBAL VARIABLES
FGS_ZERO_POINT = 29.057
J_ZERO_POINT = 0.90
CONVERSION_FACTOR = 3.1418185

class NormalizeToCounts(object):
    '''
    Input the user-defined value and unit (FGS Counts, FGS Magnitude or J Magnitude)
    and convert into FGS Counts.
    '''
    def __init__(self, value, unit, guider):
        self.value = value
        self.unit = unit
        self.guider = guider

    def to_counts(self):
        if self.unit == 'FGS Counts':
            return self.value
        elif self.unit == 'FGS Magnitude':
            return fgs_mag_to_counts(self.value, self.guider)

    def to_fgs_mag(self):
        if self.unit == 'FGS Magnitude':
            return self.value
        elif self.unit == 'FGS Counts':
            return counts_to_fgs_mag(self.value, self.guider)

# Convert between FGS Counts and Electrons
def counts_to_electrons(counts, guider):
    '''Convert count rate to electrons/s'''
    return counts * find_gain(guider)

def electrons_to_counts(electrons, guider):
    '''Convert electrons/s to count rate'''
    return electrons / find_gain(guider)

# Convert between FGS Counts andf FGS Magnitude
def counts_to_fgs_mag(counts, guider):
    '''Convert FGS counts to FGS magnitude '''
    electrons = counts_to_electrons(counts, guider)
    fgs_mag = -2.5 * np.log10(electrons/CONVERSION_FACTOR) + FGS_ZERO_POINT
    return  fgs_mag

def fgs_mag_to_counts(fgs_mag, guider):
    '''Convert FGS Magnitude to FGS counts '''
    electrons = 10**((fgs_mag - FGS_ZERO_POINT)/-2.5) * CONVERSION_FACTOR
    counts = electrons_to_counts(electrons, guider)
    return  counts

def fgs_mag_to_j_mag(fgs_mag):
    '''Convert FGS Magnitude to J Magnitude '''
    j_mag = fgs_mag - (FGS_ZERO_POINT + J_ZERO_POINT)
    return j_mag

def j_mag_to_fgs_mag(j_mag):
    '''Convert J Magnitude to FGS Magnitude '''
    fgs_mag = j_mag + (FGS_ZERO_POINT + J_ZERO_POINT)
    return fgs_mag


def fgs_counts_to_j_mag(counts, guider):
    fgs_mag = counts_to_fgs_mag(counts, guider)
    j_mag = fgs_mag_to_j_mag(fgs_mag)
    return j_mag

def j_mag_to_fgs_counts(j_mag, guider):
    fgs_mag = j_mag_to_fgs_mag(j_mag)
    counts = fgs_mag_to_counts(fgs_mag, guider)
    return counts


def find_gain(guider):
    '''Find the gain for each guider to convert from ADU/sec to e-/sec

    Parameters
    ----------
    guider : int
        Guider number (1 or 2)

    Returns
    -------
    conversion : float
        Appropriate conversion factor from ADU/sec to e-/sec
    '''
    if guider == 1:
        gain = 1.55
    else:
        gain = 1.81
    return gain
