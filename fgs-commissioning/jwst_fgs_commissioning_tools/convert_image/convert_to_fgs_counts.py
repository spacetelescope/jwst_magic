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

# Third Party
import matplotlib.pyplot as plt
from matplotlib import cycler
import matplotlib
import numpy as np

# Set plots to default to matplotlib 2.0 colors
matplotlib.rcParams['axes.prop_cycle'] = cycler(u'color', ['#1f77b4',
                                                           '#ff7f0e',
                                                           '#2ca02c',
                                                           '#d62728',
                                                           '#9467bd',
                                                           '#8c564b',
                                                           '#e377c2',
                                                           '#7f7f7f',
                                                           '#bcbd22',
                                                           '#17becf'])

def find_conversion(guider):
    '''Find the conversion factor from ADU/sec to e-/sec for a given
    guider

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
        conversion = 1.55
    else:
        conversion = 1.81
    return conversion

def to_e(countrate, guider):
    '''Convert count rate to electrons/s'''
    conversion = find_conversion(guider)
    return countrate * conversion

def to_counts(electrons, guider):
    '''Convert electrons to counts'''
    conversion = find_conversion(guider)
    return electrons/conversion

def to_jmag(value, unit):
    '''Convert Jab to J magnitude'''
    return jab_mag-0.90

def to_fgs_mag(value, unit):
    '''Convert J magnitude to  FGS magnitude '''
    if unit == 'j_mag':
        fgs_mag = jmag+0.90
    return jmag+0.90

def fgs_mag(value, unit, guider):
    '''Convert FGS counts to FGS magnitude given the count rate'''
    electrons = countrate_to_e(countrate, guider)
    fgs_mag = -2.5 * np.log10(electrons/3.1418185) + 29.057
    return  jab_mag_to_jmag(fgs_mag)



def jmag_to_fgs_counts(jmag, guider):
    '''Convert J magnitude to FGS counts given the J magnitude'''
    fgs_mag = jmag_to_jab_mag(jmag)
    electrons = 10**((fgs_mag - 29.057)/-2.5) * 3.1418185
    return e_to_counts(electrons, guider)

def plot_jmag_v_counts(jmags=None):
    '''For a range of J magnitudes, plot the relation with FGS
    counts (defaults to Jmags = (5,18)) for both guider 1 and 2.

    Parameters
    ----------
    jmags : tuple, optional
        Start and stop of J magnitude array
    '''
    if jmags is None:
        jmag_start = 5
        jmag_end = 18
    else:
        jmag_start = jmags[0]
        jmag_end = jmags[1]

    jmags = np.arange(jmag_start, jmag_end, 0.1)

    crs1 = []
    crs2 = []
    for j in jmags:
        crs1.append(jmag_to_fgs_counts(j, 1))
        crs2.append(jmag_to_fgs_counts(j, 2))

    plt.figure(figsize=(10, 10))
    plt.plot(jmags, crs1, label='Guider 1')
    plt.plot(jmags, crs2, label='Guider 2')
    plt.grid(True)
    plt.xlabel('J mag')
    plt.ylabel('Total FGS Counts')
    plt.title('J mag vs FGS counts')
    plt.legend()
    plt.show()

class NormalizeToFgs(object):
    '''
    Input the user-defined value and unit (FGS Counts, FGS Magnitude or J Magnitude)
    and convert into FGS Counts.

    #TODO turn into module that convert from anything to anything.
    '''
    def __init__(self, value, unit, guider):
        self.value = value
        self.unit = unit
        self.guider = guider

        self.fgs_zero_point = 29.057
        self.jab_zero_point = 0.90
        self.conversion = find_conversion(guider)


    def to_counts(self):
        if self.unit == 'FGS Counts':
            return self.value
        elif self.unit == 'FGS Magnitude':
            return self.fgs_mag_to_counts()
        elif self.unit == 'J Magnitude':
            return self.jmag_to_fgs_counts()

    def to_fgs_mag(self):
        if self.unit == 'FGS Magnitude':
            return self.value
        elif self.unit == 'FGS Counts':
            return self.counts_to_fgs_mag()
        elif self.unit == 'J Magnitude':
            return self.jmag_to_fgs_mag()

    def to_j_mag(self):
        if self.unit == 'J Magnitude':
            return self.value
        elif self.unit == 'FGS Counts':
            return self.counts_to_jmag()
        elif self.unit == 'FGS Magnitude':
            return self.fgs_mag_to_jmag()

    # Convert between FGS Counts and Electrons
    def counts_to_electrons(self, counts):
        '''Convert count rate to electrons/s'''
        return counts * self.conversion

    def electrons_to_counts(self, electrons):
        '''Convert electrons/s to count rate'''
        return electrons / self.conversion

    # Convert between FGS Counts andf FGS Magnitude
    def counts_to_fgs_mag(self):
        '''Convert FGS counts to FGS magnitude given the count rate'''
        electrons = self.counts_to_electrons(self.value)
        fgs_mag = -2.5 * np.log10(electrons/3.1418185) + self.fgs_zero_point
        return  fgs_mag

    def fgs_mag_to_counts(self):
        '''Convert FGS counts to FGS magnitude given the count rate'''
        electrons = 10**((fgs_mag - self.fgs_zero_point)/-2.5) * 3.1418185
        counts = self.electrons_to_counts(electrons)
        return  counts

    # Convert between counts and J magnitude
    def counts_to_jmag(self):
        '''Convert FGS counts to FGS magnitude given the count rate'''
        electrons = self.counts_to_electrons(self.value)
        jab = -2.5 * np.log10(electrons/3.1418185)
        jmag = self.jab_to_jmag
        return  jmag

    def jmag_to_counts(self):
        '''Convert FGS counts to FGS magnitude given the count rate'''
        jab = self.jmag_to_jab(self.value)
        electrons = 10**(jab/-2.5) * 3.1418185
        counts = self.electrons_to_counts(electrons)
        return  counts

    # Convert between J AB and FGS Magnitude
    def jmag_to_fgs_mag(self):
        jab = self.jmag_to_jab(self.value)
        fgs_mag = jab + self.fgs_zero_point
        return fgs_mag

    def fgs_mag_to_jmag(self):
        jab = self.value - self.fgs_zero_point
        jmag = self.jab_to_jmag(jab)
        return jmag

    # Convert between J AB and J Magnitude
    def jab_to_jmag(self, jab):
        '''Convert Jab to J magnitude'''
        return jab - self.jab_zero_point

    def jmag_to_jab(self, jmag):
        '''Convert J magnitude to  FGS magnitude '''
        return jmag + self.jab_zero_point

    # Convert between FGS Counts andf J Magnitude
    def counts_to_jmag(self):
        '''Convert FGS counts to FGS magnitude given the count rate'''
        electrons = self.counts_to_electrons(self.value)
        jab = -2.5 * np.log10(electrons/3.1418185)
        jmag = self.jab_to_jmag(jab)
        return  jmag

    def jmag_to_fgs_counts(self):
        '''Convert J magnitude to FGS counts given the J magnitude'''
        jab = self.jmag_to_jab(self.value)
        electrons = 10**(jab/-2.5) * 3.1418185
        counts = self.electrons_to_counts(electrons)
        return counts
