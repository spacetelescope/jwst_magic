'''For given FGS counts, J mag or Jab mag, get the others.

Convert between FGS counts, J magnitude, and J_ab magnitude. This
module is very basic and is only for a single bandpass; look at
count_rate.f in Sherie Holfeltz's code for the procedure for more
bandpasses.

Authors
-------
    - Keira Brooks
    - Lauren Chambers

Use
---
    This module can be imported in a Python shell as such:
    ::
        from jwst_magic.convert_image import counts_to_jmag

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

def countrate_to_e(countrate, guider):
    '''Convert count rate to electrons/s'''
    conversion = find_conversion(guider)
    return countrate * conversion

def e_to_counts(electrons, guider):
    '''Convert electrons to counts'''
    conversion = find_conversion(guider)
    return electrons/conversion

def jab_to_jmag(jab):
    '''Convert J ab to J magnitude'''
    return jab-0.90

def jmag_to_jab(jmag):
    '''Convert J magnitude to  J ab '''
    return jmag+0.90

def fgs_counts_to_jmag(countrate, guider):
    '''Convert FGS counts to J magnitude given the count rate'''
    electrons = countrate_to_e(countrate, guider)
    jab = -2.5 * np.log10(electrons/3.1418185) + 29.057
    return  jab_to_jmag(jab)

def jmag_to_fgs_counts(jmag, guider):
    '''Convert J magnitude to FGS counts given the J magnitude'''
    jab = jmag_to_jab(jmag)
    electrons = 10**((jab - 29.057)/-2.5) * 3.1418185
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
