"""Collection of unit tests to verify the correct function of the background
stars module.

Authors
-------
    - Lauren Chambers

Use
---
    ::
        pytest test_background_stars.py
"""
import os

from astropy.io import fits

from jwst_magic import background_stars

__location__ = os.path.realpath(os.path.join(os.getcwd(), os.path.dirname(__file__)))
NIRCAM_IM = os.path.join(__location__, 'data', 'nircam_data_1_ga.fits')

def test_add_background_stars():
    # Basic test to make sure add_background_stars works
    with fits.open(NIRCAM_IM) as hdulist:
        data = hdulist[1].data

    stars = {'x': [1500], 'y': [200], 'jmag': [10.0]}

    background_stars.add_background_stars(data, stars, 11.0, 'FGS Magnitude', 1)