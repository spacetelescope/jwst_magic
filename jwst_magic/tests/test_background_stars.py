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
import sys

from astropy.io import fits
JENKINS = 'jenkins' in os.getcwd()
if not JENKINS:
    from PyQt5 import QtCore
    from PyQt5.QtWidgets import QApplication
import pytest

from jwst_magic.convert_image.background_stars import add_background_stars
if not JENKINS:
    from jwst_magic.convert_image.background_stars import BackgroundStarsWindow

__location__ = os.path.realpath(os.path.join(os.getcwd(), os.path.dirname(__file__)))
NIRCAM_IM = os.path.join(__location__, 'data', 'nircam_data_1_ga.fits')

JENKINS = 'jenkins' in os.getcwd()

def test_add_background_stars():
    # Basic test to make sure add_background_stars works
    with fits.open(NIRCAM_IM) as hdulist:
        data = hdulist[1].data

    stars = {'x': [1500], 'y': [200], 'jmag': [10.0]}

    add_background_stars(data, stars, 11.0, 'FGS Magnitude', 1)

@pytest.mark.skipif(JENKINS, reason="Can't import PyQt5 on Jenkins server.")
def test_init_background_stars():
    """Make sure the background_stars GUI can launch without errors.
    """
    guider = 1
    jmag = 12
    qApp = QApplication(sys.argv)
    window = BackgroundStarsWindow(guider, jmag, qApp=qApp, in_master_GUI=False)

    # Schedule press of "Cancel" button
    QtCore.QTimer.singleShot(0, window.pushButton_cancel.clicked)
