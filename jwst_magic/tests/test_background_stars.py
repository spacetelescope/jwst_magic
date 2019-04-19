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
from PyQt5 import QtCore
from PyQt5.QtWidgets import QApplication

from jwst_magic import background_stars

__location__ = os.path.realpath(os.path.join(os.getcwd(), os.path.dirname(__file__)))
NIRCAM_IM = os.path.join(__location__, 'data', 'nircam_data_1_ga.fits')

JENKINS = 'jenkins' in os.getcwd()

def test_add_background_stars():
    # Basic test to make sure add_background_stars works
    with fits.open(NIRCAM_IM) as hdulist:
        data = hdulist[1].data

    stars = {'x': [1500], 'y': [200], 'jmag': [10.0]}

    background_stars.add_background_stars(data, stars, 11.0, 'FGS Magnitude', 1)


def test_init_background_stars():
    """Make sure the background_stars GUI can launch without errors.
    """
    # Can't currently run this test on the Jenkins server
    if JENKINS:
        return

    guider = 1
    jmag = 12
    qApp = QApplication(sys.argv)
    window = background_stars.BackgroundStarsWindow(guider, jmag, qApp=qApp, in_master_GUI=False)

    # Schedule press of "Cancel" button
    QtCore.QTimer.singleShot(0, window.pushButton_cancel.clicked)
