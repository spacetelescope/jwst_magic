'''Collection of unit tests to verify the correct function of the
convert_image modules in the JWST MAGIC package.

Authors
-------
    - Lauren Chambers

Use
---
    These tests can be run via the command line (omit the ``-s`` to
    suppress verbose output to ``stdout``):

    ::
        pytest -s test_convert_image.py
"""'''
import glob
import os

import numpy as np
import pytest
from astropy.io import fits

from jwst_magic import utils, coordinate_transforms
from jwst_magic.convert_image import convert_image_to_raw_fgs, renormalize
from jwst_magic.star_selector import select_psfs
from jwst_magic.fsw_file_writer import buildfgssteps

FGS_DATA = glob.glob('data/fgs_data*.fits')
NIRCAM_DATA = glob.glob('data/nircam_data*.fits')
TEST_DATA = glob.glob('data/*data*.fits')

OUT_PATH = os.path.join(os.getcwd())


class TestRenormalize():
    def test_NormalizeToCounts(self):

        # Initialize the object
        norm = renormalize.NormalizeToCounts(11.0, 'FGS Magnitude', 1)

        # Ensure all the attributes are assigned
        assert norm.value == 11.0
        assert norm.unit == 'FGS Magnitude'
        assert norm.guider == 1

        # Test starting with magnitude
        mag_to_counts = norm.to_counts()
        assert mag_to_counts == 30160035.219588336
        mag_to_mag = norm.to_fgs_mag()
        assert mag_to_mag == 11.0

        # Test starting with counts
        norm = renormalize.NormalizeToCounts(1e7, 'FGS Counts', 1)
        counts_to_counts = norm.to_counts()
        assert counts_to_counts == 1e7
        counts_to_mag = norm.to_fgs_mag()
        assert counts_to_mag == 12.198579610870993

    def test_electrons(self):

        assert renormalize.counts_to_electrons(1e7, 1) == 17400000.0
        assert renormalize.counts_to_electrons(1e7, 2) == 15700000.0
        assert renormalize.electrons_to_counts(1e7, 1) == 5747126.436781609
        assert renormalize.electrons_to_counts(1e7, 2) == 6369426.751592357

    def test_fgs_mag(self):

        assert renormalize.counts_to_fgs_mag(1e7, 1) == 12.198579610870993
        assert renormalize.counts_to_fgs_mag(1e7, 2) == 12.310203600554409
        assert renormalize.fgs_mag_to_counts(11.0, 1) == 30160035.219588336
        assert renormalize.fgs_mag_to_counts(11.0, 2) == 33425771.517250765

    def test_jmag(self):

        assert renormalize.fgs_mag_to_j_mag(11.0) == -18.956999999999997
        assert renormalize.fgs_mag_to_j_mag(11.0) == -18.956999999999997
        assert renormalize.fgs_counts_to_j_mag(1e7, 1) == -17.758420389129004
        assert renormalize.fgs_counts_to_j_mag(1e7, 2) == -17.646796399445588
        assert renormalize.j_mag_to_fgs_counts(11.0, 1) == 3.137847582229968e-05
        assert renormalize.j_mag_to_fgs_counts(11.0, 2) == 3.4776145178854426e-05

    def test_gain(self):

        assert renormalize.find_gain(1) == 1.74
        assert renormalize.find_gain(2) == 1.57

class TestConvertImageToRawFGS():
    with fits.open(FGS_DATA[0]) as hdulist:
        data = hdulist[1].data

    def test_apply_coarse_pointing_filter(self):
        data_gauss = convert_image_to_raw_fgs.apply_coarse_pointing_filter(self.data, 0, 0.031)
        assert data_gauss == self.data

        data_gauss = convert_image_to_raw_fgs.apply_coarse_pointing_filter(self.data, 0.15, 0.031)
        assert data_gauss != self.data
