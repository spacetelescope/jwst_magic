'''Collection of unit tests to verify the correct function of the FGS
Commissioning Tools

Authors
-------
    - Lauren Chambers

Use
---
    These tests can be run via the command line (omit the ``-s`` to
    suppress verbose output to ``stdout``):

    ::
        pytest -s test_tool.py
"""'''
import glob
import os

import numpy as np
import pytest
from astropy.io import fits

from jwst_magic import utils, coordinate_transforms
from jwst_magic.convert_image import convert_image_to_raw_fgs
from jwst_magic.star_selector import select_psfs
from jwst_magic.fsw_file_writer import buildfgssteps

FGS_DATA = glob.glob('fgs_data*.fits')
NIRCAM_DATA = glob.glob('nircam_data*.fits')
TEST_DATA = glob.glob('*data*.fits')

OUT_PATH = os.path.join(os.getcwd())

class RunningTheTool():
    def __init__(self, input_im):
        print('\n', input_im)

        # Get parameters to run through functions
        if 'ga' in input_im.lower():
            self.global_alignment = True
        else:
            self.global_alignment = False

        self.root = utils.make_root(None, input_im)

        header = fits.getheader(input_im, ext=0)
        self.guider = utils.get_guider(header)

        print('* * * GENERATING DATA * * *')
        # If provided a NIRCam image, convert it to FGS
        if 'nircam' in input_im:
            nircam = True
        else:
            nircam = False

        data = convert_image_to_raw_fgs.convert_im(input_im, self.guider,
                                                   nircam=nircam, out_dir=OUT_PATH)

        # Find peaks in data and create a reg file with randomly selected PSFs
        select_psfs.create_reg_file(data, self.root, self.guider,
                                    global_alignment=self.global_alignment,
                                    testing=True, out_dir=OUT_PATH)

        # Write all steps accordingly
        steps = ['ID', 'ACQ1', 'ACQ2', 'TRK', 'LOSTRK']
        for step in steps:
            buildfgssteps.BuildFGSSteps(data, self.guider, self.root, step,
                                        out_dir=OUT_PATH)

    def prc_coordinates_matches_peak_locations(self):
        """
        Compares the strips image and the ID prc file from one data set. Test to ensure
        that the location of the peaks in the strips image match the commanded
        locations in the prc.
        """
        print('* * * TESTING DATA * * *')
        # Create full frame array from strips fits
        strips_filename = os.path.join('out', self.root, 'dhas',
                                       '{}_G{}_IDstrips.fits'.format(self.root,
                                                                     self.guider))
        strips_data = fits.open(strips_filename)[0].data
        ff = np.zeros((4, 2048, 2048))
        for i in range(36):
            for j in range(4):
                strip = strips_data[j + i * 4]
                offset = 0
                overlap = 8
                height = 64
                y_start = offset + (i * (height - overlap))
                y_end = offset + (i * (height - overlap)) + height
                ff[j, y_start:y_end] = strip
        self.strips_data_cds = ff[1] - ff[0]

        self.images_dont_have_unit_problems()

        # (Only need the coords but reading out the whole thing anyway)
        strips_params = select_psfs.manual_star_selection(self.strips_data_cds,
                                                          self.global_alignment,
                                                          testing=True)
        strips_coords = strips_params[1]

        # Read what was commanded to the prc
        prc_file = os.path.join('out', self.root, 'dhas',
                                '{}_G{}_ID.prc'.format(self.root, self.guider))
        prc_coords = []
        with open(prc_file) as f:
            lines = f.readlines()
            for line in lines:
                if line.startswith('@IFGS_GUIDESTAR') or \
                   line.startswith('@IFGS_REFSTAR'):
                    x_dhas = float(line.split(',')[2])
                    y_dhas = float(line.split(',')[3])
                    x_raw, y_raw = coordinate_transforms.DHAS2Raw(x_dhas, y_dhas,
                                                                  self.guider)
                    # Incorporate the 12-pixel offset of the strips (DHAS default param)
                    prc_coords.append((x_raw, y_raw - 12))

        print('Coordinates from strips.fits: ', strips_coords)
        print('Coordinates from ID.prc: ', prc_coords)

        # Make sure that each commanded prc coordinate is within ~5 pixels
        # of a peak from the strips image
        test_radius = 5
        for (x_prc, y_prc) in prc_coords:
            found = False
            for (x_strips, y_strips) in strips_coords:
                if abs(x_prc - x_strips) < test_radius and \
                   abs(y_prc - y_strips) < test_radius:
                        found = True
                        break
            assert found, 'Commanded coordinate in ID.prc, ({}, {}), does not fall within {} pixels of a PSF in the strips image.'.format(x_prc, y_prc, test_radius)

    def images_dont_have_unit_problems(self):
        '''For now, just test of the mean of a CDS image is above 1000 counts.
        (cases with unsigned int problems have many pixels > 1e5.)
        In the future this should probably be more sophisticated.'''

        assert np.mean(self.strips_data_cds.flatten()) < 1e3, 'The average value of pixels is sufficiently high to suggest that there is a unsigned integer problem happening (negative values becoming arbitrarily high values).'

def test_tool():
    '''Run the whole tool for each test .fits in the tests/ directory'''
    for input_im in TEST_DATA:
        print(input_im)
        run = RunningTheTool(input_im)
        run.prc_coordinates_matches_peak_locations()
