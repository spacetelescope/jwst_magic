"""Collection of unit tests to verify the correct function of the Image
Conversion tool.

TODO
(In the future, this will) Runs the following tests:
    1. Test all of convert_image

Authors
-------
    - Keira Brooks
    - Lauren Chambers

Use
---
    ::
        pytest test_convert_image.py

Notes
-----

At this moment, convert image can take in a NIRCam or FGS image either created
with WebbPSF postage stamps or created with the ITM simulator.

The following are ways to run convert_image:

    from jwst_magic.convert_image import convert_image_to_raw
    convert_image_to_raw.convert_im(input_im, guider, root)

"""

import os

from astropy.io import fits
import numpy as np
from photutils import find_peaks

from jwst_magic.convert_image import convert_image_to_raw_fgs


__location__ = os.path.realpath(os.path.join(os.getcwd(), os.path.dirname(__file__)))
ROOT = "test_convertim"
NIRCAM_IM = os.path.join(__location__, 'data', 'nircam_data_1_ga.fits')
FGS_GA_IM = os.path.join(__location__, 'data', 'fgs_data_1_beforelos2ga.fits')
FGS_CMIMF_IM = os.path.join(__location__, 'data', 'fgs_data_2_cmimf.fits')


def test_nircam_itm():
    # Test that the normalization with FGS Mag = 12.0 returns a specific value
    # for a NIRCam image

    # Guider 1
    data_g1 = convert_image_to_raw_fgs.convert_im(NIRCAM_IM, guider=1, root=ROOT,
                                                  nircam=True, norm_value=12.0,
                                                  norm_unit="FGS Magnitude", itm=True)


    assert round(np.max(data_g1), 2) == 11963.36

    # Guider 2
    data_g2 = convert_image_to_raw_fgs.convert_im(NIRCAM_IM, guider=2, root=ROOT,
                                                  nircam=True, norm_value=12.0,
                                                  norm_unit="FGS Magnitude", itm=True)

    assert round(np.max(data_g2), 2) == 13258.76


def test_fgs_itm():
    # Test that the normalization with FGS Mag = 12.0 returns a specific value
    # for an FGS image

    # Guider 1
    data_g1 = convert_image_to_raw_fgs.convert_im(FGS_GA_IM, guider=1, root=ROOT,
                                                  nircam=True, norm_value=12.0,
                                                  norm_unit="FGS Magnitude", itm=True)

    assert round(np.max(data_g1), 2) == 29333.71

    # Guider 2
    data_g2 = convert_image_to_raw_fgs.convert_im(FGS_GA_IM, guider=2, root=ROOT,
                                                  nircam=True, norm_value=12.0,
                                                  norm_unit="FGS Magnitude", itm=True)

    assert round(np.max(data_g2), 2) == 32509.97

# This test needs to be created after we have a better understanding of how
# to test logging in out code
# def test_logging():
#     header = fits.getheader(FGS_IM)
#     assert header['ORIGIN'].strip() == 'ITM'
#     guider = 2
#     fgs_mag = 12.0
#     norm = renormalize.NormalizeToCountrate(value=fgs_mag, unit='FGS Magnitude', guider=guider)
#     fgs_countrate = norm.to_countrate()
#     _ = convert_image_to_raw_fgs.convert_im(FGS_IM, guider=guider, root=ROOT,
#                                             nircam=False, norm_value=12.0,
#                                             norm_unit="FGS Magnitude", itm=False)
#     l = LogCapture()
#
#     # with LogCapture() as l:
#     #     logger = logging.getLogger()
#     #     logger.warning("Deprecation Warning: This is an ITM image, setting itm flag to 'True'")
#     #mock_logger = logging.getLogger('mock_logger')
#     #mock_logger.warning.assert_called_with("Deprecation Warning: This is an ITM image, setting itm flag to 'True'")
#
#     l.check(('jwst_magic.convert_image.convert_image_to_raw_fgs', 'INFO', "Image Conversion: " +
#                 "Beginning image conversion to guider {} FGS image".format(guider)))
#     l.check(('jwst_magic.convert_image.convert_image_to_raw_fgs', 'INFO', "Image Conversion: Input image is expected to be in units of ADU/sec (countrate)"))
#     l.check(('jwst_magic.convert_image.convert_image_to_raw_fgs', 'INFO', "Image Conversion: Data provided in science/DMS frame; rotating to raw FGS frame."))
#     l.check(('jwst_magic.convert_image.convert_image_to_raw_fgs', 'INFO', "Image Conversion: This is an ITM image."))
#     l.check(('jwst_magic.convert_image.convert_image_to_raw_fgs', 'INFO', "Image Conversion: Normalizing to FGS Magnitude of {:.1f} ({} FGS Countrate)".format(fgs_mag, fgs_countrate)))
#     l.check(('jwst_magic.convert_image.convert_image_to_raw_fgs', 'WARNING', "Deprecation Warning: This is an ITM image, setting itm flag to 'True'"))
#     l.uninstall()

def test_transform_sci_to_fgs_raw():
    with fits.open(FGS_CMIMF_IM) as hdulist:
        hdulist.info()
        data = hdulist[1].data

    sources = find_peaks(data, np.max(data)*.05, box_size=5)
    print(sources)
    coords = [(x, y) for (x, y) in sources['x_peak', 'y_peak']]

    print(coords)

    # Transform NRCA3 raw to FGS1 raw

    # Transform NRCB1 raw to FGS1 raw

    # Transform NRCA5 raw to FGS2 raw

    # Transform NRCB1 raw to FGS2 raw

    # Transform NRCA3 sci to FGS1 raw

    # Transform NRCB3 sci to FGS2 raw

    # Transform FGS1 sci to FGS1 raw

    # Transform FGS1 sci to FGS2 raw

    # Transform FGS2 sci to FGS1 raw

    # Transform FGS2 sci to FGS2 raw
