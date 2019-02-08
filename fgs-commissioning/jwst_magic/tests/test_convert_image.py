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
        data = hdulist[1].data

    # Transform sci to FGS1 raw
    image = convert_image_to_raw_fgs.transform_sci_to_fgs_raw(data, 1)
    sources = find_peaks(image, np.max(image) * .05, box_size=5)
    coords = [(x, y) for (x, y) in sources['x_peak', 'y_peak']]
    assert coords == [(993, 975), (1021, 975), (1050, 975), (1007, 999), (1061, 999), (1035, 1000), (981, 1002), (1216, 1023), (965, 1024), (993, 1025), (1078, 1026), (982, 1049), (1007, 1050), (1062, 1050), (1035, 1051)], \
        'Incorrect transformation from DMS/science frame to raw FGS1 frame'

    # Transform sci to FGS2 raw
    image = convert_image_to_raw_fgs.transform_sci_to_fgs_raw(data, 2)
    sources = find_peaks(image, np.max(image) * .05, box_size=5)
    coords = [(x, y) for (x, y) in sources['x_peak', 'y_peak']]
    assert coords == [(997, 975), (1026, 975), (1054, 975), (986, 999), (1040, 999), (1012, 1000), (1066, 1002), (831, 1023), (1082, 1024), (1054, 1025), (969, 1026), (1065, 1049), (985, 1050), (1040, 1050), (1012, 1051)], \
        'Incorrect transformation from DMS/science frame to raw FGS2 frame'


def test_transform_nircam_raw_to_fgs_raw():
    with fits.open(FGS_CMIMF_IM) as hdulist:
        data = hdulist[1].data

    # Transform NRCA3 raw to FGS1 raw
    image = convert_image_to_raw_fgs.transform_nircam_raw_to_fgs_raw(data, 'A3', 1)
    sources = find_peaks(image, np.max(image) * .05, box_size=5)
    coords = [(x, y) for (x, y) in sources['x_peak', 'y_peak']]
    assert coords == [(1035, 996), (1007, 997), (1062, 997), (982, 998), (1078, 1021), (993, 1022), (965, 1023), (1216, 1024), (981, 1045), (1035, 1047), (1007, 1048), (1061, 1048), (993, 1072), (1021, 1072), (1050, 1072)], \
        'Incorrect transformation from raw NRCA3 (thus also A1, A5, B2, B4) frame to raw FGS1 frame'

    # Transform NRCB1 raw to FGS1 raw
    image = convert_image_to_raw_fgs.transform_nircam_raw_to_fgs_raw(data, 'B1', 1)
    sources = find_peaks(image, np.max(image) * .05, box_size=5)
    coords = [(x, y) for (x, y) in sources['x_peak', 'y_peak']]
    assert coords == [(997, 975), (1026, 975), (1054, 975), (986, 999), (1040, 999), (1012, 1000), (1066, 1002), (831, 1023), (1082, 1024), (1054, 1025), (969, 1026), (1065, 1049), (985, 1050), (1040, 1050), (1012, 1051)], \
        'Incorrect transformation from raw NRCB1 (thus also A2, A4, B3, B5) frame to raw FGS1 frame'

    # Transform NRCA5 raw to FGS2 raw
    image = convert_image_to_raw_fgs.transform_nircam_raw_to_fgs_raw(data, 'A5', 2)
    sources = find_peaks(image, np.max(image) * .05, box_size=5)
    coords = [(x, y) for (x, y) in sources['x_peak', 'y_peak']]
    assert coords == [(1012, 996), (985, 997), (1040, 997), (1065, 998), (969, 1021), (1054, 1022), (1082, 1023), (831, 1024), (1066, 1045), (1012, 1047), (986, 1048), (1040, 1048), (997, 1072), (1026, 1072), (1054, 1072)], \
        'Incorrect transformation from raw NRCA5 (thus also A1, A3, B2, B4) frame to raw FGS2 frame'

    # Transform NRCB1 raw to FGS2 raw
    image = convert_image_to_raw_fgs.transform_nircam_raw_to_fgs_raw(data, 'B1', 2)
    sources = find_peaks(image, np.max(image) * .05, box_size=5)
    coords = [(x, y) for (x, y) in sources['x_peak', 'y_peak']]
    assert coords == [(993, 975), (1021, 975), (1050, 975), (1007, 999), (1061, 999), (1035, 1000), (981, 1002), (1216, 1023), (965, 1024), (993, 1025), (1078, 1026), (982, 1049), (1007, 1050), (1062, 1050), (1035, 1051)], \
        'Incorrect transformation from raw NRCB1 (thus also A2, A4, B3, B5) frame to raw FGS2 frame'

