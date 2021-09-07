"""
This module contains helper functions to read in different values
related to the brightness of stars in FGS (countrate, FGS magnitude,
etc.) and converts them using the FGS Count Rate Module found here:
https://github.com/spacetelescope/jwst-fgs-countrate

Authors
-------
    - Keira Brooks
    - Lauren Chambers
    - Kathryn St.Laurent
    - Shannon Osborne

Use
---
    This module can be imported in a Python shell as such:
    ::
        from jwst_magic.convert_image import renormalize

"""
# Standard Library Imports
import logging

# Third Party Imports
import fgscountrate
import numpy as np

# Constants
FGS_ZERO_POINT = 29.057
J_ZERO_POINT = 0.90
CONVERSION_FACTOR = 3.1418185

# Start logger
LOGGER = logging.getLogger(__name__)

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# MAIN FUNCTION
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

def convert_to_countrate_fgsmag(value, unit, guider, gs_catalog=None):
    """
    Convert an input normalization unit and value (e.g. a guide star id)
    into a count rate and an FGS magnitude.

    value : str or float
        Specifies the Guide Star ID or the count rate/magnitude to which
        to normalize.
    unit : str,
        Specifies the unit of norm_value ("FGS Magnitude", "FGS countrate",
        or "Guide Star ID")
    guider : int
        Guider number (1 or 2)
    gs_catalog : str, optional
        Guide star catalog version to query. E.g. 'GSC242'. None will use
        the default catalog as defined in teh FGS Count Rate Module.
    """
    expected_type_dict = {'fgs countrate': (float, int),
                          'fgs magnitude': (float, int),
                          'guide star id': (str,)
                          }

    if not isinstance(value, expected_type_dict[unit.lower()]):
        types = ' or '.join([d.__name__ for d in expected_type_dict[unit.lower()]])
        raise TypeError(
            'Mismatch: Normalization value for unit of "{}" should be of type(s) "{}". User passed: {}'.format(
                unit, types, value))

    if unit.lower() == 'fgs countrate':
        if value < 1000:
            LOGGER.warning("Normalize Calculation: Small count rate has been requested. This may be an error.")
        fgs_countrate = value
        fgs_mag = fgscountrate.convert_cr_to_fgs_mag(fgs_countrate, guider)

    elif unit.lower() == 'fgs magnitude':
        fgs_mag = value
        fgs_countrate = fgscountrate.convert_fgs_mag_to_cr(fgs_mag, guider)

    elif unit.lower() == 'guide star id':
        LOGGER.info("Normalize Calculation: Using GSID {} to normalize image.".format(value))
        fgs = fgscountrate.FGSCountrate(guide_star_id=value, guider=guider)
        fgs_countrate, _, fgs_mag, _ = fgs.query_fgs_countrate_magnitude(catalog=gs_catalog)

    return fgs_countrate, fgs_mag


def query_guide_star_catalog(gs_id, gs_catalog=None):
    """Query the quide star catalog using a
    guide star ID and return the RA and DEC of the
    guide star
    """
    # Use Guide Star ID to get RA/DEC using default GSC in fgscountrate module
    data_frame = fgscountrate.query_gsc(gs_id=gs_id, catalog=gs_catalog)

    # Check there's only 1 line in the GSC with this GS ID
    if len(data_frame) == 1:
        gsc_series = data_frame.iloc[0]
    else:
        raise ValueError("This Guide Star ID points to multiple lines in catalog")

    # Pull RA and DEC
    ra = gsc_series['ra']
    dec = gsc_series['dec']

    return ra, dec
