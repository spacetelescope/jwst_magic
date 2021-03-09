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
    expected_type_dict = {
        'fgs countrate': 'float or int',
        'fgs magnitude': 'float or int',
        'guide star id': 'str',
    }

    if unit.lower() == 'fgs countrate' and isinstance(value, (float, int)):
        if value < 1000:
            LOGGER.warning("Normalize Calculation: Small count rate has been requested. This may be an error.")
        fgs_countrate = value
        fgs_mag = fgscountrate.convert_cr_to_fgs_mag(fgs_countrate, guider)

    elif unit.lower() == 'fgs magnitude' and isinstance(value, (float, int)):
        fgs_mag = value
        fgs_countrate = fgscountrate.convert_fgs_mag_to_cr(fgs_mag, guider)

    elif unit.lower() == 'guide star id' and isinstance(value, str):
        LOGGER.info("Normalize Calculation: Using GSID {} to normalize image.".format(value))
        fgs = fgscountrate.FGSCountrate(guide_star_id=value, guider=guider)
        fgs_countrate, _, fgs_mag, _ = fgs.query_fgs_countrate_magnitude(catalog=gs_catalog)

    else:
        raise TypeError('Mismatch: Normalization value for unit of "{}" should be of type "{}". User passed: {}'.format(
            unit, expected_type_dict[unit.lower()], value))

    return fgs_countrate, fgs_mag
