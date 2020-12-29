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
    if unit.lower() == 'fgs countrate':
        if isinstance(value, (float, int)):  # If already a count rate, do nothing
            if value < 1000:
                LOGGER.warning("Normalize Calculation: Small count rate has been requested. This may be an error.")
            fgs_countrate = value
            fgs_mag = fgscountrate.convert_cr_to_fgs_mag(fgs_countrate, guider)
        else:
            LOGGER.error("Normalize Calculation: Value type does not match expectation for unit type.")
            raise TypeError("Mismatched normalization value ({}) and expected type ({}) for unit of {}".format(
                value, "str", unit))

    else:  # If a GSID or magnitude, pass back gsid
        # Sneaky string replacement for FGS Magnitude
        unit = 'fgs magnitude' if unit.lower() == 'fgs mag (10, 11, 12, 13, or 14)' else unit.lower()
        gsid = check_norm_value_unit(value, unit)
        LOGGER.info("Normalize Calculation: Using GSID {} to normalize image.".format(gsid))
        fgs = fgscountrate.FGSCountrate(guide_star_id=gsid, guider=guider)
        fgs_countrate, _, fgs_mag, _ = fgs.query_fgs_countrate_magnitude(catalog=gs_catalog)

    return fgs_countrate, fgs_mag


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# HELPER FUNCTIONS
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

def check_norm_value_unit(norm_value, norm_unit):
    """
    For the case when the norm_unit is not "FGS countrate", get back the
    guide star ID for the norm_value given.

    Parameters
    ----------
    norm_value : str or float, optional
        Specifies the Guide Star ID or the count rate/magnitude to which
        to normalize.
    norm_unit : str, optional
        Specifies the unit of norm_value ("FGS Magnitude", "FGS countrate",
        or "Guide Star ID")

    Returns
    -------
    gsid: str
        The guide star ID to be used to get the FGS count rate
    """
    if norm_unit.lower() == 'guide star id' and isinstance(norm_value, str):
        gsid = norm_value
    elif norm_unit.lower() == 'fgs magnitude' and isinstance(norm_value, (float,int)):
        gsid = query_jmag_mini_library(norm_value)
    else:
        LOGGER.error("Normalize Calculation: Value type does not match expectation for unit type.")
        raise TypeError("Mismatched normalization value ({}) and expected type ({}) for unit of {}".format(
            norm_value, ["float or int" if norm_unit.lower() == 'fgs magnitude' else 'str'][0], norm_unit))
    return gsid

def query_jmag_mini_library(fgsmag):
    """
    Since the count rate module requires a Guide Star ID (GSID), we have
    pre-defined GSIDs for Guide Stars with magnitudes between (approximately)
    FGS Magnitude of 10 and 14.
    For a FGS Mag of 10, 11, 12, 13, or 14, return the following GSID:
    10: N135000314 (FGS Mag (G1) = 10.09767, J Mag = 9.99699)
    11: N13A000006 (FGS Mag (G1) = 11.16569, J Mag = 11.00100)
    12: N13A000158 (FGS Mag (G1) = 12.16066, J Mag = 12.03899)
    13: N13A000125 (FGS Mag (G1) = 13.07528, J Mag = 13.00800)
    14: N13A002729 (FGS Mag (G1) = 13.95341, J Mag = 13.98499)

    Parameters
    ----------
    fgsmag: int
        Can only be int between (inclusive)10 and 14.

    Raises
    ------
    ValueError
        The input is not an integer between (inclusive) 10 and 14.
    """
    accepted_mags = np.arange(10, 15)
    if fgsmag not in accepted_mags:
        LOGGER.error("Normalize Calculation: FGS Magnitude not value. Please enter an integer between 10 and 14.")
        raise ValueError("Unacceptable FGS Magnitude value")
    else:
        LOGGER.info("Normalize Calculation: Using pre-defined Guide Star ID for FGS Magnitude of: {}".format(fgsmag))

    fgsmag_dict = {
        10: 'N135000314',
        11: 'N13A000006',
        12: 'N13A000158',
        13: 'N13A000125',
        14: 'N13A002729'
    }

    return fgsmag_dict[fgsmag]
