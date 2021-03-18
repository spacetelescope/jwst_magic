"""Add background stars to an FGS image by repeating the current image

Authors
-------
    - Lauren Chambers
    - Shannon Osborne

"""

# Standard Library Imports
import logging
import os
import random

# Third Party Imports
import fgscountrate
import numpy as np

# Local Imports
from jwst_magic.convert_image import renormalize
from jwst_magic.utils import utils

__location__ = os.path.realpath(os.path.join(os.getcwd(), os.path.dirname(__file__)))
PACKAGE_PATH = os.path.split(__location__)[0]
OUT_PATH = os.path.split(PACKAGE_PATH)[0]  # Location of out/ and logs/ directory

# Start logger
LOGGER = logging.getLogger(__name__)


def add_background_stars(image, stars, norm_value, norm_unit, guider,
                         save_file=False, root=None, out_dir=None):
    """Add artificial copies of the input PSF to the image to mimic
    background stars.

    Parameters
    ----------
    image : 2-D numpy array
        The input data with the original star
    stars : bool or dict
        Either a boolean or a dictionary that defines how to add the
        additional stars to the image
    norm_value : float
        Specifies the value to which to normalize.
    norm_unit : str
        Specifies the unit of norm_value (FGS Magnitude or FGS Counts)
    guider : int
        The number of the guider used to take the data
    save_file : bool, optional
        Save ASCII file containing y, x, and magnitude values for each
        background star. Default is False
    root: str, optional
        The root desired for output images if different than root in image
    out_dir : str, optional
        Where output FGS image(s) will be saved. If not provided, the
        image(s) will be saved to ../out/{root}.

    Returns
    -------
    add_data : 2-D numpy array
        Image array with original star and additional stars

    Raises
    ------
    TypeError
        Incompatible value passed to stars
    ValueError
        Invalid dictionary passed to stars
    """

    # Randomly create 5 locations on the image
    size = 2048
    nstars_random = 5

    # Determine fgs_mag and fgs_countrate of guide star
    fgs_countrate, fgs_mag = renormalize.convert_to_countrate_fgsmag(norm_value, norm_unit, guider)

    # If the flag is simply set to "True", randomly place 5 stars on the image
    if stars is True:
        x_back = random.sample(range(size), nstars_random)
        y_back = random.sample(range(size), nstars_random)
        # Create the new stars 5 mags or more dimmer
        fgs_mags_back = random.sample(set(np.linspace(fgs_mag + 7, fgs_mag + 4, 100)), nstars_random)

    # If users passed a dictionary to the bkgd_stars argument, add stars
    # according the dictionary
    elif type(stars) == dict:
        input_lengths = [len(stars[key]) for key in stars.keys()]
        if len(set(input_lengths)) != 1:
            raise ValueError('Invalid dictionary provided for background star '
                             'positions and magnitudes. Ensure the same number '
                             'of entries is provided for all fields.')

        x_back = stars['x']
        y_back = stars['y']
        fgs_mags_back = stars['fgs_mag']

    else:
        raise TypeError(
            'Unfamiliar value passed to bkgd_stars: {} Please pass boolean or dictionary of background '
            'star x, y, fgs_mag.'.format(stars)
        )

    # Add stars to image
    # Copy original data array
    add_data = np.copy(image)

    # (Try to) only use the data for added stars, not the noise
    mean = np.mean(image)
    image[image < mean] = 0
    for x, y, fgs_mag_back in zip(x_back, y_back, fgs_mags_back):
        if fgs_mag_back != 0:  # should have already removed all "bad" values marked with 0
            star_fgs_countrate = fgscountrate.convert_fgs_mag_to_cr(fgs_mag_back, guider)
            scale_factor = star_fgs_countrate / fgs_countrate

            star_data = image * scale_factor
            psfx = psfy = 2048

            # Crop the image to fit on the array at the specified location
            x1 = max(0, int(x) - int(psfx / 2))
            x2 = min(2048, int(x) + int(psfx / 2) + 1)
            y1 = max(0, int(y) - int(psfy / 2))
            y2 = min(2048, int(y) + int(psfy / 2) + 1)
            if y > 1024:
                star_data = star_data[:y2 - y1]
            else:
                star_data = star_data[2048 - (y2 - y1):]
            if x > 1024:
                star_data = star_data[:, :x2 - x1]
            else:
                star_data = star_data[:, 2048 - (x2 - x1):]

            LOGGER.info(
                'Background Stars: Adding background star with magnitude {:.1f} at location ({}, {}).'.
                format(fgs_mag_back, x, y))
            add_data[y1:y2, x1:x2] += star_data

    if save_file:
        # Set up out dir
        out_dir = utils.make_out_dir(out_dir, OUT_PATH, root)
        utils.ensure_dir_exists(out_dir)

        # Write catalog of all background PSFs
        background_psfs_file = os.path.join(out_dir, 'unshifted_background_psfs_{}_G{}.txt'.format(root, guider))
        all_cols = utils.create_cols_for_coords_counts(x_back, y_back, fgs_mags_back, val=None,
                                                       labels=None, inds=range(len(x_back)))
        utils.write_cols_to_file(background_psfs_file,
                                 labels=['y', 'x', 'fgs_mag'],
                                 cols=all_cols, log=LOGGER)

    return add_data
