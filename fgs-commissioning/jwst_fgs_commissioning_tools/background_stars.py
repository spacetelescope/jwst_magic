# Add background stars to an FGS image by copying current image
import random
import numpy as np

from jwst_fgs_commissioning_tools.convert_image import counts_to_jmag
from jwst_fgs_commissioning_tools import log

def add_background_stars(image, stars, jmag, fgs_counts, guider):
    # Randomly create 5 locations on the image
    size = 2048
    nstars_random = 5

    # Determine jmag and fgs_counts of guide star
    if not fgs_counts:
        if not jmag:
            log.warning('No counts or J magnitude given, setting to default')
            jmag = 11
        fgs_counts = counts_to_jmag.jmag_to_fgs_counts(jmag, guider)
    else:
        jmag = counts_to_jmag.fgs_counts_to_jmag(fgs_counts, guider)

    # If the flag is simply set to "True", randomly place 5 stars on the image
    if stars == True:
        x_back = random.sample(range(size), nstars_random)
        y_back = random.sample(range(size), nstars_random)
        # Create the new stars 5 mags or more dimmer
        jmags_back = random.sample(set(np.linspace(jmag + 7, jmag + 4, 100)), nstars_random)

    # If users passed a dictionary to the bkgd_stars argument, add stars
    # according the dictionary
    elif type(stars) == dict:
        input_lengths = [len(stars[key]) for key in stars.keys()]
        if len(set(input_lengths)) != 1:
            raise ValueError('Invalid dictionary provided for background star '
                'positions and magnitudes. Ensure the same number of entries is '
                'provided for all fields.')

        x_back = stars['x']
        y_back = stars['y']
        jmags_back = stars['jmag']

    else:
        raise TypeError('Unfamiliar value passed to bkgd_stars: {} Please pass boolean or dictionary of background star x, y, jmag.'.format(stars))

    # Add stars to image
    add_data = np.copy(image)
    for x, y, jmag in zip(x_back, y_back, jmags_back):
        star_fgs_counts = counts_to_jmag.jmag_to_fgs_counts(jmag, guider)
        scale_factor = star_fgs_counts / fgs_counts

        star_data = image * scale_factor
        psfx = psfy = 2048

        x1 = max(0, int(x) - int(psfx / 2))
        x2 = min(2048, int(x) + int(psfx / 2) + 1)
        y1 = max(0, int(y) - int(psfy / 2))
        y2 = min(2048, int(y) + int(psfy / 2) + 1)

        if x > 1024:
            star_data = star_data[:x2 - x1]
        else:
            star_data = star_data[2048 - (x2 - x1):]
        if y > 1024:
            star_data = star_data[:, :y2 - y1]
        else:
            star_data = star_data[:, 2048 - (y2 - y1):]

        print('Adding background star with magnitude {:.1f} at location ({}, {}).'.format(jmag,
                                                                                          x, y))
        add_data[x1:x2, y1:y2] += star_data

    return add_data
