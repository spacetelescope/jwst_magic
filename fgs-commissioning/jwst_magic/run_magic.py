'''Run JWST MaGIC end-to-end

~Description

Authors
-------
    - Keira Brooks
    - Lauren Chambers

Use
---
    This module can be executed in a Python shell as such:
    ::
        from jwst_magic.convert_image import convert_image_to_raw_fgs
        convert_image_to_raw_fgs.convert_im(input_im, guider, root,
            nircam=True, nircam_det=None, fgs_counts=None, jmag=None,
            out_dir=None):

    Required arguments:
        ``input_image`` - filepath for the input (NIRCam or FGS) image
        ``guider`` - number for guider 1 or guider 2
        ``root`` - will be used to create the output directory, ./out/{root}
    Optional arguments:
        ``nircam`` - denotes if the input_image is an FGS or NIRCam
            image. If True, the image will be converted to FGS format.
            Unless out_dir is specified, the FGS-formatted image will
            be saved to ../out/{root}/FGS_imgs/{root}_binned_pad_norm.fits
        ``nircam_det`` - used to specify the detector of a provided
            NIRCam image. If left blank, the detector will be extracted
            from the header of the NIRCam FITS file.
        ``fgs_counts`` and ``jmag`` - used to normalize the input
            NIRCam image, either to a desired J magnitude or to a
            desired number of FGS counts. The jmag parameter can also
            be used to normalize an FGS image.
        ``out_dir`` - where output FGS image(s) will be saved. If not
            provided, the image(s) will be saved to ../out/{root}.
'''

# STDLIB
import os
import shutil
import logging

# THIRD PARTY
import matplotlib
if matplotlib.get_backend() != 'Qt5Agg':
    matplotlib.use('Qt5Agg')  # Make sure that we are using Qt5
print('Using backend: ', matplotlib.get_backend())
import numpy as np
from astropy.io import fits

# LOCAL
from . import utils, background_stars
from .convert_image import convert_image_to_raw_fgs
from .star_selector import select_psfs
from .fsw_file_writer import buildfgssteps

# Define paths
PACKAGE_PATH = os.path.dirname(os.path.realpath(__file__))
OUT_PATH = os.path.split(PACKAGE_PATH)[0]  # Location of out/ and logs/ directory

# Start logger
LOGGER = logging.getLogger(__name__)


def run_all(image, guider, root=None, fgs_counts=None, jmag=None,
            nircam_det=None, nircam=True, global_alignment=False, steps=None,
            in_file=None, bkgd_stars=False, out_dir=None, convert_im=True,
            star_selection=True, star_selection_gui=True, file_writer=True,
            masterGUIapp=None, copy_original=True, normalize=True,
            coarse_pointing=False, jitter_rate_arcsec=None):
    """
    This function will take any FGS or NIRCam image and create the outputs needed
    to run the image through the DHAS or other FGS FSW simulator. If no incat or
    reg file are provided, the user will be prompted to use a GUI for selecting
    the guide star and reference stars necessary for the FSW.

    Parameters
    ==========
    image: str
        The path to the image.
    guider: int
        Which guider is being used: 1 or 2
    root: str
        The root desired for output images if different than root in image
    fgs_counts: float
        If the user wants to specify the FGS counts, they do so here
    jmag: float
        Like fgs_counts but for the J magnitude (NOT JAB)
    nircam_det: str
        The NIRCam detector used for this observation. Only applicable for NIRCam
        images and if cannot be parsed from header.
    nircam: bool
        If this is a FGS image, set this flag to False
    global_alignment: bool
        If this is not a global_alignment image, set this flag to False
    steps: list of strings
        List of the steps to be completed
    in_file: str
        If this image comes with an incat or reg file, the file path
    """

    # Determine filename root
    root = utils.make_root(root, image)

    # Set up logging
    utils.create_logger_from_yaml(__name__, root=root, level='DEBUG')

    # Determine output directory
    out_dir_root = utils.make_out_dir(out_dir, OUT_PATH, root)
    utils.ensure_dir_exists(out_dir_root)

    LOGGER.info("Package directory: {}".format(PACKAGE_PATH))
    LOGGER.info("Processing request for {}. \nAll data will be saved in: {}".format(root,
                                                                                    out_dir_root))
    LOGGER.info("Input image: {}".format(os.path.abspath(image)))

    # Copy input image into out directory
    if copy_original:
        try:
            shutil.copy(os.path.abspath(image), out_dir_root)
        except shutil.SameFileError:
            pass

    # Either convert provided NIRCam image to an FGS image...
    if convert_im:
        fgs_im = convert_image_to_raw_fgs.convert_im(image, guider, root,
                                                     nircam=nircam,
                                                     fgs_counts=fgs_counts,
                                                     jmag=jmag,
                                                     nircam_det=nircam_det,
                                                     out_dir=out_dir,
                                                     logger_passed=True,
                                                     normalize=normalize,
                                                     coarse_pointing=coarse_pointing,
                                                     jitter_rate_arcsec=jitter_rate_arcsec)
        if bkgd_stars:
            fgs_im = background_stars.add_background_stars(fgs_im, bkgd_stars, jmag, fgs_counts, guider)

        # Write converted image
        convert_image_to_raw_fgs.write_FGS_im(fgs_im, out_dir, root, guider)
        LOGGER.info("*** Image Conversion COMPLETE ***")
    # Or, if an FGS image was provided, use it!
    else:
        fgs_im = image
        LOGGER.info("Assuming that the input image is a raw FGS image")

    # create reg file
    if star_selection:
        if star_selection_gui:
            select_psfs.create_reg_file(fgs_im, root, guider,
                                        return_nref=False,
                                        global_alignment=global_alignment,
                                        in_file=in_file, out_dir=out_dir,
                                        logger_passed=True, masterGUIapp=masterGUIapp)
        else:
            LOGGER.info("Star Selection: Reading guide and reference star positions from {}".format(in_file))
        LOGGER.info("*** Star Selection: COMPLETE ***")

    # create all files for FSW/DHAS/FGSES/etc.
    if file_writer:
        if steps is None:
            steps = ['ID', 'ACQ1', 'ACQ2', 'LOSTRK']
        for step in steps:
            buildfgssteps.BuildFGSSteps(fgs_im, guider, root, step,
                                        out_dir=out_dir, logger_passed=True)
        LOGGER.info("*** FSW File Writing: COMPLETE ***")