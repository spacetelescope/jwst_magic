"""Run JWST MaGIC end-to-end

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
            nircam=True, nircam_det=None, norm_value=None, norm_unit=None,
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
        ``norm_value`` and ``norm_unit`` - used to normalize the input
            NIRCam image, to a desired number of FGS counts.
        ``out_dir`` - where output FGS image(s) will be saved. If not
            provided, the image(s) will be saved to ../out/{root}.
"""

# Standard Library Imports
import logging
import os
import shutil

# Third Party Imports
import matplotlib
if matplotlib.get_backend() != 'Qt5Agg':
    matplotlib.use('Qt5Agg')  # Make sure that we are using Qt5
print('Using backend: ', matplotlib.get_backend())
import numpy as np

# Local Imports
from jwst_magic.convert_image import background_stars, convert_image_to_raw_fgs
from jwst_magic.fsw_file_writer import buildfgssteps, write_files
from jwst_magic.star_selector import select_psfs
from jwst_magic.utils import utils

# Define paths
PACKAGE_PATH = os.path.dirname(os.path.realpath(__file__))
OUT_PATH = os.path.split(PACKAGE_PATH)[0]  # Location of out/ and logs/ directory

# Start logger
LOGGER = logging.getLogger(__name__)


def run_all(image, guider, root=None, norm_value=None, norm_unit=None,
            nircam_det=None, nircam=True, smoothing='default', steps=None,
            guiding_selections_file=None, bkgd_stars=False, out_dir=None, convert_im=True,
            star_selection=True, star_selection_gui=True, file_writer=True,
            masterGUIapp=None, copy_original=True, normalize=True,
            coarse_pointing=False, jitter_rate_arcsec=None, itm=False,
            shift_id_attitude=True, crowded_field=False, threshold=0.6):
    """
    This function will take any FGS or NIRCam image and create the outputs needed
    to run the image through the DHAS or other FGS FSW simulator. If no incat or
    reg file are provided, the user will be prompted to use a GUI for selecting
    the guide star and reference stars necessary for the FSW.

    Parameters
    ----------
    image: str
        The path to the image.
    guider: int
        Which guider is being used: 1 or 2
    root: str, optional
        The root desired for output images if different than root in image
    norm_value: float, optional
        The value to be used for normalization (depends on norm_unit)
    norm_unit: str, optional
        The unit to be used for normalization (expecting "FGS Counts" or "FGS Magnitude")
    nircam_det: str, optional
        The NIRCam detector used for this observation. Only applicable for NIRCam
        images and if cannot be parsed from header.
    nircam: bool, optional
        If this is a FGS image, set this flag to False
    smoothing: str, optional
        Options are "low" for minimal smoothing (e.g. MIMF), "high" for large
        smoothing (e.g. GA), or "default" for medium smoothing for other cases
    steps: list of strings, optional
        List of the steps to be completed
    guiding_selections_file: list of str, optional
        If this image comes with an incat or reg file, the file path
    bkgd_stars : boolean, optional
        Add background stars to the image?
    out_dir : str, optional
        Where output FGS image(s) will be saved. If not provided, the
        image(s) will be saved to ../out/{root}.
    convert_im : boolean, optional
        Run the convert_image module?
    star_selection : boolean, optional
        Run the  star_selector module?
    star_selection_gui : boolean, optional
        Show the GUI for the star_selector module?
    file_writer : boolean, optional
        Run the fsw_file_writer module?
    masterGUIapp : PyQt5.QtCore.QCoreApplication instance, optional
        The QApplication instance of the master GUI, if it is already
        open.
    copy_original : boolean, optional
        Copy the original data to {out_dir}/{root}?
    normalize : boolean, optional
        Normalize the provided image during convert_image?
    coarse_pointing : boolean, optional
        Apply jitter to simulate coarse pointing?
    jitter_rate_arcsec : float, optional
        If coarse_pointing is true, the rate of jitter in arcseconds
        per second to apply in the form of a Gaussian filter.
    itm : bool, Optional
        If this image come from the ITM simulator (important for normalization).
    threshold : float
        The threshold (aka count rate uncertainty factor) to use when writing
        FSW files.
    """

    # Determine filename root
    root = utils.make_root(root, image)

    # Set up logging
    utils.create_logger_from_yaml(__name__, root=root, level='DEBUG')

    # Determine output directory
    out_dir_root = utils.make_out_dir(out_dir, OUT_PATH, root)
    utils.ensure_dir_exists(out_dir_root)

    LOGGER.info("Package directory: {}".format(PACKAGE_PATH))
    LOGGER.info("Processing request for {}.".format(root))
    LOGGER.info("All data will be saved in: {}".format(out_dir_root))
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
                                                     nircam_det=nircam_det,
                                                     normalize=normalize,
                                                     norm_value=norm_value,
                                                     norm_unit=norm_unit,
                                                     coarse_pointing=coarse_pointing,
                                                     jitter_rate_arcsec=jitter_rate_arcsec,
                                                     logger_passed=True,
                                                     itm=itm)

        if bkgd_stars:
            if not normalize and not itm:
                norm_value = np.sum(fgs_im[fgs_im > np.median(fgs_im)])
                norm_unit = "FGS Counts"
            fgs_im = background_stars.add_background_stars(fgs_im, bkgd_stars,
                                                           norm_value, norm_unit,
                                                           guider)

        # Write converted image
        convert_image_to_raw_fgs.write_fgs_im(fgs_im, out_dir, root, guider)
        LOGGER.info("*** Image Conversion COMPLETE ***")
    # Or, if an FGS image was provided, use it!
    else:
        fgs_im = image
        LOGGER.info("Assuming that the input image is a raw FGS image")

    # Select guide & reference PSFs
    if star_selection:
        guiding_selections_path_list, all_found_psfs, psf_center_file = select_psfs.select_psfs(
            fgs_im, root, guider,
            smoothing=smoothing,
            guiding_selections_file=guiding_selections_file,
            out_dir=out_dir,
            logger_passed=True, masterGUIapp=masterGUIapp)
        LOGGER.info("*** Star Selection: COMPLETE ***")

    # Create all files for FSW/DHAS/FGSES/etc.
    if file_writer:
        out_dir = utils.make_out_dir(out_dir, OUT_PATH, root)

        # Shift the image and write out new fgs_im, guiding_selections, all_found_psfs, and psf_center files
        for i, guiding_selections_file in enumerate(guiding_selections_path_list):

            # Change out_dir to write data to guiding_config_#/ sub-directory next to the selections file
            if 'guiding_config' in guiding_selections_file:
                out_dir_fsw = os.path.join(out_dir, 'guiding_config_{}'.format(
                    guiding_selections_file.split('guiding_config_')[1].split('/')[0]))
            else:
                out_dir_fsw = out_dir

            if shift_id_attitude:
                fgs_im_fsw, guiding_selections_file_fsw, psf_center_file_fsw = buildfgssteps.shift_to_id_attitude(
                    fgs_im, root, guider, out_dir_fsw, guiding_selections_file=guiding_selections_file,
                    all_found_psfs_file=all_found_psfs, psf_center_file=psf_center_file,
                    crowded_field=crowded_field, logger_passed=True)
            else:
                fgs_im_fsw = fgs_im
                guiding_selections_file_fsw = guiding_selections_file
                psf_center_file_fsw = psf_center_file

            if steps is None:
                steps = ['ID', 'ACQ1', 'ACQ2', 'LOSTRK']

            for step in steps:
                fgs_files_obj= buildfgssteps.BuildFGSSteps(
                    fgs_im_fsw, guider, root, step, out_dir=out_dir_fsw, threshold=threshold,
                    logger_passed=True, guiding_selections_file=guiding_selections_file_fsw,
                    psf_center_file=psf_center_file_fsw, shift_id_attitude=shift_id_attitude,
                )
                write_files.write_all(fgs_files_obj)
            LOGGER.info("*** Finished FSW File Writing for Selection #{} ***".format(i+1))

        LOGGER.info("*** FSW File Writing: COMPLETE ***")

    LOGGER.info("*** Run COMPLETE ***")
