# STDLIB
import os
import shutil

# THIRD PARTY
import matplotlib
if matplotlib.get_backend() != 'Qt5Agg':
    matplotlib.use('Qt5Agg')  # Make sure that we are using Qt5
print('Using backend: ', matplotlib.get_backend())
import numpy as np
from astropy.io import fits

# LOCAL
from . import log, utils, background_stars
from .convert_image import convert_image_to_raw_fgs
from .star_selector import select_psfs
from .fsw_file_writer import buildfgssteps

# Because Jupyter Notebook cannot open a matplotlib object, I have copied what is
# done in Run FGS Commissioning Tool.ipynb into this script that should be run in
# IPython
PACKAGE_PATH = os.path.dirname(os.path.realpath(__file__))
OUT_PATH = os.path.split(PACKAGE_PATH)[0]  # Location of out/ and logs/ directory
TASKNAME = 'run_all'

def run_all(image, guider, root=None, fgs_counts=None, jmag=None,
            nircam_det=None, nircam=True, global_alignment=False, steps=None,
            in_file=None, bkgd_stars=False, out_dir=None):

    # Determine filename root
    root = utils.make_root(root, image)

    # Set up logging
    taskname = '_'.join([TASKNAME, root])
    log_path = os.path.join(OUT_PATH, 'logs')
    logname = utils.get_logname(log_path, taskname)

    @log.logtofile(logname)
    def run_all_with_logging(image, guider, root=None, fgs_counts=None, jmag=None,
                             nircam_det=None, nircam=True, global_alignment=False,
                             steps=None, in_file=None, bkgd_stars=False,
                             out_dir=None):
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

        # Determine output directory
        out_dir_root = utils.make_out_dir(out_dir, OUT_PATH, root)

        log.info("Package directory: {}".format(PACKAGE_PATH))
        log.info("Processing request for {}. \nAll data will be saved in: {}".format(root, out_dir_root))
        log.info("Input image: {}".format(os.path.abspath(image)))
        utils.ensure_dir_exists(out_dir_root)

        # Either convert provided NIRCam image to an FGS image...
        fgs_im = convert_image_to_raw_fgs.convert_im(image, guider, root,
                                                     nircam=nircam,
                                                     fgs_counts=fgs_counts,
                                                     jmag=jmag,
                                                     nircam_det=nircam_det,
                                                     out_dir=out_dir)
        log.info("*** Image Conversion COMPLETE ***")

        if steps is None:
            steps = ['ID', 'ACQ1', 'ACQ2', 'TRK', 'LOSTRK']

        if bkgd_stars:
            fgs_im = background_stars.add_background_stars(fgs_im, jmag, fgs_counts, guider)

        # create reg file
        select_psfs.create_reg_file(fgs_im, root, guider,
                                    return_nref=False,
                                    global_alignment=global_alignment,
                                    in_file=in_file, out_dir=out_dir)
        log.info("*** Star Selection: COMPLETE ***")

        # create all files for FSW/DHAS/FGSES/etc.
        for step in steps:
            buildfgssteps.BuildFGSSteps(fgs_im, guider, root, step, out_dir=out_dir)
        log.info("*** FSW File Writing: COMPLETE ***")

    run_all_with_logging(image, guider, root=root, fgs_counts=fgs_counts,
                         jmag=jmag, nircam_det=nircam_det, nircam=nircam,
                         global_alignment=global_alignment,
                         steps=steps, in_file=in_file, bkgd_stars=bkgd_stars,
                         out_dir=out_dir)
