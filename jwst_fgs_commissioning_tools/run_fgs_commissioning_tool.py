# STDLIB
import os
import shutil

# THIRD PARTY
import matplotlib
if matplotlib.get_backend() != 'Qt5Agg':
    matplotlib.use('Qt5Agg')  # Make sure that we are using Qt5
import numpy as np
from astropy.io import fits
import pprint

# LOCAL
import select_psfs
import FGS_commissioning
import log
import nircam_to_fgs
import utils
import background_stars
import counts_to_jmag


# Because Jupyter Notebook cannot open a matplotlib object, I have copied what is
# done in Run FGS Commissioning Tool.ipynb into this script that should be run in
# IPython
LOCAL_PATH = os.path.dirname(os.path.realpath(__file__))
TASKNAME = 'run_all'

def run_all(image, guider, root=None, fgs_counts=None, jmag=None,
            nircam_det=None, nircam=True, global_alignment=False, in_file=None,
            bkgd_stars=False):

    taskname = '_'.join([TASKNAME, root])
    LOGNAME = utils.get_logname(os.path.join(LOCAL_PATH, 'logs'), taskname)
    print(LOGNAME)

    @log.logtofile(LOGNAME)
    def run_all_with_logging(image, guider, root=None, fgs_counts=None, jmag=None,
                             nircam_det=None, nircam=True, global_alignment=False,
                             in_file=None, bkgd_stars=False):
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
        in_file: str
            If this image comes with an incat or reg file, the file path
        """

        # print('log: ', pprint.pprint(dir(log)))

        if root is None:
            root = os.path.basename(image).split('.')[0]

        out_dir = os.path.join(LOCAL_PATH, 'out', root)

        log.info("Processing request for {}. \nAll data will be saved in: {}".format(root, out_dir))
        log.info("Input image: {}".format(os.path.abspath(image)))
        utils.ensure_dir_exists(out_dir)

        # Either convert provided NIRCam image to an FGS image...
        if nircam:
            log.info("This is a NIRCam image")
            fgs_im = nircam_to_fgs.convert_im(image, guider, fgs_counts=fgs_counts,
                                              jmag=jmag, nircam_det=nircam_det,
                                              return_im=True)

        # ... or process provided FGS image
        else:
            log.info("This is a FGS image")
            fgs_im = fits.getdata(image)

            # If J magnitude is provided, normalize the entire image to match that jmag
            if jmag:
                log.info("Normalizing to jmag = {}".format(jmag))
                fgs_counts = counts_to_jmag.jmag_to_fgs_counts(jmag, guider)
                fgs_im = fgs_im / np.sum(fgs_im) * fgs_counts

            # Correct high or low pixels
            fgs_im[fgs_im < 0] = 0.
            fgs_im[fgs_im > 65000] = 0.

            utils.ensure_dir_exists(os.path.join(out_dir, 'FGS_imgs'))
            shutil.copyfile(image, os.path.join(LOCAL_PATH, 'out', root, 'FGS_imgs',
                                                '{}.fits'.format(root)))
        if bkgd_stars:
            fgs_im = background_stars.add_background_stars(fgs_im, jmag, fgs_counts, guider)

        # create reg file
        nref = select_psfs.create_reg_file(fgs_im, root, guider, out_dir=out_dir,
                                           return_nref=True,
                                           global_alignment=global_alignment,
                                           in_file=in_file)

        # create all files for FSW/DHAS/FGSES/etc.
        FGS_commissioning.run_id(fgs_im, guider, root, out_dir=out_dir)
        FGS_commissioning.run_acq(fgs_im, guider, root, out_dir=out_dir)
        FGS_commissioning.create_lostrk(fgs_im, guider, root, nx=43, ny=43,
                                        out_dir=out_dir)
        FGS_commissioning.run_trk(fgs_im, guider, root, 5000, out_dir=out_dir, jitter=False)

    return run_all_with_logging(image, guider, root=root, fgs_counts=fgs_counts,
                                jmag=jmag, nircam_det=nircam_det, nircam=nircam,
                                global_alignment=global_alignment,
                                in_file=in_file, bkgd_stars=bkgd_stars)
