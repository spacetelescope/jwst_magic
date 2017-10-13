# STDLIB
import os
import shutil

# Third Party
from astropy.io import fits

# LOCAL
import select_psfs
import FGS_commissioning
import log
import nircam_to_fgs
import utils


# Because Jupyter Notebook cannot open a matplotlib object, I have copied what is
# done in Run RGS Commissioning Tool.ipynb into this script that should be run in
# IPython
LOCAL_PATH = os.path.dirname(os.path.realpath(__file__))
TASKNAME = 'run_all'
LOGNAME = utils.get_logname(os.path.join(LOCAL_PATH, 'logs'), TASKNAME)

@log.logtofile(LOGNAME)
def run_all(image, guider, root=None, fgs_counts=None, jmag=None,
            nircam_mod=None, nircam=True, global_alignment=False, in_file=None):
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
    nircam_mod: str
        The NIRCam module used for this observation. Only applicable for NIRCam
        images and if cannot be parsed from header.
    nircam: bool
        If this is a FGS image, set this flag to False
    global_alignment: bool
        If this is not a global_alignment image, set this flag to False
    in_file: str
        If this image comes with an incat or reg file, the file path
    """

    if root is None:
        root = os.path.basename(image).split('.')[0]

    out_dir = os.path.join(LOCAL_PATH, 'out', root)

    log.info("Processing request for {}. \nAll data will be saved in: {}".format(root, out_dir))
    utils.ensure_dir_exists(out_dir)

    # convert NIRCam image to an FGS image
    if nircam:
        log.info("This is a NIRCam image")
        fgs_im = nircam_to_fgs.convert_im(image, guider, fgs_counts=fgs_counts,
                                          jmag=jmag, nircam_mod=nircam_mod,
                                          return_im=True)
    else:
        log.info("This is a FGS image")
        fgs_im = fits.getdata(image)
        utils.ensure_dir_exists(os.path.join(out_dir, 'FGS_imgs'))
        shutil.copyfile(image, os.path.join(LOCAL_PATH, 'out', root, 'FGS_imgs',
                                            '{}.fits'.format(root)))

    # create all files for FSW/DHAS/FGSES/etc.
    FGS_commissioning.run_id(fgs_im, guider, root, out_dir=out_dir,
                             global_alignment=global_alignment)
    FGS_commissioning.run_acq(fgs_im, guider, root, out_dir=out_dir,
                              global_alignment=global_alignment)
    FGS_commissioning.create_lostrk(fgs_im, guider, root, nx=43, ny=43,
                                    out_dir=out_dir,
                                    global_alignment=global_alignment)
