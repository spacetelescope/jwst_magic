# STDLIB
import os
import shutil

# THIRD PARTY
import numpy as np
import matplotlib.pyplot as plt

# LOCAL
import select_psfs
import FGS_commissioning
import log
import nircam_to_fgs
import utils
import getbias


# Because Jupyter Notebook cannot open a matplotlib object, I have copied what is
# done in Run RGS Commissioning Tool.ipynb into this script that should be run in
# IPython
LOCAL_PATH = os.path.dirname(os.path.realpath(__file__))
TASKNAME = 'run_all'
LOGNAME = utils.get_logname(os.path.join(LOCAL_PATH, 'logs'), TASKNAME)

@log.logtofile(LOGNAME)
def run_all(im, guider, root=None, fgs_counts=None, jmag=None,
            nircam_mod=None, nircam=True, global_alignment=False, in_file=None):
    """
    This function will take any FGS or NIRCam image and create the outputs needed
    to run the image through the DHAS or other FGS FSW simulator. If no incat or
    reg file are provided, the user will be prompted to use a GUI for selecting
    the guide star and reference stars necessary for the FSW.

    Parameters
    ==========
    im: str
        The path to the image.
    guider: int
        Which guider is being used: 1 or 2
    root: str
        The root desired for output images if different than root in im
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
        root = os.path.basename(im).split('.')[0]

    out_dir = os.path.join(LOCAL_PATH, 'out', root)

    log.info("Processing request for {}. \nAll data will be saved in: {}".format(root, out_dir))
    utils.ensure_dir_exists(out_dir)

    # convert NIRCam image to an FGS image
    if nircam:
        log.info("This is a NIRCam image")
        fgs_im = nircam_to_fgs.convert_im(im, guider, fgs_counts=fgs_counts,
                                          jmag=jmag, nircam_mod=nircam_mod,
                                          return_im=True)
    else:
        log.info("This is a FGS image")
        fgs_im = utils.read_fits(im)[1]

        # # Add bias, pedestal, noise...
        # # (Assume ID stage)
        # nreads = nramps = 2
        # nx = ny = 2048
        # poissonfactor = 1
        # fgs_bias = getbias.getbias(guider, nreads, nramps, nx, ny, bzp=False)
        # time_normed_im = fgs_im * 0.338
        # fgs_bias += 0.25 * np.random.poisson(poissonfactor * time_normed_im)
        # # plt.imshow(fgs_bias[-1])  # Just take the first frame... probs wrong
        # # plt.show()
        # fgs_im += fgs_bias[0]

        utils.ensure_dir_exists(os.path.join(out_dir, 'FGS_imgs'))
        shutil.copyfile(im, os.path.join(LOCAL_PATH, 'out', root, 'FGS_imgs',
                                         '{}.fits'.format(root)))

    # create reg file
    nref = select_psfs.create_reg_file(fgs_im, root, guider, out_dir=out_dir,
                                       return_nref=True,
                                       global_alignment=global_alignment,
                                       in_file=in_file)

    # create all files for FSW/DHAS/FGSES/etc.
    FGS_commissioning.run_ID(fgs_im, guider, root, nref=nref, out_dir=out_dir)
    FGS_commissioning.run_ACQ(fgs_im, guider, root, out_dir=out_dir)
    FGS_commissioning.create_LOSTRK(fgs_im, guider, root, nx=43, ny=43, out_dir=out_dir)
