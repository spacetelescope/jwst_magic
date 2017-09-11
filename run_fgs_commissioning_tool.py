import os
import shutil

# From local
import FGS_commissioning
import select_psfs
import nircam_to_fgs
import utils
import log
import time

#Because Jupyter Notebook cannot open a matplotlib object, I have copied what is
# done in Run RGS Commissioning Tool.ipynb into this script that should be run in
# IPython
LOCAL_PAATH= os.path.dirname(os.path.realpath(__file__))
TASKNAME = 'run_all'
LOGNAME = utils.get_logname(os.path.join(local_path,'logs'), TASKNAME)

@log.logtofile(LOGNAME)
def run_all(im, guider, root=None, fgs_counts=None, jmag=None,
            nircam_mod=None, nircam=True,num_psfs=18,compact=False):
    if root is None:
        root = os.path.basename(im).split('.')[0]

    output_path = os.path.join(LOCAL_PATH,'out',root)
    utils.ensure_dir_exists(output_path)

    # convert NIRCam image to an FGS image
    if nircam:
        log.info("This is a NIRCam image")
        fgs_im = nircam_to_fgs.convert_im(im, guider, fgs_counts=fgs_counts,
                                      jmag=jmag, nircam_mod=None, return_im=True)
    else:
        log.info("This is a FGS image")
        fgs_im = utils.read_fits(im)[1]
        utils.ensure_dir_exists(os.path.join(output_path,'FGS_imgs'))
        if not os.path.exists(os.path.join(output_path,'FGS_imgs')):
            print(os.path.join(output_path,'FGS_imgs'))
            print('didnt work')

        shutil.copyfile(im,os.path.join(local_path,'out',root,'FGS_imgs','{}.fits'.format(root)))

    exit
    # create reg file
    utils.ensure_dir_exists(output_path)
    nref = select_psfs.create_reg_file(fgs_im,root,guider,output_path=output_path,
                                       return_nref=True, num_psfs=num_psfs,compact=compact)
    # create all files for FSW/DHAS/FGSES/etc.
    FGS_commissioning.run_ID(fgs_im, guider, root, nref=nref)
    FGS_commissioning.run_ACQ(fgs_im, guider, root)
    FGS_commissioning.create_LOSTRK(fgs_im, guider, root, nx=43, ny=43)
