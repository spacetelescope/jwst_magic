** See Run FGS Commissioning Tool.ipynb for examples on how to run nircam_to_fgs
and FGS_commissioning **

nircam_to_fgs.py
---------------
nircam_to_fgs takes a (simulated) NIRCAM image and converts it to an FGS image.
It does this by correcting bad pixels, re-binning, and normalizing the image for
a specific magnitude. It is optimized currently for the Global Alignment phase
and will roughly find each PSF in the image. This algorithm is a work in progress.

Running nircam_to_fgs.py will create an ‘out’ directory in the directory where
you are running this script. Inside the 'out' directory will be directories for
every image that you run the script for with different forms of files with
coordinate/countrate information and a png with the bin_norm_img with the center
of each PSF marked in blue.

Note: The image that you will use for FGS_commissioning.py is
FGS_imgs/{root}_bin_norm_img. Guider is indicated in the file name.



FGS_commissioning.py
--------------------
The FGS_commissioning tool (name in progress) will take in an FGS image and create
the necessary files to run the image through ID, ACQ, and TRK.

To run FGS_commissioning.py, open IPython and do the following:
In[1]: import FGS_commissioning
In[2]: root = ‘image_root’
     ex:  'jw00000_100_001_01100_00000_NRCA3_img'
In[3]: guider = 1
     Or guider = 2
In[4]: im = ‘full/path/to/bin_norm_imgs/image.fits’
     ex: '/Users/kbrooks/Documents/tel/FGS/out/
          {}/bin_norm_imgs/
          {}_NRCA3_img_G{}_binned_pad_norm.fits'.format(root,root,guider)


You can specify out_dir and template_path, but it is not necessary to do so, as
long as you don’t mind the files going to the same ‘out' directory as where
FGS_bin_tool.py put things, and you have a ‘templates’ directory with all the
.prc templates in the same directory as your scripts.

Now, you can run ID, ACQ, and TRK with:

id0 = FGS_commissioning.run_ID(im, guider, root, interactive=True)
acq1, acq2 = FGS_commissioning.run_ACQ(im, guider, root, interactive=True)
trk = FGS_commissioning.run_TRK(im, guider, root, num_frames=5000, interactive=True)

With id0, acq1, acq2, and trk being objects that have all the information needed
to run these through DHAS as attributes. You don’t need these since these
functions will do all the things you will need for now, it will allow you to
play around with what you can do with that information. Note, that num_frames
for track, can be whatever you want.
