#! /usr/bin/env python

#setenv WEBBPSF_PATH /itar/jwst/tel/perrin/webbpsf-data
#setenv PYSYN_CDBS /grp/hst/cdbs/
#ipython
import sys, getopt
import os
import webbpsf
import pylab as plt
import astropy.io
from  astropy.io import ascii
from  astropy.io import fits
import numpy as np

def mkimgs(root,nx=2048,ny=2048,pad=200,templatesdir=None):
    '''
    Given the root, open an existing <root>.incat file with x,y, count rate, and
    flags indicating if the PSFs are in the image and in the catalog
    (1 = yes, 0 = no).

    From this list, create place PSFs from a WebbPSF PSF library in a 2048 x 2048
    image, creating a simulated NIRCam image. Boom.
    '''
    #fakenorm1=1  # if fakenorm=1 then psfs norm'd to 1
    print(fakenorm1)
    xscale=0.069
    yscale=0.069

    iopd = np.random.randint(0,10)  # includes low, excludes high
    fgs = webbpsf.FGS()
    fgs.pupilopd = (fgs.opd_list[0],iopd)
    if(fakenorm1 == 1):
         iopd = 6


    # read input catalog [ADU/sec]
    if os.access(root+'.incat', os.F_OK):
         print("reading ",root+'.incat')
         indat = ascii.read(root+'.incat')
         y = indat['xreal']  # x & y axes are
         x = indat['yreal']  # swapped in python
         cps   = indat['countrate']  # ADU/sec
         inimg = indat['img']
         incat = indat['cat']
    else:
         print("Trouble reading incat")
         exit()

    # read WebbPSF template
    # old PSFs had 19 wavelengths
    #templatesdir = '/users/holfeltz/simdata/dev/srd/mar2015/webb_templates/'
    # new PSFs have 50 wavelengths
    if templatesdir is None:
        templatesdir = '/ifs/jwst/wit/fgs/comsims/PSFs/'
    psf_counts_tab = templatesdir+'webb_template_psfs.dat'
    if os.access(psf_counts_tab, os.F_OK):
         print("reading ",psf_counts_tab)
         psf_ct = ascii.read(psf_counts_tab)
    else:
         print("Trouble reading psf_counts_tab")
         exit()

    # define truth image [ADU/sec]
    # make it too big so it's easy to lay down PSFs
    # then trim later
    truth = np.zeros((nx+pad,ny+pad))
    x0 = int((nx+pad)/2)
    y0 = int((nx+pad)/2)

    xx=np.round(x-1+(pad)/2)
    yy=np.round(y-1+(pad)/2)
    xx.astype(int)
    yy.astype(int)

    # accumulate stars onto truth img
    eol="\n"
    fout1 = open(root+'.truecat',"w")
    fout1.write("     xreal      yreal       c3x3         cps       img   cat\n")

    nstar=np.size(x)
    #kpsf=[0]*len(x)
    for i in range(nstar):
         if(inimg[i] == 1):
              ipsf = np.random.randint(0,1000)
              #kpsf[i] = ipsf
              if(fakenorm1 == 0):
                   dum = 'M0V_OPD{}_{}.fits'.format(iopd,str(ipsf).zfill(3))
              if(fakenorm1 == 1):
                   dum = 'M0V_OPD6_316.fits'
                   ipsf = 316

              hdulist = fits.open(os.path.join(templatesdir,dum))
              print(psf_ct[np.where(psf_ct["PSFfile"] == dum)])
              #print(y[i], x[i],dum,' * ',cps[i])
              psf = hdulist[0].data
              hdulist.close()
              # these already come out of webbPSF normalized
              curctot = psf_ct[np.where(psf_ct["PSFfile"] == dum)]["total_adu"][0]
              cur3x3 = psf_ct[np.where(psf_ct["PSFfile"] == dum)]["pk3x3_adu"][0]
              curpk = psf_ct[np.where(psf_ct["PSFfile"] == dum)]["pk_adu"][0]
              print(cps[i],curctot,cps[i]*curctot)
              curctot = cps[i]*curctot
              cur3x3 = cps[i]*cur3x3
              print(psf.max(),cps[i]*curctot.max())
              # truth img in counts per sec

              truth[xx[i]-50:xx[i]+51,yy[i]-50:yy[i]+51] = truth[xx[i]-50:xx[i]+51,yy[i]-50:yy[i]+51]+cps[i]*psf
              fout1.write ("%10.3f %10.3f %12i %12i %5i %5i %s"% (y[i], x[i],cur3x3, curctot,inimg[i], incat[i], eol))

    # trim to 2048x2048
    dx=int(nx/2)
    dy=int(ny/2)
    truth = truth[x0-dx:x0+dx,y0-dy:y0+dy]
    out=root+'_truth.fits'
    fitsobj = fits.HDUList()
    hdu = fits.PrimaryHDU()
    hdu.data = truth
    fitsobj.append(hdu)
    fitsobj.writeto(out,clobber=True)
    fitsobj.close()

    # write output catalog
    #eol="\n"
    #fout1 = open(root+'.truecat',"w")
    #fout1.write("     xreal      yreal       c3x3       cps      img   cat\n")
    #for ii in range(len(x)):
    #     fout1.write ("%10.3f %10.3f %10.3f %10.3f %5i %5i %s"% (y[ii], x[ii],
    #          tmp_c3x3[kpsf[ii]]*cps[ii], tmp_ctot[kpsf[ii]]*cps[ii],inimg[ii], incat[ii], eol))
    fout1.close()


if __name__ == "__main__":
     len(sys.argv)
     if (len(sys.argv) == 2):
          fakenorm1 = 0
          mkimgs(sys.argv[1])
     if (len(sys.argv) == 3):
          if(sys.argv[2]) == "-norm1":
             fakenorm1 = 1
          mkimgs(sys.argv[1])
     if (len(sys.argv) == 1):
          print("Usage:  mk_imgs.py root")
