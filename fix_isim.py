#!/usr/bin/env python
#import pyfits
import astropy.io.fits as pyfits
import argparse
import glob, os, sys


def update_one_file(filename, outdir):
    # UPDATE SOME HEADER INFORMATION FOR ONE FILE AT A TIME:

    print ("Updating file "+filename+"...")
    fits = pyfits.open(filename, mode='update')

    print (fits.info())

    # FIRST: UPDATE A FEW KEYWORDS FOR THE WAS:
    fits[0].header.set("APERNAME", "APERNAME")

    if fits[1].header["EXTNAME"]!="SCI": fits[1].header.set("EXTNAME", "SCI")

    # FIX SUBARRAY keyword if not correct (CV-3 data looks alright; CV-2 was hardcoded to False):
    if fits[0].header["INSTRUME"]=="MIRI": full_dim=1032
    else: full_dim=2048

    if fits['SCI'].data.shape[0]==full_dim and fits['SCI'].data.shape[1]==full_dim:
        if fits[0].header["SUBARRAY"] != False:
            fits['SCI'].header.update("SUBARRAY", False)
            print ("  ---> fixed subarray keyword")
    else:
         if fits[0].header["SUBARRAY"] != True:
            fits[0].header.update("SUBARRAY", True)
            print ("  ---> fixed subarray keyword")



    # SECOND: ADD A FEW REQUIRED KEYWORDS FOR THE WAS:
    if ("PUPIL" in fits[0].header) and (fits[0].header["PUPIL"]=="Weak Lens 1"): fits[0].header.update("PUPIL","WLP8")

    fits[0].header.set("SUBSTRT1",int(fits[0].header["COLCORNR"]))
    fits[0].header.set("SUBSTRT2",int(fits[0].header["ROWCORNR"]))
    fits[0].header.set("SUBSIZE1",fits['SCI'].data.shape[0])
    fits[0].header.set("SUBSIZE2",fits['SCI'].data.shape[1])


    # Add FGS-specific keywords:
    if fits[0].header["DETECTOR"]=="GUIDER1": fits[0].header.set("FOCUSPOS", fits[0].header["G1FOCPOS"])
    if fits[0].header["DETECTOR"]=="GUIDER2": fits[0].header.set("FOCUSPOS", fits[0].header["G2FOCPOS"])


    # Add NIRISS-specific keywords:
    if fits[0].header["DETECTOR"]=="NIRISS":
        fits[0].header.set("FOCUSPOS", fits[0].header["FMCCUPOS"])


    # Add NIRSpec-specific keywords:
    if (fits[0].header["DETECTOR"]=="NRS1") or (fits[0].header["DETECTOR"]=="NRS2"):
        fits[0].header.set("SLIT","SLIT")
        fits[0].header.set("FCSRLPOS", fits[0].header["RMA_STEP"])
        fits[0].header.set("RMA_POS", fits[0].header["RMA_MIC"])
        fits[0].header.set("FOCUSPOS", fits[0].header["CIMGFSOF"])


    # SAVE CHANGES TO FILE:
    fits.flush()
    fits.close()
    print ("OK")




def update_all(dir, string, outdir):
    # Grab the FITS files in directory and loop through each of them:

    os.chdir(dir)
    print (os.getcwd()+"/"+dir)
    for filename in glob.glob('*'+string+'*.fits'):
        update_one_file(filename, outdir)



if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Fix FITS file headers for the WAS")
    parser.add_argument("dir", metavar="dir", nargs=1, help="directory to process FITS files")
    parser.add_argument("-o", help="output directory")
    parser.add_argument("-s",  default='', metavar="string"  , help="string to match anywhere in the filenames")
    args = parser.parse_args()

    dir=args.dir[0]
    string=args.s
    outdir = args.o

    update_all(dir, string, outdir)
