#!/usr/bin/env python
import numpy as np
import astropy.io.fits as pyfits
import argparse
import glob, os, sys


############################################################
# TO MAKE EXECUTABLE:
# 1. make sure the shebang is on top: #!/usr/bin/env python
# 2. chmod ugo+x fix_isim_orient.py
# 3. In .mysetenv or .bashrc, add script directory to PATH
#    e.g. setenv PATH ${PATH}:$HOME/Desktop/CV3-WFSC/
############################################################


def update_one_file(filename, outdir):
    # UPDATE SOME HEADER INFORMATION FOR ONE FILE AT A TIME:

    fits = pyfits.open(filename)
    print '\n',fits.info()
    outhdu = pyfits.HDUList()

    ####################################
    # GET IMAGE SIZE FROM SCI EXTENSION:
    ####################################
    xdim = fits[1].header["NAXIS1"]
    ydim = fits[1].header["NAXIS2"]
    det  = fits[0].header["DETECTOR"]
    print det,

    
    ####################################
    # DETERMINE IF IMAGE FLIP NECESSARY:
    ####################################
    flipud = False
    fliplr = False
    if det in ['NRCA2', 'NRCA4', 'NRCB3', 'NRCB1', 'NRCBLONG']:
        flipud = True
        print 'Data flipped up/down'
    elif det in ['NRCA1', 'NRCA3', 'NRCB2', 'NRCB4', 'NRCALONG', 'GUIDER2']:
        fliplr = True
        print 'Data flipped left/right'
    elif det in ['GUIDER1', 'NIRISS', 'NRS2']:
        flipud = True
        fliplr = True
        print 'Data flipped up/down, left/right'
    elif det in ['MIRIMAGE', 'NRS1']:
        print 'No need to flip data'
    else: print 'Detector {:s} not recognized; ignoring image flips'.format(det)


        
    ##################################
    # LOOP THROUGH ALL THE EXTENSIONS:
    ##################################
    iprint = 0
    for iext in xrange( len(fits) ):
        hdr = fits[iext].header
        data= fits[iext].data
        
        #---------------------
        # UPDATE HEADER FIRST:
        #---------------------
        if iext==1 and hdr["EXTNAME"]!="SCI": hdr.set("EXTNAME", "SCI")
        if iext==2 and hdr["EXTNAME"]!="ERR": hdr.set("EXTNAME", "ERR")
        if iext==3 and hdr["EXTNAME"]!="DQ" : hdr.set("EXTNAME", "DQ")
        
        if iext==0:
            
            # FIRST: UPDATE A FEW KEYWORDS FOR THE WAS:
            hdr.set("APERNAME", "APERNAME",'Added by {:s}'.format(os.path.basename(__file__)) )
            
            if hdr["INSTRUME"]=="MIRI": full_dim=1032
            else: full_dim=2048

                
            # FIX SUBARRAY keyword if not correct (CV-3 data looks alright; CV-2 was hardcoded to False):
            if xdim==full_dim and ydim==full_dim:
                hdr["SUBARRAY"]=False
            else:
                hdr["SUBARRAY"]=True


            # SECOND: ADD A FEW REQUIRED KEYWORDS FOR THE WAS:
            #if ("PUPIL" in hdr) and (hdr["PUPIL"]=="Weak Lens 1"): hdr.update("PUPIL","WLP8")
            
            hdr.set("SUBSTRT1", int(hdr["COLCORNR"]),'Added by {:s}'.format(os.path.basename(__file__)))
            hdr.set("SUBSTRT2", int(hdr["ROWCORNR"]),'Added by {:s}'.format(os.path.basename(__file__)))
            hdr.set("SUBSIZE1",ydim,'Added by {:s}'.format(os.path.basename(__file__)))
            hdr.set("SUBSIZE2",xdim,'Added by {:s}'.format(os.path.basename(__file__)))
            
            #------------------------------------------------
            # BELOW MIGHT BE RELIC OF EARLIER TEST CAMPAIGNS:
            #------------------------------------------------
            # Add FGS-specific keywords:
            if hdr["DETECTOR"]=="GUIDER1": hdr.set("FOCUSPOS", hdr["G1FOCPOS"],'Added by {:s}'.format(os.path.basename(__file__)))
            if hdr["DETECTOR"]=="GUIDER2": hdr.set("FOCUSPOS", hdr["G2FOCPOS"],'Added by {:s}'.format(os.path.basename(__file__)))     


            # Add NIRISS-specific keywords:
            if hdr["DETECTOR"]=="NIRISS":
                hdr.set("FOCUSPOS", hdr["FMCCUPOS"],'Added by {:s}'.format(os.path.basename(__file__)))


            # Add NIRSpec-specific keywords:
            if (hdr["DETECTOR"]=="NRS1") or (hdr["DETECTOR"]=="NRS2"):
                hdr.set("SLIT","SLIT",'Added by {:s}'.format(os.path.basename(__file__)))
                hdr.set("FCSRLPOS", hdr["RMA_STEP"],'Added by {:s}'.format(os.path.basename(__file__)))
                hdr.set("RMA_POS", hdr["RMA_MIC"],'Added by {:s}'.format(os.path.basename(__file__)))
                hdr.set("FOCUSPOS", hdr["CIMGFSOF"],'Added by {:s}'.format(os.path.basename(__file__)))

            hdr['HISTORY']=''
            hdr['HISTORY']='Header information and image orientation fixed by {:s}'.format(os.path.basename(__file__))


            
        #-------------------
        # UPDATE DATA:
        #-------------------        
        if data is not None:
            if flipud: data = np.flipud(data)
            if fliplr: data = np.fliplr(data)
            
            if hdr["EXTNAME"]: extname = hdr["EXTNAME"]             
            outhdu.append(pyfits.ImageHDU(data, header=hdr, name=extname))   
        else: outhdu.append(pyfits.ImageHDU(header=hdr))


        
    #-------------------------
    # WRITE THE FITS FILE OUT:
    #-------------------------    
    if outdir:
        if outdir[-1]!="/": outdir+="/"
        if not os.path.exists(outdir): os.makedirs(outdir)
        filename = outdir + filename
        
    outfile = filename.split('.fits')[0]+"_FIXED.fits"
    outhdu.writeto(outfile, clobber=True)  
    outhdu.close()
    fits.close()



    

def update_all(dir, string, outdir):
    # Grab the FITS files in directory and loop through each of them:

    os.chdir(dir)
    print os.getcwd()
    for filename in glob.glob('*'+string+'*cal.fits'):
        update_one_file(filename, outdir)


        


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Fix FITS file headers and image orientation for the WAS, e.g.: fix_isim_orient.py ~/Desktop/OTIS_Data -s NRCA -o ~/Desktop/OTIS_Data_Fixed")
    parser.add_argument("dir"               , metavar="indir" , help="directory to process FITS files", nargs=1)
    parser.add_argument("-o",  default='./' , metavar="outdir", help="output directory")
    parser.add_argument("-s",  default=''   , metavar="string", help="string to match anywhere in the filenames")

    args = parser.parse_args()

    dir=args.dir[0]
    string=args.s
    outdir = args.o

    update_all(dir, string, outdir)

