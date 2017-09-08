import sys
import numpy as np
from astropy.io import ascii as asc


# I think this is like a function definition? ~name, args~
# name of file must equal name of procedure
def use_imgs_Lauren(root, overlap, biasZeroPt=0, biasKTC=0, biasPed=0, poissonNoise=0):

   # Print user input
   print(root)
   print('          overlap =',overlap)
   print('       biasZeroPt =',biasZeroPt)
   print('          biasKTC =',biasKTC)
   print('          biasPed =',biasPed)
   print('     poissonNoise =',poissonNoise)

   yoffset=12

   for guider in [1,2]:
      gid='_G' + str(guider)

      # Define useful constants
      tcdsID=0.338    # neil 18 may 12
      tcdsTRK=0.0256  # neil 18 may 12
      tcdsFG=0.0512   # neil 18 may 12
      tcdsACQ1=0.17286 # 13 dec 13 vicki balzano
      tcdsACQ2=0.05016 # 20 apr 17 sth from TNO_CSA_50084_1113_RevP1 (oct 2015)
      #tcdsACQ1=2.*0.16384 # 13 dec 13 sth
      #tcdsACQ2=0.01024 # 13 dec 13 sth   

      # predicted countrate
      #readcol,root+'.incat',ctot,/silent,format='(x,x,f)'  # counts per tcdsID

      # observation parameters
      nx = 2048
      ny = 2048
      h = 64  
      nreads = 2
      nramps = 2
      nz = nreads * nramps
      #incat=1

      # Read in .incat file
      incat = asc.read(root + '.incat', names=['x', 'y', 'ctot', 'inimg', 'incat'])
      incat['x'] = incat['x'].astype(float)
      incat['y'] = incat['y'].astype(float)

      #c3x3=ctot*0.66
      countrate = incat['ctot'] #*tcdsID # counts/sec NOT ADU
      norm = np.sum(countrate)
      #y=y-12.0 # shift stars y-12 for ernie 

      # Other parameters...
      nstar = len(incat['x'])
      gs=0  # index of guide star
      xgs = incat['x'][gs]
      ygs = incat['y'][gs]
      #brightgs=bright[gs]

      if_overlap_then_nstrips = [[0,  32],
                                 [2,  33],
                                 [4,  34],
                                 [6,  35],
                                 [8,  36],
                                 [10, 37],
                                 [12, 39],
                                 [14, 40],
                                 [16, 42]]
      for o, n in if_overlap_then_nstrips:
         if overlap == o:
            nstrips = n      

      # define truth image [ADU/sec]
      print('reading truth file')
      #truth=float(mrdfits(root+'_truth.fits',/unsigned,/silent))
      truth=float(mrdfits(root+'.fits',/silent))
      #truth=float(mrdfits(root+'.fits',1,/silent))  # ITM img data in ext 1
      truth=norm*truth/total(truth) # make sure it's normalized

      sky=replicate(0.,nx,ny)
      sky[*,*]=truth[*,*]*tcdsID      # ADU

      print('baseline',max(sky))
      writefits,root+gid+'_IDsky.fits',sky

      print('truth',max(truth))
      #writefits,root+gid+'_truth.fits',truth  # ADU/sec

      # star catalog in real pixs
      get_lun,outlun
      print('writing gssscat')
      openw,outlun,root+gid+'_ID.gssscat'
      for i=0,n_elements(x)-1 do printf,outlun,x[i],y[i]
      close,outlun

      # write stc files using offset, rotated catalog
      # printing 0.66*countrate 24 apr 17 sth (and smaller thresholds)
      print('writing ID stc file')
      if(guider eq 1) then begin
         openw,outlun,root+'_G1_ID.stc'
         i=where(incat eq 1, n)
         file_delete,'ideal.tmp',/allow_nonexistent
         g1RPtoIA,x[i],y[i]
         readcol,'ideal.tmp',inum,xi,yi,format='(i,f,f)',/silent
         for j=0,n-1 do printf,outlun,inum[j],xi[j],yi[j],0.66*countrate[i[j]],format='(i,d,d,1x,e)'
         close,outlun
      endif
      if(guider eq 2) then begin
         openw,outlun,root+'_G2_ID.stc'
         i=where(incat eq 1, n)
         file_delete,'ideal.tmp',/allow_nonexistent
         g2RPtoIA,x[i],y[i]
         readcol,'ideal.tmp',inum,xi,yi,format='(i,f,f)'
         for j=0,n-1 do printf,outlun,inum[j],xi[j],yi[j],0.66*countrate[i[j]],format='(i,d,d,1x,e)'
         close,outlun
      endif

      nx=128
      ny=128
      if(guider eq 1) then begin
         file_delete,'ideal.tmp',/allow_nonexistent
         g1RPtoIA,x[0]-(nx/2),y[0]-(nx/2)
         readcol,'ideal.tmp',inum,xi,yi,format='(i,f,f)',/silent
         openw,outlun,root+'_G1_ACQ1.stc'
         printf,outlun,inum[0],xi[0],yi[0],0.66*countrate[0],format='(i,d,d,1x,e)'
         close,outlun
      endif
      if(guider eq 2) then begin
         file_delete,'ideal.tmp',/allow_nonexistent
         g2RPtoIA,x[0]-(nx/2),y[0]-(nx/2)
         readcol,'ideal.tmp',inum,xi,yi,format='(i,f,f)',/silent
         openw,outlun,root+'_G2_ACQ1.stc'
         printf,outlun,inum[0],xi[0],yi[0],0.66*countrate[0],format='(i,d,d,1x,e)'
         close,outlun
      endif

      nx=32
      ny=32
      if(guider eq 1) then begin
         file_delete,'ideal.tmp',/allow_nonexistent
         g1RPtoIA,x[0]-(nx/2),y[0]-(nx/2)
         readcol,'ideal.tmp',inum,xi,yi,format='(i,f,f)',/silent
         openw,outlun,root+'_G1_ACQ2.stc'
         printf,outlun,inum[0],xi[0],yi[0],0.66*countrate[0],format='(i,d,d,1x,e)'
         close,outlun
      endif
      if(guider eq 2) then begin
         file_delete,'ideal.tmp',/allow_nonexistent
         g2RPtoIA,x[0]-(nx/2),y[0]-(nx/2)
         readcol,'ideal.tmp',inum,xi,yi,format='(i,f,f)',/silent
         openw,outlun,root+'_G2_ACQ2.stc'
         printf,outlun,inum[0],xi[0],yi[0],0.66*countrate[0],format='(i,d,d,1x,e)'
         close,outlun
      endif
      free_lun,outlun
      nx=2048
      ny=2048

      # bias img
      if(keyword_set(biasZeroPt) eq 0) then bzp = 0
      if(keyword_set(biasKTC) eq 0) then bktc = 0
      if(keyword_set(biasPed) eq 0) then bp = 0
      if(keyword_set(biasZeroPt) eq 1) then bzp = 1
      if(keyword_set(biasKTC) eq 1) then bktc = 1
      if(keyword_set(biasPed) eq 1) then bp = 1

      IDbias = getIDbias(guider, bzp, bktc, bp)
      writefits,root+gid+'_IDbias.fits',IDbias

      # put it all together
      nz=nreads*nramps
      id=replicate(0.,nx,ny,nz)
      for iz=0,nz-1 do begin
         if(poissonNoise =1) then tmp=poidev(sky[*,*])
         if(poissonNoise =0) then tmp=sky[*,*]
         id[*,*,iz]=id[*,*,iz]+IDbias[*,*,iz]+tmp[*,*] # ADU
      endfor
      for iz=1,nz-1,nreads do begin
         if(poissonNoise =1) then tmp=poidev(sky[*,*])
         if(poissonNoise =0) then tmp=sky[*,*]
         id[*,*,iz]=id[*,*,iz]+tmp[*,*]
      endfor
      isat=where(id gt 65000, nsat)
      if(nsat gt 0) then id[isat] = 65000
      i=where(finite(id) eq 0, n)
      if(n gt 0) then id[i]=0.0
      writefits,root+gid+'_IDff.fits',id

      # full frame cds img (does not include drifts)
      cds=replicate(0.,nx,ny,nz/2)
      j=0
      for i=1,nz-1,2 do begin
         tmp=id[*,*,i]-id[*,*,i-1]
         cds[*,*,j]=tmp
         j=j+1
      endfor
      cds[isat]=25000.
      writefits,root+gid+'_IDcds.fits',cds

      #yoffset=12 # for ernie
      #yoffset=0
      x1=fltarr(nstrips*nz)           # strip coords in final full frame img
      x2=fltarr(nstrips*nz)
      y1=fltarr(nstrips*nz)
      y2=fltarr(nstrips*nz)
      k=0
      for i=0,nstrips-1 do begin
         for j=0,nz-1 do begin
            x1[k]=0
            x2[k]=nx-1
            y1[k]=i*(h-overlap)+yoffset
            y2[k]=y1[k]+h-1 
            k=k+1
         endfor
      endfor

      # extract strips from ff img
      #ii=ii-x0 & jj=jj-x0
      #kk=kk-y0 & ll=ll-y0
      strips=replicate(0.,nx,h,nstrips*nreads*nramps)
      nn=0
      for i=0,nstrips-1 do begin
         for iz=0,nz-1 do begin
            strips[*,*,nn]=id[x1[nn]:x2[nn],y1[nn]:y2[nn],iz]
            nn=nn+1
         endfor
      endfor
      ineg=where(strips lt 0, nneg) & if(nneg ge 1) then strips[ineg]=0
      isat=where(strips ge 65000, nsat) & if(nsat ge 1) then strips[isat]=65000
      i=where(finite(strips) eq 0, n)
      if(n gt 0) then strips[i]=0.0
      strips=uint(strips)
      print('reading magic hdr imgs')  
      if(guider eq 1) then dum=mrdfits('G1magicHdrImg.fits',0,hdr0,/unsigned,/silent)
      if(guider eq 2) then dum=mrdfits('G2magicHdrImg.fits',0,hdr0,/unsigned,/silent)
      writefits,root+gid+'_IDstrips.fits',strips,hdr0
      print('converting fits2dat')
      print('nz = ',nz)
      fits2dat=1                      # convert to ascii hex dat format for FGSES
      if(fits2dat eq 1) then begin
         print('converting ID strips from fits to dat')
         outfile=root+gid+'_IDstrips.dat'
         get_lun, outlun666 
         openw,outlun666,outfile,width=10240,/swap_if_little_endian
         nz=nstrips*nreads*nramps
         for iz=0,nz-1 do begin
            tmp=strips[*,*,iz]
            nn=n_elements(tmp)-1L
            for ipix=0L,nn do begin
               pix=strupcase(string(tmp[ipix],format='(z4.4)'))
               printf,outlun666,pix,format='($,a4," ")'
            endfor
         endfor
         fits2dat=0
         close,outlun666
         free_lun,outlun666
      endif


      # ACQ1 = repeat 6x(reset, drop, read, drop, read)
      # cds time = 2*(128*128*10.e-6) = 0.32768s
      nx=128 & ny=128                 # set img size
      nreads=2 & nramps=6 & nz=nreads*nramps
      acq=fltarr(nx,ny,nz)
      if(keyword_set(background) eq 1) then bkg=500+10.*randomn(seed,[nx,ny,nz],/normal) # add randomized background
      if(keyword_set(background) eq 0) then bkg=replicate(0.,nx,ny,nz)
      x1=fix(xgs)-nx/2 & x2=fix(xgs)+nx/2  # select appropriate region of 0th read bias
      y1=fix(ygs)-ny/2 & y2=fix(ygs)+ny/2
      if((x2-x1) eq 128) then x2=x2-1
      if((y2-y1) eq 128) then y2=y2-1

      # star catalog in real pixs
      get_lun,outlun
      openw,outlun,root+gid+'_ACQ1.cat'
      printf,outlun,xgs,ygs
      close,outlun
      free_lun,outlun

      # make bias frames
      ACQ1bias = getACQ1bias(guider, bzp, bktc, bp, x[0], y[0])
      writefits,root+gid+'_ACQ1bias.fits',ACQ1bias

      xx=nx/2
      yy=ny/2
      sky=replicate(0.,nx,ny,nz)
      if(keyword_set(SCdrift) eq 0) then begin
         for i=0,nz-1 do sky[*,*,i]=truth[x1:x2,y1:y2]*tcdsACQ1  # -> ADU
      print('ACQ1',max(sky))
      endif
      writefits,root+gid+'_ACQ1sky.fits',sky


      # ACQ1 = repeat 6x(reset, drop, read, drop, read)
      if(keyword_set(poissonNoise) eq 1) then begin
         for iz=0,nz-1 do acq[*,*,iz]=ACQ1bias[*,*,iz]+poidev(2.*sky[*,*,iz])
         for iz=1,nz-1,nreads do acq[*,*,iz]=acq[*,*,iz]+poidev(2.*sky[*,*,iz])
      endif
      if(keyword_set(poissonNoise) eq 0) then begin
         for iz=0,nz-1 do acq[*,*,iz]=ACQ1bias[*,*,iz]+2.*sky[*,*,iz]
         for iz=1,nz-1,nreads do acq[*,*,iz]=acq[*,*,iz]+2.*sky[*,*,iz]
      endif
      ineg=where(acq lt 0, nneg) & if(nneg ge 1) then acq[ineg]=0            # neg pixs -> 0
      isat=where(acq ge 65000, nsat) & if(nsat ge 1) then acq[isat]=65000    # sat'd pixs -> 65K
      writefits,root+gid+'_ACQ1.fits',uint(acq)
      #fits2dat,root+gid+'_ACQ1.fits','ACQ1' # 18oct12
      fits2dat=1                      # convert to ascii hex dat format for FGSES
      if(fits2dat eq 1) then begin
         print('converting ACQ1 from fits to dat')
         outfile=root+gid+'_ACQ1.dat'
         get_lun, outlun666  
         openw,outlun666,outfile,width=10240,/swap_if_little_endian
         for iz=0,nz-1 do begin
            tmp=acq[*,*,iz]
            nn=n_elements(tmp)-1L
            for i=0L,nn do begin      # convert to ascii hex dat format
               pix=strupcase(string(tmp[i],format='(z4.4)'))
               printf,outlun666,pix,format='($,a4," ")'
            endfor
         endfor
         close,outlun666
         free_lun,outlun666
         fits2dat=0
      endif

      cds=replicate(0.,nx,ny,nz/2)          # make cds
      j=0
      for i=0,nz-1,2 do begin
         tmp=acq[*,*,i+1]-acq[*,*,i]
         cds[*,*,j]=tmp
         j=j+1
      endfor
      writefits,root+gid+'_ACQ1cds.fits',cds

      # ACQ2 = repeat 5x(reset, drop, read, drop, drop, drop, read, drop)
      # cds time = 4*(32*32*10e-6) = 0.04096s
      nx=32 & ny=32                   # set img size
      nreads=2 & nramps=5 & nz=nreads*nramps
      acq=fltarr(nx,ny,nz)
      if(keyword_set(background) eq 1) then bkg=500+10.*randomn(seed,[nx,ny,nz],/normal) # add randomized background
      if(keyword_set(background) eq 0) then bkg=replicate(0.,nx,ny,nz)
      x1=fix(xgs)-nx/2 & x2=fix(xgs)+nx/2      # select appropriate region of 0th read bias
      y1=fix(ygs)-ny/2 & y2=fix(ygs)+ny/2
      if((x2-x1) eq 32) then x2=x2-1
      if((y2-y1) eq 32) then y2=y2-1

      # star catalog in real pixs
      get_lun,outlun
      openw,outlun,root+gid+'_ACQ2.cat'
      printf,outlun,xgs,ygs
      close,outlun
      free_lun,outlun

      # make bias frames
      ACQ2bias = getACQ2bias(guider, bzp, bktc, bp, x[0], y[0])
      writefits,root+gid+'_ACQ2bias.fits',ACQ2bias

      xx=nx/2
      yy=ny/2
      sky=replicate(0.,nx,ny,nz)
      if(keyword_set(SCdrift) eq 0) then begin
            for i=0,nz-1 do sky[*,*,i]=truth[x1:x2,y1:y2]*tcdsACQ2  # -> ADU
      print('ACQ2',max(sky)
      endif
      writefits,root+gid+'_ACQ2sky.fits',sky


      # ACQ2 = repeat 5x(reset, drop, read, drop, drop, drop, read, drop)
      # second read has 3x the signal of the 1st read
      if(keyword_set(poissonNoise) eq 1) then begin
         for iz=0,nz-1 do acq[*,*,iz]=ACQ2bias[*,*,iz]+poidev(2.*sky[*,*,iz])
         for iz=1,nz-1,nreads do acq[*,*,iz]=acq[*,*,iz]+poidev(4.*sky[*,*,iz]) # 13 dec 13 sth
      endif
      if(keyword_set(poissonNoise) eq 0) then begin
         for iz=0,nz-1 do acq[*,*,iz]=ACQ2bias[*,*,iz]+2.*sky[*,*,iz]
         for iz=1,nz-1,nreads do acq[*,*,iz]=acq[*,*,iz]+4.*sky[*,*,iz]
      endif
      ineg=where(acq lt 0, nneg) & if(nneg ge 1) then acq[ineg]=0         # neg pixs -> 0
      isat=where(acq ge 65000, nsat) & if(nsat ge 1) then acq[isat]=65000 # sat'd pixs -> 60K
      writefits,root+gid+'_ACQ2.fits',uint(acq)
      #fits2dat,root+gid+'_ACQ2.fits','ACQ2' # 18oct12
      fits2dat=1
      if(fits2dat eq 1) then begin
         print('converting ACQ2 from fits to dat')
         outfile=root+gid+'_ACQ2.dat'
         get_lun, outlun666
         openw,outlun666,outfile,width=10240,/swap_if_little_endian
         for iz=0,nz-1 do begin
            tmp=acq[*,*,iz]
            nn=n_elements(tmp)-1L
            for i=0L,nn do begin      # convert to ascii hex dat format
               pix=strupcase(string(tmp[i],format='(z4.4)'))
               printf,outlun666,pix,format='($,a4," ")'
            endfor
         endfor
         close,outlun666
         free_lun,outlun666
         fits2dat=0
      endif

      cds=replicate(0.,nx,ny,nz/2)       # make cds
      j=0
      for i=0,nz-1,2 do begin
         tmp=acq[*,*,i+1]-acq[*,*,i]
         cds[*,*,j]=tmp
         j=j+1
      endfor
      writefits,root+gid+'_ACQ2cds.fits',cds

      #if(guider eq 2) then begin
      # LOS PSF
         nx=43 & ny=43
         x1=fix(xgs)-nx/2 & x2=fix(xgs)+nx/2
         y1=fix(ygs)-ny/2 & y2=fix(ygs)+ny/2
         xx=fix(nx/2)-1
         yy=fix(ny/2)-1
         sky=truth[x1:x2,y1:y2]
         tmp=rebin(sky,258,258)
         tmp=tmp[1:255,1:255]
         tmp=tmp/total(tmp)

         print('converting LOS PSF from fits to dat')
         get_lun,outlun
         openw,outlun,root+'_LOSpsf.dat',width=4080,/swap_if_little_endian
         for i=0,ny-1 do printf,outlun,tmp[*,i],format='(255e16.7)'
         close,outlun
         free_lun,outlun
      #endif

   endfor


   #close,/all
   #retall

   end

if __name__ == '__main__':
   print('Running __name__ as __main__')

