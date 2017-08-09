#pro mkTRKdat,inroot
import numpy as np
import os
from astropy.io import fits

#convert 256x256 section of TRK image to FGSES .dat file

def rebin(a, *args):
    '''rebin ndarray data into a smaller ndarray of the same rank whose dimensions
    are factors of the original dimensions. eg. An array with 6 columns and 4 rows
    can be reduced to have 6,3,2 or 1 columns and 4,2 or 1 rows.
    example usages:
    >>> a=rand(6,4); b=rebin(a,3,2)
    >>> a=rand(6); b=rebin(a,2)
    From: http://scipy-cookbook.readthedocs.io/items/Rebinning.html
    '''
    shape = a.shape
    lenShape = len(shape)
    factor = asarray(shape)/asarray(args)
    evList = ['a.reshape('] + \
             ['args[%d],factor[%d],'%(i,i) for i in range(lenShape)] + \
             [')'] + ['.sum(%d)'%(i+1) for i in range(lenShape)] + \
             ['/factor[%d]'%i for i in range(lenShape)]
    print ''.join(evList)
    return eval(''.join(evList))

# 256x256 section of TRK image
  data=float(mrdfits(inroot+'.fits',/unsigned,/silent))

  ; oversample by factor of 6
  data=rebin(data,1536,1536)

  ;extract 256x256 section of image
  x=768
  y=768
  trk=data[x-127:x+128,y-127:y+128]
  ;writefits,inroot+'.fits',trk

  ; write ASCII float FGSES .dat file
  get_lun,outlun
  openw,outlun,inroot+'.dat',width=4080,/swap_if_little_endian
  ;for j=0,256-1 do printf,outlun,trk[*,j],format='(256e16.7)'
  nx=256 & ny=256
  for ix=0,nx-1 do begin
     tmp=data[ix,*]
     nn=n_elements(tmp)-1L
     for i=0L,nx-1 do begin     ; convert to ascii float format
        printf,outlun,tmp[i],format='($,255e16.7," ")'
     endfor
  endfor
  close,outlun
  free_lun,outlun

end
