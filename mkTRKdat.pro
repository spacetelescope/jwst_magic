pro mkTRKdat,inroot
  ; convert 256x256 section of TRK image to FGSES .dat file

  ; 256x256 section of TRK image
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
