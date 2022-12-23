@thres_delta_nbr.pro

PRO TDNBR

  ;Input data
  Path = './dataset/MODIS/SAL_deltaNBR_ago-out20.tif' ;DeltaNBR
  Path_base = './dataset/MODIS/SAL_jun19-jun20.tif'   ;SPAM/Median
  bands = [0,1,2,3,4,5,6,7,8,9,10]                    ;Selected bands/instants
  bandNBR_spam = 13                                   ;Band containing the reference NBR 
  xi = 0.185                                          ;Non-change threshold (xi)
  path_out = './dataset/output/'                      ;Output path
  ;--------------------------------------------------------------------------

  print,'Data reading stage...'
  
  imgBase = read_tiff(Path_base)
  imgBaseNBR = FLTARR(n_elements(imgBase[0,*,0]),n_elements(imgBase[0,0,*]))
  imgBaseNBR[*,*] = imgBase[bandNBR_spam,*,*]

  img = read_tiff(Path)
  img = img[bands,*,*]

  nb = n_elements(img[*,0,0])
  nc = n_elements(img[0,*,0])
  nl = n_elements(img[0,0,*])

  maxImage = fltarr(nc,nl)

  t0 = systime(/seconds)

  for i = 0, nc-1 do begin
    for j = 0, nl-1 do begin
      maxImage[i,j] = max(abs(img[*,i,j]))
    endfor
  endfor

  ;Dimension adjusts, if necessary
  Result = QUERY_TIFF(Path_base, Info, GEOTIFF=geoVar)
  dimRef = info.DIMENSIONS
  maxImage_CG = CONGRID(maxImage, dimRef[0], dimRef[1])
  maxRelativeImage_CG = CONGRID(maxImage_CG[*,*]/sqrt(abs(imgBaseNBR[*,*])), dimRef[0], dimRef[1])

  print,'Applying the TDNBR rule...'
  thres_maxImage_CG = THRES_DELTA_NBR(maxRelativeImage_CG,xi)
  t1 = systime(/seconds)

  ;Saving the results...
  print,'Saving the results...'
  write_tiff, path_out + 'Res_ThresRelativeNBR__'+ (Path.split('/'))[-1], geotiff=geoVar, thres_maxImage_CG
  write_tiff, path_out + 'Res_MaxDeltaNBR__' + (Path.split('/'))[-1], geotiff=geoVar, maxImage_CG, /Float

  ;--------------------------------------------------------------

  print, 'TDNBR run-time: ', t1-t0
  print, 'End of process!'

END