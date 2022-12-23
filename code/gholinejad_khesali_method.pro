PRO Gholinejad_Khesali_method

  ;Input data
  Path = './dataset/MODIS/SAL_deltaNBR_ago-out20.tif' ;DeltaNBR
  Path_base = './dataset/MODIS/SAL_jun19-jun20.tif'   ;SPAM/Median
  bands = [0,1,2,3,4,5,6,7,8,9,10]                    ;Selected bands/instants
  bandNBR_spam = 13                                   ;Band containing the reference NBR
  xi = 0.185                                          ;Non-change threshold (xi)
  path_out = './dataset/output/'                      ;Output path
  ;--------------------------------------------------------------------------
  
  bandsDNBR  = [0,1,2,3]                            ;SAL-DNBR
  bandsNBR  = [4,5,6,7]                             ;SAL-NBR
  bandsNDVI = [8,9,10,11]                           ;SAL-NDVI
  bandsNDWI = [12,13,14,15]                         ;SAL-NDWI

  band_baselineNBR_NDVI_NDWI = [19,20,21]           ;SAL-DeltaRef

  tresNDWI=0.2
  tresNDVI=-0.1
  tresDNBR=0.4;

  ;--------------------------------------------------------------------------
   
  t0 = systime(/seconds)

  imgBase = read_tiff(Path_base)
  imgBaseRef = imgBase[band_baselineNBR_NDVI_NDWI,*,*]
  
  img = read_tiff(Path)
  imgNBR = img[bandsNBR,*,*]
  imgNDVI = img[bandsNDVI,*,*]
  imgDNBR = img[bandsDNBR,*,*]
  imgNDWI = img[bandsNDWI,*,*]

  nb = n_elements(imgNBR[*,0,0])
  nc = n_elements(imgNBR[0,*,0])
  nl = n_elements(imgNBR[0,0,*])
  
  maxImage = fltarr(nc,nl)
  maskWater = INTARR(nc,nl)
  maskVegeta = INTARR(nc,nl)
  imgWater = FLTARR(nc,nl)
  imgVegeta = FLTARR(nb,nc,nl)

  for i = 0, nc-1 do begin
    for j = 0, nl-1 do begin
      
      if min(imgNDWI[*,i,j]) GT tresNDWI then maskWater[i,j] = 1 
      for k = 0, nb-1 do imgVegeta[k,i,j] = imgNDVI[k,i,j] - imgBaseRef[1,i,j]
      if min(imgVegeta[*,i,j]) LT tresNDVI then maskVegeta[i,j] = 1
      
      maxImage[i,j] = max(imgDNBR[*,i,j])
    endfor
  endfor

  ;...
  pos = where(maskVegeta*(~maskWater) ne 0)
  imgFil = maxImage[pos]
  
  imgOrd = imgFil[ SORT(imgFil) ]
  m = imgOrd[ UNIQ(imgOrd) ]
  posBurn = where(maxImage GT tresDNBR)
  GKMap = (maskWater[*,*]*0)   &   GKMap[posBurn] = 1
  
  plot, m
  tvscl, rotate(~GKMap,7)
 
  ;Dimension adjusts, if necessary
  Result = QUERY_TIFF(Path_base, Info, GEOTIFF=geoVar)
  dimRef = info.DIMENSIONS
  Image_GK = CONGRID(~GKMap, dimRef[0], dimRef[1])
  
  ;Saving the results...
  write_tiff, path_out + 'Res_GKM__' + strtrim(string(tresNDWI))+'-'+strtrim(string(tresNDVI))+'-'+strtrim(string(tresDNBR)) + (Path.split('/'))[-1], geotiff=geoVar, Image_GK
  ;--------------------------------------------------------------

  t1 = systime(/seconds)

  print, 'tempo gkm: ', t1-t0
  print, 'fim...'
  stop

END