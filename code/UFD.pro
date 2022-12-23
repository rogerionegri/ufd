@functions_logit.pro

@ICM.pro
@ML.pro
@ICM_ADDS.pro

PRO UFD

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

  transfImage = fltarr(nc,nl)
  maxImage = fltarr(nc,nl)
  minImage = fltarr(nc,nl)

  t0 = systime(/seconds)

  print,'Computing the divergence values...'
  
  ;Point-to-line distance (divergence measure)
  W = fltarr(nb) + 1
  W *= 1/norm(W)
  WT = transpose(W)
  nW = norm(W)

  for i = 0, nc-1 do begin
    for j = 0, nl-1 do begin
      projEscalar = ( img[*,i,j] ## WT ) / nW
      proj =  W # projEscalar
      d = norm(img[*,i,j] - proj)

      transfImage[i,j] = d

      minImage[i,j] = min(abs(img[*,i,j]))
      maxImage[i,j] = max(abs(img[*,i,j]))
    endfor
  endfor

  seedLow = where(maxImage le xi)
  vecLowDist = transfImage[seedLow]

  seedHigh = where(minImage gt xi)
  vecHighDist = transfImage[seedHigh]

  ;Lower/Upper envelope limit for selecting non-burned pixels/examples
  Lub = max(vecLowDist)  + 0.5*stddev(vecLowDist)
  Hlb = Lub


  print,'Adjusting the Logit model...'

  posOne = where(transfImage le Lub)
  posZero = where(transfImage ge Hlb)

  ;preparar o dado...
  x_one = transfImage[posOne]
  y_one = intarr(n_elements(x_one)) +1

  x_zero = transfImage[posZero]
  y_zero = intarr(n_elements(x_zero))

  x_train = [x_one , x_zero]
  y_train = [y_one , y_zero]

  modLogit = train_logit(x_train, y_train, 1.5, 10000)

  classModLogit = class_logit_prob(transfImage, modLogit)

  print,'Applying the ICM algorithm...'
  
  imgc = ICM(classModLogit.class+1,1,10,2,classModLogit.prob_img,'NULL')

  t1 = systime(/seconds)

  ;Dimension adjusts, if necessary
  Result = QUERY_TIFF(Path_base, Info, GEOTIFF=geoVar)
  dimRef = info.DIMENSIONS
  imgc_CG = CONGRID(imgc, dimRef[0], dimRef[1])
  classModLogit_CG = CONGRID(classModLogit.class, dimRef[0], dimRef[1])
  transfImage_CG = CONGRID(transfImage, dimRef[0], dimRef[1])
  
  maxImage_CG = CONGRID(maxImage, dimRef[0], dimRef[1])
  maxRelativeImage_CG = CONGRID(maxImage_CG[*,*]/sqrt(abs(imgBaseNBR[*,*])), dimRef[0], dimRef[1])

  ;Saving the results...
  print,'Saving the results...'
  write_tiff, path_out + 'UFD+MRF__' +STRTRIM(STRING(xi),1)+ (Path.split('/'))[-1], geotiff=geoVar, imgc_CG-1
  write_tiff, path_out + 'UFD__' + (Path.split('/'))[-1], geotiff=geoVar, classModLogit_CG  ;classModLogit.class
  write_tiff, path_out + 'TransfImage__' + (Path.split('/'))[-1], geotiff=geoVar, transfImage_CG, /FLOAT

  print, 'UFD run-time: ', t1-t0
  print, 'End of process!'

END