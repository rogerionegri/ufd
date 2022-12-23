FUNCTION THRES_DELTA_NBR, img, xi

  thres = fix(img[*,*] * 0)
  pos = where(img le xi)
  thres[pos] = 1

  Return, thres
END