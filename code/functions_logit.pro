PRO functions_logit
   print, 'importing...'
END


;---------------------------------
FUNCTION sigmoid, x
  return, 1.0/(1.0 + EXP(-x))
END


;---------------------------------
FUNCTION init_pars, dim
  w = dblarr(dim)
  b = 0D
  return, {w:w, b:b} 
END

;---------------------------------
FUNCTION accuracy, y_pred, y_true
  
  len = n_elements(y_pred)
  vecPred = fltarr(len)
  for i = 0, len-1 do begin
    if y_pred[i] gt 0.5 then vecPred[i] = 1.0 else vecPred[i] = 0.0 
  endfor

  res = 1.0 - total(abs(vecPred - y_true))/len

  return, res
END

;---------------------------------
FUNCTION predict, x, w, b
  ;z = X ## W + b
  z = (x[*] * w[0]) + b
  y_hat = sigmoid(z)
  return, y_hat
END


;---------------------------------
FUNCTION train_logit, x, y, rate, iters

  if n_elements(size(x,/dimensions)) eq 1 then begin
    dim = 1
    len = N_ELEMENTS(x)
  endif else begin
    dim = N_ELEMENTS(x[*,0])
    len = N_ELEMENTS(x[0,*])
  endelse
  
  
  stru = init_pars(dim)
  W = stru.w
  b = stru.b
  
  vecLoss = [-1]
  
  for it = 0L, iters do begin
    
    vecHats = dblarr(len)
    loss = 0D

    z = (x[*] * w[0]) + b
    vecHats = sigmoid(z)
    loss = (-1.0/len) * total( y * alog(vecHats) + (- y[*] + 1) * alog(- vecHats[*] + 1 ) )    
    
    vecLoss = [vecLoss , loss]
    
    dW = total( x * (vecHats - y) )/len
    db = mean(vecHats - y)
    
    W -= rate * dW
    b -= rate * db
    
    ;plot, vecHats
    if it mod 100 eq 0 then begin
      print, 'it/loss:', it, loss
      xx = (indgen(1000)/1000.0) & yy = xx*0D &  for i = 0, 999 do yy[i] = sigmoid((xx[i] * w[0]) + b) ;& plot, xx,yy
    endif
    
  endfor

  return, {w:w, b:b, loss:vecLoss}
END




;---------------------------------
FUNCTION class_logit_prob, transfImage, modLogit
  cla = FIX(transfImage*0)
  prob = (transfImage*0)
  nc = n_elements(transfImage[*,0])
  nl = n_elements(transfImage[0,*])
  prob_img = fltarr(nc,nl,2)
  
  len = n_elements(transfImage)
  
  for i = 0, nc-1 do begin
    for j = 0, nl-1 do begin
      prob[i,j] = sigmoid( (transfImage[i,j] * modLogit.w[0]) + modLogit.b )
      prob_img[i,j,1] = prob[i,j]
      prob_img[i,j,0] = 1.0 - prob[i,j] 
       
      if prob[i,j] gt 0.5 then cla[i,j] = 1.0 else cla[i,j] = 0.0
    endfor
  endfor

  return, {class: cla, prob: prob, prob_img: prob_img}
END