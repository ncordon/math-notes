function value = TPower(center, t, degree)
  if t < 0
    value = 0;
  elseif t > center
    value = 0;
  else
    value = (center - t)^degree;
  endif
endfunction