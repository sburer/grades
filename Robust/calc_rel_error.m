function relerr = calc_rel_error(v1,v2)

  if abs(v1) == Inf | abs(v2) == Inf
    relerr = Inf;
  else
    relerr = abs(v1-v2)/max(1,0.5*(abs(v1)+abs(v2)));
  end

return
