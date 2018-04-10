function num_pseudobor = calc_num_pseudobor(dat,f,a,b)

% Save info for convenience

m = dat.m;
e = dat.eps_bor;
d = dat.delta;

% Calculate number of pseudo borderline

num_pseudobor = 0.0;

for i = 1:m

  if f(i) <= b - e - d
    num_pseudobor = num_pseudobor + 0;
  elseif f(i) <= b - e
    num_pseudobor = num_pseudobor + (1/d)*(f(i) - (b - e - d));
  elseif f(i) <= b - d
    num_pseudobor = num_pseudobor + 1;
  elseif f(i) <= b
    num_pseudobor = num_pseudobor + 1 - (1/d)*(f(i) - (b - d));
  elseif f(i) <= a - e - d
    num_pseudobor = num_pseudobor + 0;
  elseif f(i) <= a - e
    num_pseudobor = num_pseudobor + (1/d)*(f(i) - (a - e - d));
  elseif f(i) <= a - d
    num_pseudobor = num_pseudobor + 1;
  elseif f(i) <= a
    num_pseudobor = num_pseudobor + 1 - (1/d)*(f(i) - (a - d));
  elseif f(i) <= 100
    num_pseudobor = num_pseudobor + 0;
  else
    error('get_num_pseudobor: f(i) was above 100');
  end

end
