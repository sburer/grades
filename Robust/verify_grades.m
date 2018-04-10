clear

%load mat/output3_norm0_eps1.0_unc0.2.mat
%load mat/tighter_tol_int/output3_norm0_eps0.5_unc0.2.mat
%load mat/tighter_tol_int/output3_norm0_eps0.5_unc0.2.mat
%load mat/old/output3_norm0_eps0.5_unc0.2.mat
%load mat/old/output3_norm0_eps1.0_unc0.2.mat
load mat/old/results3_norm1_eps1_unc0.2.mat
load mat/results3_norm0_eps0.5_unc0.2.mat

for iter = 30 : 30 % 1 : length(output_us)

  dat = data_us{iter};

  a = output_us{iter}(1);
  b = output_us{iter}(2);
  gub = output_us{iter}(3);
  maxgap = output_us{iter}(4);

  mod = build_model_rect_disj(dat);
  mod = extract_gurobi_model(dat,mod);
  [bd,f] = calc_bound_rect_disj(dat,mod,a,a,b,b,dat.fL,dat.fU,'upper');
  tmp = abs(gub - calc_num_pseudobor(dat,f,a,b));
  calc_num_pseudobor(dat,f,a,b)
  if tmp > dat.tol_gen
    fprintf('iter = %d did not verify GUB (%e)\n', iter, tmp);
  end

  a = a - dat.tol_gen/10;
  b = b - dat.tol_gen/10;

  [gapcheck,~,~] = get_maxgap_ubV(a,a,b,b,dat.f_nom);
  tmp2 = abs( maxgap - gapcheck );
  if tmp2 > dat.tol_gen
    fprintf('iter = %d did not verify MAXGAP (%.3f  %.3f  %e)\n',...
      iter, gapcheck, maxgap, tmp2);
  end

end
