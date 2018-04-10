function dat = preprocess(dat);

  % Save data for convenience

  w       = dat.w;
  m       = dat.m;     % Will likely be updated
  f_nom   = dat.f_nom; % Will likely be updated
  G_nom   = dat.G_nom; % Will likely be updated
  GL      = dat.GL;    % Will likely be updated
  GU      = dat.GU;    % Will likely be updated
  l_a     = dat.l_a;
  l_b     = dat.l_b;
  u_a     = dat.u_a;
  u_b     = dat.u_b;
  eps_bor = dat.eps_bor;

  % Calculate fL = f_nom

  fL = G_nom*w;

  % Calculate a priori upper bounds on each of the f(i)'s.

  fU = get_fU(dat);

  % So that the maxgap can be calculated correctly in all cases, we will
  % make sure to keep the two (nominal) students directly below l_a and
  % l_b.

  keep_special_a = find(fL < l_a);
  [~,ind]        = max(fL(keep_special_a));
  keep_special_a = keep_special_a(ind);

  keep_special_b = find(fL < l_b);
  [~,ind]        = max(fL(keep_special_b));
  keep_special_b = keep_special_b(ind);

  if numel(keep_special_a) == 0 | numel(keep_special_b) == 0
    error('Problem in preprocessing.');
  end
  %Veronica's question: what if we have more than one?Shouldn't we keep all
  %of them?s
  keep_special = union(keep_special_a,keep_special_b);

  % Find those students who could never be borderline wrt a *and* who
  % are always above u_b. These are the students that definitely get a
  % B. Remove them.

  remove = intersect( find(fU + eps_bor < l_a) , find(fL >= u_b) );
  keep = setdiff(1:m,remove);
  keep = union(keep, keep_special);
  G_nom = G_nom(keep,:);
  GL = GL(keep,:);
  GU = GU(keep,:);
  fL = fL(keep);
  fU = fU(keep);
  m = size(G_nom,1);

  % Recalculate keep_special since G_nom has changed

  keep_special_a = find(fL < l_a);
  [~,ind] = max(fL(keep_special_a));
  keep_special_a = keep_special_a(ind);

  keep_special_b = find(fL < l_b);
  [~,ind] = max(fL(keep_special_b));
  keep_special_b = keep_special_b(ind);

  if numel(keep_special_a) == 0 | numel(keep_special_b) == 0
    error('Problem in preprocessing.');
  end

  keep_special = union(keep_special_a,keep_special_b);

  % Find those students who could never be borderline wrt b. Remove them.
  % BUT keep at least one below l_b so that gaps can be calculated.

  remove = find(fU + eps_bor < l_b);
  keep = setdiff(1:m,remove);
  keep = union(keep, keep_special);
  G_nom = G_nom(keep,:);
  GL = GL(keep,:);
  GU = GU(keep,:);
  fL = fL(keep);
  fU = fU(keep);
  m = size(G_nom,1);

  % Recalculate keep_special since G_nom has changed

  keep_special_a = find(fL < l_a);
  [~,ind] = max(fL(keep_special_a));
  keep_special_a = keep_special_a(ind);

  keep_special_b = find(fL < l_b);
  [~,ind] = max(fL(keep_special_b));
  keep_special_b = keep_special_b(ind);

  if numel(keep_special_a) == 0 | numel(keep_special_b) == 0
    error('Problem in preprocessing.');
  end

  keep_special = union(keep_special_a,keep_special_b);

  % Find those students who definitely get an A. Remove them.

  remove = find(fL >= u_a);
  keep = setdiff(1:m,remove);
  keep = union(keep, keep_special);
  G_nom = G_nom(keep,:);
  GL = GL(keep,:);
  GU = GU(keep,:);
  fL = fL(keep);
  fU = fU(keep);
  m = size(G_nom,1);

  % Recalculate keep_special since G_nom has changed

  keep_special_a = find(fL < l_a);
  [~,ind] = max(fL(keep_special_a));
  keep_special_a = keep_special_a(ind);

  keep_special_b = find(fL < l_b);
  [~,ind] = max(fL(keep_special_b));
  keep_special_b = keep_special_b(ind);

  if numel(keep_special_a) == 0 | numel(keep_special_b) == 0
    error('Problem in preprocessing.');
  end

  keep_special = union(keep_special_a,keep_special_b);

  % For each student, determine whether they *could* be borderline with
  % respect to a or b

  status = [];
  for i = 1:m
    if fL(i) <= u_b && l_b <= fU(i)+eps_bor
      status = [status;0];
    elseif fL(i) <= u_a && l_a <= fU(i)+eps_bor
      status = [status;1];
    elseif i == keep_special_a
      status = [status;1];
    elseif i == keep_special_b
      status = [status;0];
    else
      status = [status;1];
      [fL(i),fU(i)]
      [l_b,u_b]
      [l_a,u_a]
      % keyboard
      % warning('Preprocess got unexpected');
      error('Preprocess got unexpected');
    end
  end

  % Save everything before returning

  dat.m      = m;
  dat.f_nom  = fL;
  dat.G_nom  = G_nom;
  dat.GL     = GL;
  dat.GU     = GU;
  dat.fL     = fL;
  dat.fU     = fU;
  dat.f_min  = min(fL);
  dat.f_max  = max(fU);
  dat.status = status;

return

function fU = get_fU(dat)

  % Build basic model compenents for f = G*w, GL <= G <= GU, # changed
  % entries in G <= f_norm_rhs

  mod_fGw = build_model_fGw(dat);

  % Initialize fU to f_max

  fU = dat.f_max*ones(dat.m,1);

  settings = set_solver_options;

  for i = 1:dat.m
    
    diagnostics = solvesdp( mod_fGw.con , -mod_fGw.f(i) , settings );

    if diagnostics.problem ~= 0
      if diagnostics.problem == 1
        str = 'In pre_process/fU, CPLEX appears to have encountered an infeasible problem.';
        warning(str);
        error(str);
      else 
        str = 'In pre_process/fU, CPLEX did not solve subproblem within time limit.';
        warning(str);
        error(str);
      end
    end

    fU(i) = double( mod_fGw.f(i) );

  end

return
