clear

cd mat/

addpath ..

infiles = dir('output3*.mat');

for ind = 1:length(infiles)

  load(infiles(ind).name);

  results = {};

  for iter = 1 : length(data_us)

    iter

    dat = data_us{iter};

    % Traditional

    mod = build_model_rect_disj(dat);
    mod = extract_gurobi_model(dat,mod);

    a = dat.u_a;
    b = dat.u_b;
    numbor     = calc_bound_rect_disj(dat,mod,a,a,b,b,dat.fL,dat.fU,'upper');
    numbor_nom = calc_bound_rect_disj(dat,mod,a,a,b,b,dat.fL,dat.fL,'upper');
    aa = a-dat.tol_gen/10;
    bb = b-dat.tol_gen/10;
    maxgap     = get_maxgap_ubV(aa,aa,bb,bb,dat.fL);

    output_trad{iter} = [a,b,numbor,maxgap,numbor_nom];

    % Maxgap

    [maxgap,a,b] = get_maxgap_ubV(dat.l_a, dat.u_a, dat.l_b, dat.u_b, dat.fL);
    numbor       = calc_bound_rect_disj(dat,mod,a,a,b,b,dat.fL,dat.fU,'upper');
    numbor_nom   = calc_bound_rect_disj(dat,mod,a,a,b,b,dat.fL,dat.fL,'upper');

    output_maxgap{iter} = [ a, b, numbor, maxgap, numbor_nom ];

    % Basic

    % [numbor_nom,~,retcode,a,b] = calc_bound_rect_disj(...
    %   dat, mod, dat.l_a, dat.u_a, dat.l_b, dat.u_b, dat.fL, dat.fL, 'lower');
    % if retcode == 1
    %   numbor_nom = -Inf;
    %   maxgap = -Inf;
    %   numbor = -Inf;
    % else
    %   aa = a-dat.tol_gen/10;
    %   bb = b-dat.tol_gen/10;
    %   maxgap = get_maxgap_ubV(aa,aa,bb,bb,dat.fL);
    %   numbor = calc_bound_rect_disj(dat,mod,a,a,b,b,dat.fL,dat.fU,'upper');
    % end

    fprintf('Trying to run basic with maxgap...\n');
    save_dat_fU = dat.fU;
    dat.fU = dat.fL;
    [gub,a_opt,b_opt,gap_opt,feedback] = run_algo_ab_new(dat);

    % [numbor_nom,gub]
    % [maxgap,gap_opt]

    % maxgap = gap_opt;
    % numbor_nom = gub;

    a = a_opt(1); % Because we were getting multiple solutions from above
    b = b_opt(1);
    dat.fU = save_dat_fU;
    numbor = calc_bound_rect_disj(dat,mod,a,a,b,b,dat.fL,dat.fU,'upper');
    maxgap = gap_opt;
    numbor_nom = gub;

    output_basic{iter} = [a,b,numbor,maxgap,numbor_nom];

    % numbor_nom for our method

    a = output_us{iter}(1);
    b = output_us{iter}(2);
    numbor_nom   = calc_bound_rect_disj(dat,mod,a,a,b,b,dat.fL,dat.fL,'upper');
    output_us{iter}(5) = numbor_nom;

    % Summarize

    results{iter} = ...
      [output_trad{iter};output_maxgap{iter};output_basic{iter};[output_us{iter}]];

  end

  save(strrep(infiles(ind).name,'output','results'),...
    'data_us', 'output_us', 'feedback_us', ...
    'output_trad', 'output_maxgap', 'output_basic', ...
    'results');

end 

cd ..
