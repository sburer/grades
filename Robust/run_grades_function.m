function run_grades_function(DATA, NORM, EPS, UNC, INSTANCE)

% Load data, which includes matrix tot_G_nom
% Sam: I think we mean "load mat/Gnom3"

if DATA == 3
    load mat/Gnom3
elseif DATA == 4
    load mat/Gnom4
end

% Initialize the data structures for saving the output

output_us     = [];
output_basic  = [];
output_90     = [];
output_maxgap = [];

% Solve desired instance in tot_G_nom

for iter = INSTANCE:INSTANCE % length(tot_G_nom)

  % Clear dat, the main structure for a data instance 

  clear dat

  % Set tolerances

  dat.tol_gen = 1.0e-6;     % Generic tolerance
  dat.tol_bor = 1.0e-3;     % Tolerance for borderlines, which
                            % assumes Gurobi's MIPGap is 1.0e-4 by default!
  dat.tol_gap = 1.0e-6;     % Maxgap tolerance
  dat.tol_int = 0.9*1.0e-1; % Interval tolerance

  % Set lower and upper limits on a and b

  dat.l_a = 88;
  dat.u_a = 90;
  dat.l_b = 78;
  dat.u_b = 80;

  % Set default f_min and f_max values. Will be updated later in the
  % precprocess subroutine

  dat.f_min = 0;
  dat.f_max = 100;

  % Set borderline tolerance. We have been manually setting either 1.0
  % or 0.5

  % dat.eps_bor = 1.0;
  %dat.eps_bor = 0.5;
  dat.eps_bor = EPS;
  %if dat.eps_bor == 0.5
    %warning('Note: eps_bor is now 0.5');
  %end

  % Set "pseudo-borderline" approximation tolerance. We have *not*
  % tested different values for this.

  dat.delta = 0.001;

  % Set nominal data
  
  dat.G_nom = tot_G_nom{iter};
  dat.w     = ones(size(dat.G_nom, 2), 1)/size(dat.G_nom, 2);
  dat.f_nom = dat.G_nom*dat.w;

  % Sort nominal data in terms of descending f_nom

  [dat.f_nom,perm] = sort(dat.f_nom,'descend');
  dat.G_nom        = dat.G_nom(perm,:);

  % Sanity check: largest entry of G_nom should be <= 100

  if max(max(dat.G_nom)) > dat.f_max
    error('f_nom has entry greater than 100');
  end

  % Sanity check: largest entry of f_nom should be <= 100

  if max(dat.f_nom) > dat.f_max
    error('f_nom has entry greater than 100');
  end

  % Set primary dimensions

  dat.n = length(dat.w);
  dat.m = length(dat.f_nom);

  % Calculate and set lower and upper bounds on G

  dat.GL = dat.G_nom;
  dat.GU = dat.G_nom+5; % Hard-coded. Have always been using +5
  for j = 1:dat.n
    dat.GU(:,j) = min(dat.GU(:,j),100); % But max can logically be 100
  end

  % Set the type of norm used in the uncertainty set. We have only
  % experimented with 0 and 1

  dat.norm = NORM;
  %dat.norm = 0; % 1 or 0
  %if dat.norm == 1
    %warning('Uncertainty set defined by 1-norm');
  %end

  % Calculate and set the number (or amount) of entries that will be
  % allowed to change in G. Relates to choice of norm

  if dat.norm == 0
    dat.f_norm_rhs = ceil(UNC*dat.m*dat.n);
  else
    dat.f_norm_rhs = UNC;
    % dat.f_norm_rhs = 0.0001;
    % warning('f_norm_rhs is much smaller than usual');
  end

  % Preprocess data

  dat = preprocess(dat);

  % Run algorithm

  [gub,a_opt,b_opt,gap_opt,feedback] = run_algo_ab_new(dat);

  % Sam removed a lot of old, commented code just below this line. Can
  % pull archive if necessary. 2015.04.30

  % Store and save output for our method. The 'Inf' means that the value
  % needs to be calculated later; see run_others.m

  output_us{iter} = [ a_opt(1), b_opt(1), gub , gap_opt , Inf ]; 
  data_us{iter} = dat;
  feedback_us{iter} = feedback;

  %warning('Not saving results to disk')
  filename = strcat('mat/output', num2str(DATA), '_norm', num2str(NORM), '_eps', num2str(EPS), '_unc', num2str(UNC), '_instance', num2str(INSTANCE), '.mat');
  save(filename, 'output_us', 'feedback_us', 'data_us')

  % For Sam, send text message with update on code's progress.

  if mod(iter,5) == 0
    !sendtxt.sh
  end

end
