function settings = set_cplex_options;

settings = sdpsettings;

settings.verbose = 0;
settings.solver = 'gurobi';

settings.gurobi.Threads = 3;
settings.gurobi.TimeLimit = 300;
settings.gurobi.OutputFlag = 0;
% settings.gurobi.PrePasses = 5;
settings.gurobi.Cuts = 3;

settings.cplex.timelimit = 60;
settings.cplex.threads = 3;
% settings.cplex.mip.tolerances.integrality = 0.1;
% settings.cplex.emphasis.mip = 4; % This is important!

% Clear YALMIP Gurobi option that is no longer available in Gurobi

settings.gurobi = rmfield(settings.gurobi, 'PreMIQPMethod');
