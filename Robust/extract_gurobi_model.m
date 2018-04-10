function gur = extract_gurobi_model(dat,mod)

a_rand = 100*rand;
b_rand = 100*rand;
f_rand = 100*rand(dat.m,1);

mod.con = [ mod.con ; a_rand <= mod.a <= a_rand ];
mod.con = [ mod.con ; b_rand <= mod.b <= b_rand ];
mod.con = [ mod.con ; f_rand <= mod.f <= f_rand ];

settings = set_solver_options;
settings.savesolverinput = 1;
settings.gurobi.TimeLimit = 5;

diagnostics = solvesdp(mod.con,-mod.obj,settings);

gur = diagnostics.solverinput.model;
gur.a_lb = find(gur.rhs == -a_rand);
gur.a_ub = find(gur.rhs ==  a_rand);
gur.b_lb = find(gur.rhs == -b_rand);
gur.b_ub = find(gur.rhs ==  b_rand);
gur.f_lb = [];
gur.f_ub = [];
for i = 1:dat.m
  gur.f_lb = [ gur.f_lb ; find(gur.rhs == -f_rand(i)) ];
  gur.f_ub = [ gur.f_ub ; find(gur.rhs ==  f_rand(i)) ];
end
gur.obj = -gur.obj;
[~,gur.f] = find(gur.A(gur.f_lb,:));
[~,gur.a] = find(gur.A(gur.a_lb,:));
[~,gur.b] = find(gur.A(gur.b_lb,:));

return
