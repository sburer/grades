function mod_fGw = build_model_fGw(dat)

% See run_algo_ab_new for dictionary and commenting conventions

%% Set (local) constants for the convenience of less typing

m          = dat.m;
n          = dat.n;
w          = dat.w;
tol        = dat.tol_gen;
GL         = dat.GL;
GU         = dat.GU;
f_norm_rhs = dat.f_norm_rhs;

%% Setup variables f and G

f = sdpvar(m,1);
G = sdpvar(m,n);

%% Initialize constraint structure

con = [];

%% Add f = G*w constraint

con = [ con ; f == G*w ];

%% Add entry-wise lower and upper bounds on G

con = [ con ; GL <= G <= GU ]; % S: For YALMIP, is GL(:) <= G(:) <= GU(:) safer?

%% If the chosen norm is 1, then...

if dat.norm == 1

    %% Add constraint to limit the total amount that G can change above
    %% GL

    con = [ con ; sum(sum(G - GL)) <= f_norm_rhs*sum(sum(GU - GL)) ];

%% Else, if the chosen norm is 0, then...

elseif dat.norm == 0

    %% Add binary variables and constraints to limit the total number of
    %% entries that G can change compared to GL

    above_GL = binvar(m,n);
    con = [ con ; tol*above_GL <= G - GL <= above_GL .* (GU - GL) ];
    con = [ con ; sum(sum(above_GL)) <= f_norm_rhs ]; 

    %% Also add limits on the number of changes per row

    %con = [ con ; sum(above_GL, 2) <= 3 ];
    

else

    error('Unexpected dat.norm value');

end

%% Return YALMIP variables

mod_fGw.f = f;
mod_fGw.G = G;
if dat.norm == 0
  mod_fGw.above_GL = above_GL;
end

% Return YALMIP constraints

mod_fGw.con = con;
