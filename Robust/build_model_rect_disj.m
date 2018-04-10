function mod = build_model_rect_disj(dat)

% See run_algo_ab_new for dictionary and commenting conventions

%% Set (local) constants for the convenience of less typing

m          = dat.m;
n          = dat.n;
w          = dat.w;
GL         = dat.GL;
GU         = dat.GU;
f_norm_rhs = dat.f_norm_rhs;
e          = dat.eps_bor;
d          = dat.delta;
f_min      = dat.f_min;
f_max      = dat.f_max;
fL         = dat.fL;
fU         = dat.fU;
status     = dat.status;
u_a        = dat.u_a;
l_a        = dat.l_a;
u_b        = dat.u_b;
l_b        = dat.l_b;

%% Build basic model (variables and constraints) that relate f and G

mod_fGw = build_model_fGw(dat);

%% Pull the variables from the basic model into the current workspace

f = mod_fGw.f;
G = mod_fGw.G;
if dat.norm == 0
    above_GL = mod_fGw.above_GL;
end

%% Also pull in the constraints

con = mod_fGw.con;

%% Declare a and b variables, which are our primary variables that we do
%% B&B (branch and bound) over

a = sdpvar(1);
b = sdpvar(1);

%% Declare the overall tau and tau_eps variables (as described in the
%% 2016-07-14 paper) but with the following change: tau -> pi

pi  = sdpvar(m,1);
pie = sdpvar(m,1);

%% Declare the variables that are in the homegenization of the three
%% polyhedra

hom_a   = sdpvar(m,3); % Corresponds to a in paper
hom_ae  = sdpvar(m,3); % Corresponds to a_eps
hom_f   = sdpvar(m,3); % Corresponds to \varphi in paper
hom_fe  = sdpvar(m,3); % Corresponds to \varphi_\eps
hom_pi  = sdpvar(m,3); % Corresponds to tau in paper
hom_pie = sdpvar(m,3); % Corresponds to tau_\eps
hom_t   = binvar(m,3); % The homogenization variable in paper
hom_te  = binvar(m,3); % The homogenization variable for \eps

%% Setup constraints (lots of them!)

% S: Can we also add pie >= pi? Think so

% Note union(ind_a, ind_b) is all students

ind_a = find(status == 1); % Students that might be borderline wrt a
ind_b = find(status == 0); % Students that might be borderline wrt b

% All a related variables sum to a
% All b related variables sum to b
% All f related variables sum to f or f+\epsilon
% All tau related variables sum to tau or tau_\epsilon

con = [ con ; a   == sum(hom_a (ind_a,:),2) ];
con = [ con ; a   == sum(hom_ae(ind_a,:),2) ];
con = [ con ; b   == sum(hom_a (ind_b,:),2) ];
con = [ con ; b   == sum(hom_ae(ind_b,:),2) ];
con = [ con ; f   == sum(hom_f  ,2) ];
con = [ con ; f+e == sum(hom_fe ,2) ];
con = [ con ; pi  == sum(hom_pi ,2) ];
con = [ con ; pie == sum(hom_pie,2) ];

%% Sum artificial t variables 1

con = [ con ; sum(hom_t,2) == 1 ; sum(hom_te,2) == 1 ];

%% Add some valid inequalities

% If the third tau variables is 1, then so is the third tau_\epsilon
% variable

% If the third f variable is positive, then so is the third f+\epsilon
% variable

% If the third artificial t variable is positive, so is the third
% artificial t_\epsilon variable

% S: These are experimental, but I think they are OK

con = [ con ; hom_pi(:,3) <= hom_pie(:,3) ];
con = [ con ; hom_f(:,3)  <= hom_fe(:,3)  ];
con = [ con ; hom_t(:,3)  <= hom_te(:,3)  ];

% We repeat the preceding logic, except now a positive second variable
% implies that at least one of the second or third variables are
% positive

% S: These are experimental, but I think they are OK

con = [ con ; hom_pi(:,2) <= hom_pie(:,2) + hom_pie(:,3) ];
con = [ con ; hom_f(:,2)  <= hom_fe(:,2)  + hom_fe(:,3)  ];
con = [ con ; hom_t(:,2)  <= hom_te(:,2)  + hom_te(:,3)  ];

%% Build the defining polyhedral relationships between the variables using
%% function build_model_pi defined below

for i = 1:m
    if status(i) == 1
        con = build_model_pi(con, l_a, hom_a(i,:), u_a, fL(i), ...
            hom_f(i,:), fU(i), hom_pi(i,:), hom_t(i,:), d);
        con = build_model_pi(con, l_a, hom_ae(i,:), u_a, fL(i) + e, ...
            hom_fe(i,:), fU(i) + e, hom_pie(i,:), hom_te(i,:), d);
    else
        con = build_model_pi(con, l_b, hom_a(i,:), u_b, fL(i), ...
            hom_f (i,:), fU(i), hom_pi(i,:), hom_t(i,:), d);
        con = build_model_pi(con, l_b, hom_ae(i,:), u_b, fL(i) + e, ...
            hom_fe(i,:), fU(i) + e, hom_pie(i,:), hom_te(i,:), d);
    end
end

%% Set the objective to be the sum of the differences between the tau's
%% (except that the tau's are called pi)

obj = sum(pie - pi);

%% Store model in the 'mod' structure for output to calling function

mod.a   = a;
mod.b   = b;
mod.pi  = pi;
mod.pie = pie;

mod.hom_a   = hom_a  ;
mod.hom_ae  = hom_ae ;
mod.hom_f   = hom_f  ;
mod.hom_fe  = hom_fe ;
mod.hom_pi  = hom_pi ;
mod.hom_pie = hom_pie;
mod.hom_t   = hom_t  ;
mod.hom_te  = hom_te ;

mod.f = f;
mod.G = G;
if dat.norm == 0
  mod.above_GL = above_GL;
end

mod.con = con;
mod.obj = obj;

return

function con = build_model_pi(con, al, a, au, fl, f, fu, pi, t, d)

    % These are the three polyhedra

    con = [ con ; al*t(1) <= a(1) <= au*t(1) ];
    con = [ con ; fl*t(1) <= f(1) <= fu*t(1) ];
    con = [ con ; pi(1) == 0 ];
    con = [ con ; f(1) <= a(1) - d*t(1) ];
    
    con = [ con ; al*t(2) <= a(2) <= au*t(2) ];
    con = [ con ; fl*t(2) <= f(2) <= fu*t(2) ];
    con = [ con ; pi(2) == (1/d)*(f(2) - a(2)) + t(2) ];
    con = [ con ; a(2) - d*t(2) <= f(2) <= a(2) ];
    
    con = [ con ; al*t(3) <= a(3) <= au*t(3) ];
    con = [ con ; fl*t(3) <= f(3) <= fu*t(3) ];
    con = [ con ; pi(3) == t(3) ];
    con = [ con ; a(3) <= f(3) ];

return
