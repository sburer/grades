#this file needs to be inserted at the link:
#https://neos-server.org/neos/solvers/milp:Gurobi/AMPL.html
#into the model file

reset;
param m; #number of students
set ind_a within {1..m};#set of students that may be borderline w.r.t to a 
set ind_b within {1..m};#set of students that may be borderline w.r.t to b

param f{1..m}>=0,<=100;#nominal grades
param eps >0 default 1;
param delta>0 default 0.01;

param l_a default 88;
param u_a >=l_a default 90;
param l_b default 78;
param u_b >=l_b default 80;
param l_phi_a{ind_a};
param u_phi_a{ind_a};
param l_phieps_a{ind_a};
param u_phieps_a{ind_a};
param l_phi_b{ind_b};
param u_phi_b{ind_b};
param l_phieps_b{ind_b};
param u_phieps_b{ind_b};

var a>=l_a,<=u_a; #breakpoint variable
var alpha_a{ind_a};
var alphaeps_a{ind_a};
var b>=l_b,<=u_b; #breakpoint variable
var alpha_b{ind_b};
var alphaeps_b{ind_b};

#variables for modeling the three polyhedrons

var t_a{ind_a,k in 1..3} binary;
var teps_a{ind_a,k in 1..3} binary;
var t_b{ind_b,k in 1..3} binary;
var teps_b{ind_b,k in 1..3} binary;
var ak{ind_a,k in 1..3};
var aepsk{ind_a,k in 1..3};
var bk{ind_b,k in 1..3};
var bepsk{ind_b,k in 1..3};
var alphak_a{ind_a, 1..3}>=0,<=1;
var alphaepsk_a{ind_a, 1..3}>=0,<=1;
var phik_a{ind_a,1..3};
var phiepsk_a{ind_a,1..3};
var alphak_b{ind_b, 1..3}>=0,<=1;
var alphaepsk_b{ind_b, 1..3}>=0,<=1;
var phik_b{ind_b,1..3};
var phiepsk_b{ind_b,1..3};

minimize pseudoBorder:sum{i in ind_a}(alphaeps_a[i]-alpha_a[i])+sum{i in ind_b}(alphaeps_b[i]-alpha_b[i]);#sum (pie[i]-pi[i])


#constraints stating the union of polyhedra

s.t. sum1a{i in ind_a}:a=sum{k in 1..3}ak[i,k];
s.t. sum1aeps{i in ind_a}:a=sum{k in 1..3}aepsk[i,k];
s.t. sum1alp{i in ind_a}:alpha_a[i] = sum{k in 1..3}alphak_a[i,k];
s.t. sum1alpeps{i in ind_a}:alphaeps_a[i] = sum{k in 1..3}alphaepsk_a[i,k];

s.t. sum1b{i in ind_b}:b=sum{k in 1..3}bk[i,k];
s.t. sum1beps{i in ind_b}:b=sum{k in 1..3}bepsk[i,k];
s.t. sum1alpb{i in ind_b}:alpha_b[i] = sum{k in 1..3}alphak_b[i,k];
s.t. sum1alpepsb{i in ind_b}:alphaeps_b[i] = sum{k in 1..3}alphaepsk_b[i,k];

#constraints on phi equal to fbar
s.t. sphi{i in ind_a}:sum{k in 1..3}phik_a[i,k]=f[i];
s.t. sphieps{i in ind_a}:sum{k in 1..3}phiepsk_a[i,k]=f[i]+eps;
s.t. sti{i in ind_a}:sum{k in 1..3}t_a[i,k]=1;
s.t. stieps{i in ind_a}:sum{k in 1..3}teps_a[i,k]=1;

s.t. sphib{i in ind_b}:sum{k in 1..3}phik_b[i,k]=f[i];
s.t. sphiepsb{i in ind_b}:sum{k in 1..3}phiepsk_b[i,k]=f[i]+eps;
s.t. stib{i in ind_b}:sum{k in 1..3}t_b[i,k]=1;
s.t. stiepsb{i in ind_b}:sum{k in 1..3}teps_b[i,k]=1;

s.t. Qkl{i in ind_a,k in 1..3}:ak[i,k]>=t_a[i,k]*l_a;
s.t. Qku{i in ind_a,k in 1..3}:ak[i,k]<=t_a[i,k]*u_a;
s.t. Qkphil{i in ind_a,k in 1..3}:phik_a[i,k]>=t_a[i,k]*l_phi_a[i];
s.t. Qkphiu{i in ind_a,k in 1..3}:phik_a[i,k]<=t_a[i,k]*u_phi_a[i];
s.t. Qkleps{i in ind_a,k in 1..3}:aepsk[i,k]>=teps_a[i,k]*l_a;
s.t. Qkueps{i in ind_a,k in 1..3}:aepsk[i,k]<=teps_a[i,k]*u_a;
s.t. Qkphiepsl{i in ind_a,k in 1..3}:phiepsk_a[i,k]>=teps_a[i,k]*l_phieps_a[i];
s.t. Qkphiepsu{i in ind_a,k in 1..3}:phiepsk_a[i,k]<=teps_a[i,k]*u_phieps_a[i];

s.t. Qklb{i in ind_b,k in 1..3}:bk[i,k]>=t_b[i,k]*l_b;
s.t. Qkub{i in ind_b,k in 1..3}:bk[i,k]<=t_b[i,k]*u_b;
s.t. Qkphilb{i in ind_b,k in 1..3}:phik_b[i,k]>=t_b[i,k]*l_phi_b[i];
s.t. Qkphiub{i in ind_b,k in 1..3}:phik_b[i,k]<=t_b[i,k]*u_phi_b[i];
s.t. Qklepsb{i in ind_b,k in 1..3}:bepsk[i,k]>=teps_b[i,k]*l_b;
s.t. Qkuepsb{i in ind_b,k in 1..3}:bepsk[i,k]<=teps_b[i,k]*u_b;
s.t. Qkphiepslb{i in ind_b,k in 1..3}:phiepsk_b[i,k]>=teps_b[i,k]*l_phieps_b[i];
s.t. Qkphiepsub{i in ind_b,k in 1..3}:phiepsk_b[i,k]<=teps_b[i,k]*u_phieps_b[i];



#constraints on each (a,f[i],alpha[i]) belonging to Q1
s.t. Qk1{i in ind_a}:phik_a[i,1]<=ak[i,1]-delta*t_a[i,1];
s.t. alp10{i in ind_a}:alphak_a[i,1]=0;

#constraints on each (b,f[i],alpha_b[i]) belonging to Q1
s.t. Qk1b{i in ind_b}:phik_b[i,1]<=bk[i,1]-delta*t_b[i,1];
s.t. alp10b{i in ind_b}:alphak_b[i,1]=0;

#constraints on each (a,f[i]+eps,alphaeps[i]) belonging to Q1
s.t. Qk1eps{i in ind_a}:phiepsk_a[i,1]<=aepsk[i,1]-delta*teps_a[i,1];
s.t. alp10eps{i in ind_a}:alphaepsk_a[i,1]=0;
#constraints on each (b,f[i]+eps,alphaeps_b[i]) belonging to Q1
s.t. Qk1epsb{i in ind_b}:phiepsk_b[i,1]<=bepsk[i,1]-delta*teps_b[i,1];
s.t. alp10epsb{i in ind_b}:alphaepsk_b[i,1]=0;

#constraints on each (a,f[i],alpha[i]) belonging to Q2
s.t. alp2{i in ind_a}:alphak_a[i,2]=(1/delta)*(phik_a[i,2]-ak[i,2]+delta*t_a[i,2]);
s.t. Qk2l{i in ind_a}:ak[i,2]-delta*t_a[i,2]<=phik_a[i,2];
s.t. Qk2lbis{i in ind_a}:phik_a[i,2]<=ak[i,2];
#constraints on each (b,f[i],alpha_b[i]) belonging to Q2
s.t. alp2b{i in ind_b}:alphak_b[i,2]=(1/delta)*(phik_b[i,2]-bk[i,2]+delta*t_b[i,2]);
s.t. Qk2lb{i in ind_b}:bk[i,2]-delta*t_b[i,2]<=phik_b[i,2];
s.t. Qk2lbisb{i in ind_b}:phik_b[i,2]<=bk[i,2];

#constraints on each (a,f[i]+eps,alphaeps[i]) belonging to Q2
s.t. alp2eps{i in ind_a}:alphaepsk_a[i,2]=(1/delta)*(phiepsk_a[i,2]-aepsk[i,2]+delta*teps_a[i,2]);
s.t. Qk2leps{i in ind_a}:aepsk[i,2]-delta*teps_a[i,2]<=phiepsk_a[i,2];
s.t. Qk2lepsbis{i in ind_a}:phiepsk_a[i,2]<=aepsk[i,2];
#constraints on each (b,f[i]+eps,alphaeps_b[i]) belonging to Q2
s.t. alp2epsb{i in ind_b}:alphaepsk_b[i,2]=(1/delta)*(phiepsk_b[i,2]-bepsk[i,2]+delta*teps_b[i,2]);
s.t. Qk2lepsb{i in ind_b}:bepsk[i,2]-delta*teps_b[i,2]<=phiepsk_b[i,2];
s.t. Qk2lepsbisb{i in ind_b}:phiepsk_b[i,2]<=bepsk[i,2];

#constraints on each (a,f[i],alpha[i]) belonging to Q3
s.t. Qk3{i in ind_a}:ak[i,3]<=phik_a[i,3];
s.t. alp3{i in ind_a}:alphak_a[i,3] = t_a[i,3];
#constraints on each (b,f[i],alpha_b[i]) belonging to Q3
s.t. Qk3b{i in ind_b}:bk[i,3]<=phik_b[i,3];
s.t. alp3b{i in ind_b}:alphak_b[i,3] = t_b[i,3];

#constraints on each (a,f[i]+eps,alphaeps[i]) belonging to Q3
s.t. Qk3eps{i in ind_a}:aepsk[i,3]<=phiepsk_a[i,3];
s.t. alp3eps{i in ind_a}:alphaepsk_a[i,3] = teps_a[i,3];
#constraints on each (a,f[i]+eps,alphaeps[i]) belonging to Q3
s.t. Qk3epsb{i in ind_b}:bepsk[i,3]<=phiepsk_b[i,3];
s.t. alp3epsb{i in ind_b}:alphaepsk_b[i,3] = teps_b[i,3];

