The professor who wants to run the algorithm called Basic in the paper

    S. Burer, V. Piccialli, Three Methods for Robust Grading

can do it using the following ampl files:

    basic_grades_ab.mod = defines MILP model (no changes required by
    professor)

    example_grades_ab.dat = defines data of model (professor must insert
    the number of students and the grades of his class)

    commands_grades_ab.run = defines commands to run (no changes
    required by professor)

It is possible to solve the problem by submitting the files to the webpage

    https://neos-server.org/neos/solvers/milp:Gurobi/AMPL.html

where the solver Gurobi is used. (Any other solver for MILP is also
fine.) The file basic_grades_ab.mod has to be submitted in the field
Model File, while the file example_grades_ab.dat has to be submitted
in the field Data File and the file commands_grades_ab.run has to be
submitted in the field Commands File. Then push the button Submit to
Neos.

The server will process the file and print the answer on the webpage,
with an output as shown below. In particular, the final lines contain
the optimal solutions and optimal values.

*************************************************************

   NEOS Server Version 5.0
   Job#     : 5984368
   Password : tCJZVBWX
   User     : None
   Solver   : milp:Gurobi:AMPL
   Start    : 2018-04-09 08:35:23
   End      : 2018-04-09 08:35:43
   Host     : NEOS HTCondor Pool

   Disclaimer:

   This information is provided without any express or
   implied warranty. In particular, there is no warranty
   of any kind concerning the fitness of this
   information  for any particular purpose.
*************************************************************
File exists
You are using the solver gurobi_ampl.
Checking ampl.mod for gurobi_options...
Checking ampl.com for gurobi_options...
Executing AMPL.
processing data.
processing commands.
Executing on prod-exec-4.neos-server.org
set ind_a := 5 11 13 19 20 27 29 35 38;

set ind_b := 18 21 31;


Presolve eliminates 129 constraints and 79 variables.
Adjusted problem:
235 variables:
	57 binary variables
	178 linear variables
423 constraints, all linear; 1017 nonzeros
	119 equality constraints
	304 inequality constraints
1 linear objective; 19 nonzeros.

Gurobi 7.5.1: threads=4
outlev=1
Optimize a model with 423 rows, 235 columns and 1017 nonzeros
Variable types: 178 continuous, 57 integer (57 binary)
Coefficient statistics:
  Matrix range     [1e-03, 1e+03]
  Objective range  [1e+00, 1e+00]
  Bounds range     [1e+00, 9e+01]
  RHS range        [1e+00, 9e+01]
Presolve removed 339 rows and 164 columns
Presolve time: 0.01s
Presolved: 84 rows, 71 columns, 277 nonzeros
Variable types: 43 continuous, 28 integer (28 binary)
Found heuristic solution: objective 3.0000000
Found heuristic solution: objective 2.0000000

Root relaxation: objective -2.299925e+00, 45 iterations, 0.00 seconds

    Nodes    |    Current Node    |     Objective Bounds      |     Work
 Expl Unexpl |  Obj  Depth IntInf | Incumbent    BestBd   Gap | It/Node Time

     0     0   -2.29993    0    6    2.00000   -2.29993   215%     -    0s
H    0     0                       0.0000000   -2.29993      -     -    0s
     0     0   -1.75008    0    9    0.00000   -1.75008      -     -    0s
     0     0   -0.73645    0    4    0.00000   -0.73645      -     -    0s
     0     0   -0.12560    0    8    0.00000   -0.12560      -     -    0s
     0     0     cutoff    0         0.00000    0.00000  0.00%     -    0s

Explored 1 nodes (86 simplex iterations) in 0.04 seconds
Thread count was 4 (of 24 available processors)

Solution count 3: 0 2 3 

Optimal solution found (tolerance 1.00e-04)
Best objective 0.000000000000e+00, best bound 0.000000000000e+00, gap 0.0000%
Optimize a model with 423 rows, 235 columns and 1017 nonzeros
Coefficient statistics:
  Matrix range     [1e-03, 1e+03]
  Objective range  [1e+00, 1e+00]
  Bounds range     [1e+00, 9e+01]
  RHS range        [1e+00, 9e+01]
Iteration    Objective       Primal Inf.    Dual Inf.      Time
       0      handle free variables                          0s
     117    0.0000000e+00   0.000000e+00   0.000000e+00      0s

Solved in 117 iterations and 0.00 seconds
Optimal objective  0.000000000e+00
Gurobi 7.5.1: optimal solution; objective 0
86 simplex iterations
1 branch-and-cut nodes
plus 117 simplex iterations for intbasis
The optimal breakpoints are a = 89.101000 and b = 79.901000 
The total number of pseudoborderline students w.r.t. a and b are 0.000000
The total number of pseudoborderline students w.r.t. a  are 0.000000
The total number of pseudoborderline students w.r.t.  are 0.000000
