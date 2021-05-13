#######################
`aplc-optimization`
######################


Due to the linear properties, solutions are found by using a linear programming (LP) solver. Gurobi solver;


encode the optimization problem in matrix form using Python,  the Gurobi solver is directly called from Python using the gurobipy package.

For a given set of design survey parameters, the toolkit formulates the optimization program as a linear program

writes a linear programming script, which relies on the Gurobi solver to determine the apodizer mask solution with
maximum off-acis transmission for a given set of design constraints (namely the raw contrast goal, dark zone extent and spectral bandwidth, telescope pupil, occulting
mask, IWA, OWA, and Lyot stop profile).

For a given APLC configuration, a discrete, algebraic propagation model enables us to exactly define
the optimization objectives and constraints

we code the algebraic model and design goals as a linear program in the Python programming language.
For each of the APLC experiments, we used the Gurbo solver98 to solve the linear program and obtain the mask solution.

Gurobi python package to implement the solver algorithm