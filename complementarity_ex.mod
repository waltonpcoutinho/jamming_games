'''
Complementarity constraints example

A complementarity constraint enforces that two variables are complementary to each other; i.e., that the following conditions hold for scalar variables x and y: 

x >= 0, y >= 0, x*y = 0.

The condition above is sometimes expressed more compactly as:

0 \le x \perp y \ge 0.

Complementarity constraints are used to model a variety of problems, including variational inequalities, mathematical programs with equilibrium constraints, and bilevel programs.

source: https://www.artelys.com/app/docs/knitro/2_userGuide/complementarity.html
'''

# Reset model
reset;

# Variables
var x{j in 0..7} >= 0;

# Objective function
minimize obj:
                (x[0]-5)^2 + (2*x[1]+1)^2;

# Constraints
s.t. c0: 2*(x[1]-1) - 1.5*x[0] + x[2] - 0.5*x[3] + x[4] = 0;
s.t. c1: 3*x[0] - x[1] - 3 - x[5] = 0;
s.t. c2: -x[0] + 0.5*x[1] + 4 - x[6] = 0;
s.t. c3: -x[0] - x[1] + 7 - x[7] = 0;
s.t. c4: 0 <= x[5] complements x[2] >= 0;
s.t. c5: 0 <= x[6] complements x[3] >= 0;
s.t. c6: 0 <= x[7] complements x[4] >= 0;

# Select solver
option solver knitro;

# Solver options
option knitro_options 'presolve 0';

# Solve
solve;

# Display solution
display obj;
display x;