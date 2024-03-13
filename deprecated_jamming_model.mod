###############################################################
# AMPL model for the paper: Jamming Games for Fleets of Mobile 
# Vehicles
#
# Authors: Jörg Fliege and Walton Coutinho
# Date: 11/02/2024
# Version: 0.0
# Known bugs: Not tested
###############################################################

###############################################################
####################### Model parameters ######################
###############################################################

# Scalar parameters
param m; # Number of fleets of agents
param n_i; # Number of agents in each fleet
param varrho; # Antenna parameter
param alpha; # Path-loss coefficient
param sigmasqr; # Background noise

# Sets
# Set S is a set of fleets of agents
set S = 1..m;
# Set of agents in each fleet
set F{f in S} = 1..n_i;

# Vector parameters
# initial (x,y) position of agent i in fleet f
param x0{f in S, i in F[f]};
param y0{f in S, i in F[f]};
# Max power of agent i in fleet f
param maxpow{f in S, i in F[f]};
# Conversion factor from dist to power
param c{f in S, i in F[f]}; # dimension W/m

###############################################################
####################### Model variables #######################
###############################################################

# Variables
# Comms power received by agent j from agent i in fleet f, i != j
var pow{f in S, i in F[f], j in F[f]: i != j} >= 0;

# Jamming power received by agent i in f from agents k not in f
var pow_J{f in S, i in F[f], foes in S, k in F[foes]: foes != f} >= 0;

# (x,y) position of agent i in fleet f
var x{f in S, i in F[f]};
var y{f in S, i in F[f]};
# Lagrange multiplier for the power constraints
var lambda{f in S, i in F[f]} >= 0;
# Lagrange multiplier for the nonnegative power constraints
var mu_R{f in S, i in F[f], j in F[f]: i != j} >= 0; 
var mu_J{f in S, i in F[f], foes in S, k in F[foes]: foes != f} >= 0;

# Auxiliary variables
# Distance between agents i and j in fleet f
var intra_dist{f in S, i in F[f], j in F[f]: i != j} = sqrt((x[f,i] - x[f,j])^2 + (y[f,i] - y[f,j])^2)^(-alpha);

# Distance between agents i in fleet f and agents k not in f
var inter_dist{f in S, i in F[f], foes in S, k in F[foes]: foes != f} = sqrt((x[f,i] - x[foes,k])^2 + (y[f,i] - y[foes,k])^2)^(-alpha);

# distance between agents l and k not in fleet f
var ext_dist{f in S, foes1 in S, foes2 in S, k in F[foes1], l in F[foes2] : foes1 != f and foes2 != f and (l != k and foes1 == foes2)} = sqrt((x[foes1,k] - x[foes2,l])^2 + (y[foes1,k] - y[foes2,l])^2)^(-alpha);

# Distance from origin to agent i in fleet f
var orig_dist{f in S, i in F[f]} = sqrt((x[f,i] - x0[f,i])^2 + (y[f,i] - y0[f,i])^2);

###############################################################
####################### Model formulation #####################
###############################################################

# Objective function
minimize Objective:
   0
;

subject to

derivative_pow{f in S, i in F[f], j in F[f]: i != j}:
   (varrho * intra_dist[f,i,j]^(-alpha))/(sigmasqr + varrho * sum{foes in S, k in F[foes]: foes != f} pow_J[f, j, foes, k] * inter_dist[f, j, foes, k]^(-alpha)) + lambda[f,i] - mu_R[f,i,j] == 0;

derivative_pow_J{f in S, i in F[f], foes in S, k in F[foes]: foes != f}: sum{foes2 in S, l in F[foes2]: foes2 != f and l != k} (varrho * pow[foes, l, k] * (intra_dist[foes, l, k]^(-alpha)) * (inter_dist[f, i, foes, k]^(-alpha)))/(sigmasqr + varrho * (sum{j in F[f]: j != i} (pow_J[f, j, foes, k] * (inter_dist[f, j, foes, k]^(-alpha)))))^2 + lambda[f,i] - mu_J[f,i,foes,k] == 0;

complement_lambda{f in S, i in F[f]}:
   0 <= lambda[f,i] complements (sum{j in F[f]: j != i} pow[f, i, j] + sum{foes in S, k in F[foes]: foes != f} pow_J[f, i, foes, k] + c[f,i] * orig_dist[f,i] - maxpow[f,i]) >= 0;

complement_mu1{f in S, i in F[f], j in F[f]: i != j}:
   0 <= mu_R[f,i,j] complements pow[f,i,j] >= 0;

complement_mu2{f in S, i in F[f], foes in S, k in F[foes]: foes != f}:
   0 <= mu_J[f,i,foes,k] complements pow_J[f,i,foes,k] >= 0;

###############################################################
####################### Data declaration ######################
###############################################################

data;
param m := 2;
param n_i := 2;
param varrho := 1;
param alpha := 1;
param sigmasqr := 1;

param: x0 y0 :=
1 1 0 0
1 2 0 0
2 1 0 0
2 2 0 0
;

param: maxpow :=
1 1 1
1 2 1
2 1 1
2 2 1
;

param: c :=
1 1 1
1 2 1
2 1 1
2 2 1
;

###############################################################
######################## Solve model ##########################
###############################################################

#option solver pathampl;
option solver knitro;
option presolve 0;
solve;

###############################################################
####################### Display results #######################
###############################################################

display x0, y0, maxpow;

display x, y, pow, pow_J;

display intra_dist, inter_dist, ext_dist, orig_dist;

end;
```