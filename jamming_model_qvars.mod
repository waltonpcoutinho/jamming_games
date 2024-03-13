###############################################################
# AMPL model for the paper: Jamming Games for Fleets of Mobile 
# Vehicles
#
# Authors: Jörg Fliege and Walton Coutinho
# Date: 15/02/2024
# Version: 0.0
# Known bugs: Not tested
###############################################################

###############################################################
########################## Model sets #########################
###############################################################

# Set TEAMS is a set of fleets of agents
set TEAMS;
# Set of agents in each fleet
set FLEET{TEAMS};
# Set of all agents is the union of all FLEETs
set ALL_AGENTS := union {i in TEAMS} FLEET[i];

###############################################################
####################### Model parameters ######################
###############################################################

# Scalar parameters
param varrho; # Antenna parameter
param alpha; # Path-loss coefficient
param sigmasqr; # Background noise
param eps; # General purpose tolerance parameter

# Vector parameters
# initial (x,y) position of agent i in fleet f
param x0{i in ALL_AGENTS};
param y0{i in ALL_AGENTS};
param maxpow{i in ALL_AGENTS}; # Max power of agent i in fleet f
param c{i in ALL_AGENTS}; # conversion factor (W/m)
param delta{i in ALL_AGENTS}; # Maximum distance between agent's i origin and final position

###############################################################
####################### Model variables #######################
###############################################################

# Power sent from agent i in F to agent j in F
var q{t in TEAMS, i in FLEET[t], j in FLEET[t]} >= 0;
# Jamming power of agent i in F to jam agent j not in F
var q_J{i in ALL_AGENTS, j in ALL_AGENTS: j != i} >= 0;

# (x,y) position of agent i in fleet f
var x{i in ALL_AGENTS};
var y{i in ALL_AGENTS};

# Auxiliary variables
# Distance between agent's i origin and final position
var dist0{i in ALL_AGENTS} = sqrt((x[i] - x0[i])^2 + (y[i] - y0[i])^2);
# Distance between agent i in F and agent j in F
var dist{i in ALL_AGENTS, j in ALL_AGENTS: j != i} = sqrt((x[i] - x[j])^2 + (y[i] - y[j])^2);

#Lagrange multipliers
var lambda{i in ALL_AGENTS};
var mu{i in ALL_AGENTS, j in ALL_AGENTS: j != i} >= 0;
var mu_dist{i in ALL_AGENTS} >= 0;
var nu{i in ALL_AGENTS, j in ALL_AGENTS: j != i} >= 0;
var xcsi{i in ALL_AGENTS} >= 0;

###############################################################
######################### Aux Expressions #####################
###############################################################

var power{t in TEAMS, i in FLEET[t], j in FLEET[t]: i != j} = q[t, i, j] * dist[i, j]^alpha;

var power_J{t in TEAMS, i in FLEET[t], k in ALL_AGENTS: k not in FLEET[t]} = q_J[i, k] * dist[i, k]^alpha;

var partial_intra_L_q{t in TEAMS, i in FLEET[t], j in FLEET[t]: i != j} =
   ((-1)/((sigmasqr/varrho) + sum{k in ALL_AGENTS: k not in FLEET[t]} q_J[k, j])) + lambda[i] * (dist[i, j]^alpha) - nu[i, j];

var partial_inter_L_q{t in TEAMS, i in FLEET[t], t2 in TEAMS, l in FLEET[t2]: t2 != t} = 
   sum{k in FLEET[t2]: k != l} ((- q[t2, k, l])/((sigmasqr/varrho + sum{j in FLEET[t]} q_J[j, l])^2)) + lambda[i] * (dist[i, l]^alpha) - nu[i, l];

var partial_intra_L_t{t in TEAMS, i in FLEET[t], j in FLEET[t]: i != j} = 
   alpha * lambda[i] * q[t, i, j] * (dist[i, j]^(alpha - 1)) - mu[i, j];

var partial_inter_L_t{t in TEAMS, i in FLEET[t], k in ALL_AGENTS: k not in FLEET[t]} =
   alpha * lambda[i] * q_J[i, k] * (dist[i, k]^(alpha - 1)) - mu[i, k];

var partial_L_dist0{i in ALL_AGENTS} = c[i] * lambda[i] - mu_dist[i];

var partial_L_x{t in TEAMS, i in FLEET[t]} = 
   sum{j in FLEET[t]: j != i} ((mu[i, j] + mu[j, i])/(sqrt((x[i] - x[j])^2 + (y[i] - y[j])^2)))*(x[i] - x[j])
   + sum{k in ALL_AGENTS: k not in FLEET[t]} (mu[i, k]/(sqrt((x[i] - x[k])^2 + (y[i] - y[k])^2)))*(x[i] - x[k])
   + (mu_dist[i] + xcsi[i]) * (x[i] - x0[i])/(sqrt((x[i] - x0[i])^2 + (y[i] - y0[i])^2));

var partial_L_y{t in TEAMS, i in FLEET[t]} = 
   sum{j in FLEET[t]: j != i} ((mu[i, j] + mu[j, i])/(sqrt((x[i] - x[j])^2 + (y[i] - y[j])^2)))*(y[i] - y[j])
   + sum{k in ALL_AGENTS: k not in FLEET[t]} (mu[i, k]/(sqrt((x[i] - x[k])^2 + (y[i] - y[k])^2)))*(y[i] - y[k])
   + (mu_dist[i] + xcsi[i]) * (y[i] - y0[i])/(sqrt((x[i] - x0[i])^2 + (y[i] - y0[i])^2));

###############################################################
######################### KKT conditions ######################
###############################################################

subject to 

# Stationary points
KKT_power_stationary_1{t in TEAMS, i in FLEET[t], j in FLEET[t]: i != j}:
   partial_intra_L_q[t, i, j] == 0;
KKT_power_stationary_2{t in TEAMS, i in FLEET[t], t2 in TEAMS, l in FLEET[t2]: t2 != t}:
   partial_inter_L_q[t, i, t2, l] == 0;

KKT_dist_stationary1{t in TEAMS, i in FLEET[t], j in FLEET[t]: i != j}:
   partial_intra_L_t[t, i, j] == 0;
KKT_dist_stationary2{t in TEAMS, i in FLEET[t], k in ALL_AGENTS: k not in FLEET[t]}:
   partial_inter_L_t[t, i, k] == 0;
KKT_dist_stationary3{i in ALL_AGENTS}:
   partial_L_dist0[i] == 0;

KKT_x_stationary{t in TEAMS, i in FLEET[t]}:
   partial_L_x[t,i] == 0;
KKT_y_stationary{t in TEAMS, i in FLEET[t]}:
   partial_L_y[t,i] == 0;

# Feasibility constraints
# Power constraints
feasibility_power{t in TEAMS, i in FLEET[t]}:
   sum{j in FLEET[t]: j != i} q[t, i, j] * (dist[i, j]^alpha) 
   + sum{k in ALL_AGENTS: k not in FLEET[t]} q_J[i, k] * (dist[i, k]^alpha)
   + c[i] * dist0[i] - maxpow[i] == 0;

# Distance constraints
feasibility_intra_dist{t in TEAMS, i in FLEET[t], j in FLEET[t]: i != j}:
   sqrt((x[i] - x[j])^2 + (y[i] - y[j])^2) <= dist[i, j];
feasibility_inter_dist{t in TEAMS, i in FLEET[t], k in ALL_AGENTS: k not in FLEET[t]}:
   sqrt((x[i] - x[k])^2 + (y[i] - y[k])^2) <= dist[i, k];
feasibility_dist{i in ALL_AGENTS}:
   sqrt((x[i] - x0[i])^2 + (y[i] - y0[i])^2) <= dist0[i];
feasibility_vel{i in ALL_AGENTS}:
   sqrt((x[i] - x0[i])^2 + (y[i] - y0[i])^2) <= delta[i];

# Complemetarity constraints
complementarity_power{t in TEAMS, i in FLEET[t]}:
   lambda[i] complements (sum{j in FLEET[t]: j != i} q[t, i, j] * (dist[i, j]^alpha) 
               + sum{k in ALL_AGENTS: k not in FLEET[t]} q_J[i, k] * (dist[i, k]^alpha)
               + c[i] * dist0[i] - maxpow[i]) == 0;

complementarity_dist1{t in TEAMS, i in FLEET[t], j in FLEET[t]: i != j}:
   0 <= mu[i, j] complements (sqrt((x[i] - x[j])^2 + (y[i] - y[j])^2) - dist[i, j]) <= 0;
complementarity_dist2{t in TEAMS, i in FLEET[t], k in ALL_AGENTS: k not in FLEET[t]}:
   0 <= mu[i, k] complements (sqrt((x[i] - x[k])^2 + (y[i] - y[k])^2) - dist[i, k]) <= 0;
complementarity_dist3{i in ALL_AGENTS}:
   0 <= mu_dist[i] complements (sqrt((x[i] - x0[i])^2 + (y[i] - y0[i])^2) - dist0[i]) <= 0;
complementarity_dist4{i in ALL_AGENTS}:
   0 <= xcsi[i] complements (sqrt((x[i] - x0[i])^2 + (y[i] - y0[i])^2) - delta[i]) <= 0;

# Non-negativity complementarity
non_negativity_complementarity{t in TEAMS, i in FLEET[t], j in FLEET[t]: i != j}:
   0 <= nu[i, j] complements q[t, i, j] >= 0;
non_negativity_complementarity2{t in TEAMS, i in FLEET[t], k in ALL_AGENTS: k not in FLEET[t]}:
   0 <= nu[i, k] complements q_J[i, k] >= 0;

###############################################################
######################## Solve model ##########################
###############################################################

#load data
data toy_problem.dat;

#option solver pathampl;
option presolve 0;
option solver knitro;
option knitro_options "ms_enable=1 ms_deterministic=0 ms_numthreads=0 ms_num_to_save=5 ms_savetol=0.01";
solve;

###############################################################
########################## Initial values #####################
###############################################################

let x['a'] := x0['a'] - eps;
let y['a'] := y0['a'] - eps;
let x['b'] := x0['b'] - eps;
let y['b'] := y0['b'] - eps;
let x['A'] := x0['A'] - eps;
let y['A'] := y0['A'] - eps;
let x['B'] := x0['B'] - eps;
let y['B'] := y0['B'] - eps;

for {t in TEAMS, i in FLEET[t], j in FLEET[t]: i !=j} 
   let q[t, i ,j] := 0.0;
;

for {t in TEAMS, i in FLEET[t], k in ALL_AGENTS: i != k} 
   let q_J[i ,k] := 0.0;
;

for {i in ALL_AGENTS}
   let lambda[i] := 0.0;
;

for {i in ALL_AGENTS}
   let mu_dist[i] := 0.0;
;

for {i in ALL_AGENTS, j in ALL_AGENTS: i != j}
   let mu[i, j] := 0.0;
;

for {i in ALL_AGENTS, j in ALL_AGENTS: i != j}
   let nu[i, j] := 0.0;
;

###############################################################
########################## Display results ####################
###############################################################

display ALL_AGENTS;

display x0, y0, maxpow, c;

display x, y, q, q_J;

display power, power_J;

#display feasibility_power;

for {t in TEAMS, i in FLEET[t]}
   display sum{j in FLEET[t]: j != i} q[t, i, j] * dist[i, j]^alpha 
      + sum{k in ALL_AGENTS: k not in FLEET[t]} q_J[i, k] * dist[i, k]^alpha
      + c[i] * dist0[i];
;

display dist0, dist;

display lambda;

display mu, mu_dist;

display nu;














