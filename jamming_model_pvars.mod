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

# Vector parameters
# initial (x,y) position of agent i in fleet f
param x0{i in ALL_AGENTS};
param y0{i in ALL_AGENTS};
param maxpow{i in ALL_AGENTS}; # Max power of agent i in fleet f
param c{i in ALL_AGENTS}; # conversion factor (W/m)

###############################################################
####################### Model variables #######################
###############################################################

# Power sent from agent i in F to agent j in F
var pow{t in TEAMS, i in FLEET[t], j in FLEET[t]} >= 0;
# Jamming power of agent i in F to jam agent j not in F
var pow_J{i in ALL_AGENTS, j in ALL_AGENTS: j != i} >= 0;

# (x,y) position of agent i in fleet f
var x{i in ALL_AGENTS};
var y{i in ALL_AGENTS};

#Lagrange multipliers
var lambda{i in ALL_AGENTS};
var mu{i in ALL_AGENTS, j in ALL_AGENTS: j != i} >= 0;

###############################################################
######################### Aux variables #######################
###############################################################

# Distance between agent's i origin and final position
var dist0{i in ALL_AGENTS} = sqrt((x[i] - x0[i])^2 + (y[i] - y0[i])^2);

# Distance between agent i in F and agent j in F
var dist{i in ALL_AGENTS, j in ALL_AGENTS: j != i} = sqrt((x[i] - x[j])^2 + (y[i] - y[j])^2);

# Derivatives
var partial_f1_pow{i in ALL_AGENTS, j in ALL_AGENTS: j != i};
var partial_f2_pow_J{i in ALL_AGENTS, j in ALL_AGENTS: j != i};
var partial_f1_x{i in ALL_AGENTS};
var partial_f2_x{i in ALL_AGENTS};

###############################################################
###################### Model formulation ######################

# Objective function
minimize Objective:
   0
;

# Constraints
subject to

# Aux variables
partial_f1_partial_p{t in TEAMS, i in FLEET[t], j in FLEET[t]: i !=j}:
   partial_f1_pow[i,j] = (varrho * dist[i, j]^(-alpha))/(sigmasqr + varrho * sum{k in ALL_AGENTS: k not in FLEET[t]} pow_J[k, j] * dist[k, j]^(-alpha));

partial_f2_partial_p{t in TEAMS, i in FLEET[t], k in ALL_AGENTS: k not in FLEET[t]}:
   partial_f2_pow_J[i,k] = sum{t2 in TEAMS, l in FLEET[t2]: t2 != t and l != k}(
      (varrho * pow[t2, l, k] * (dist[l, k]^(-alpha)) * (dist[i, k]^(-alpha)))/
      (sigmasqr + varrho * sum{j in FLEET[t]: j != i} (pow_J[j, k] * dist[j, k]^(-alpha)))^2
   );

partial_f1_x_def{t in TEAMS, i in FLEET[t]}:
   partial_f1_x[i] = sum{j in FLEET[t]: j != i}(
      (-varrho * pow[t, i, j] * alpha * (dist[i,j])^(-alpha - 2))/(sigmasqr + varrho * sum{k in ALL_AGENTS: k not in FLEET[t]} pow_J[k, j] * dist[k, j]^(-alpha))
   )*(x[i] - x[j])
   +  sum{j in FLEET[t]: j != i}(
      (-varrho * pow[t, j, i] * alpha * (dist[j,i])^(-alpha - 2))/(sigmasqr + varrho * sum{k in ALL_AGENTS: k not in FLEET[t]} pow_J[k, i] * dist[k, i]^(-alpha))
   )*(x[i] - x[j])
   + sum{j in FLEET[t]: j != i}(
      (-varrho * pow[t, j, i] * (dist[j,i])^(-alpha) * (-alpha * varrho * sum{k in ALL_AGENTS: k not in FLEET[t]} pow_J[k, i] * dist[k, i]^(-alpha - 2)*(x[i] - x[k])))
      /
      (sigmasqr + varrho * sum{k in ALL_AGENTS: k not in FLEET[t]} pow_J[k, i] * dist[k, i]^(-alpha))^2
   )
   ;

partial_f2_x_def{t in TEAMS, i in FLEET[t]}:
   partial_f2_x[i] = sum{t2 in TEAMS, k in FLEET[t2], l in FLEET[t2]: t2 != t and k != l}(
      ((alpha * varrho^2 * pow[t2, k, l] * pow_J[i, l] * (dist[k,l])^(-alpha) * (dist[i, l])^(-alpha - 2))/
      (sigmasqr + varrho * sum{j in FLEET[t]: j != i} pow_J[j, l] * dist[j, l]^(-alpha))^2)*(x[i] - x[l])
   );

# Lagrange critical points
partial_L_p{t in TEAMS, i in FLEET[t], j in FLEET[t]: j != i}:
   partial_f1_pow[i,j] + lambda[i] - mu[i,j] = 0;

partial_L_pJ{t in TEAMS, i in FLEET[t], k in ALL_AGENTS: k not in FLEET[t]}:
   partial_f2_pow_J[i,k] + lambda[i] - mu[i,k] = 0;

partial_L_x{t in TEAMS, i in FLEET[t]}:
   partial_f1_x[i] - partial_f2_x[i] 
   + ((c[i]*lambda[i])/(dist0[i]))*(x[i] - x0[i]) = 0;

# Feasibility conditions
power_feasibility{t in TEAMS, i in FLEET[t]}:
   sum{j in FLEET[t]: j != i} pow[t, i, j]
      + sum{k in ALL_AGENTS: k not in FLEET[t]} pow_J[i, k]
      + c[i] * (x[i] - x0[i]) - maxpow[i] == 0;

# Complementarity conditions
complementarity1{t in TEAMS, i in FLEET[t]}:
   lambda[i] complements (
      sum{j in FLEET[t]: j != i} pow[t, i, j]
      + sum{k in ALL_AGENTS: k not in FLEET[t]} pow_J[i, k]
      + c[i] * (x[i] - x0[i]) - maxpow[i]
   ) == 0;

complementarity2{t in TEAMS, i in FLEET[t], j in FLEET[t]: j != i}:
   0 <= mu[i,j] complements pow[t, i, j] >= 0;

complementarity3{t in TEAMS, i in FLEET[t], k in ALL_AGENTS: k not in FLEET[t]}:
   0 <= mu[i,k] complements pow_J[i,k] >= 0;

###############################################################
######################## Solve model ##########################
###############################################################

#load data
data toy_problem.dat;

#option solver pathampl;
option solver knitro;
option presolve 0;
solve;

###############################################################
########################## Display results ####################
###############################################################

display ALL_AGENTS;

display x0, y0, maxpow, c;

display x, y, pow, pow_J;

display dist0, dist;

display lambda;

display mu;
