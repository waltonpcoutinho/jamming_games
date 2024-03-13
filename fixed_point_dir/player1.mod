###############################################################
# AMPL model for the paper: Jamming Games for Fleets of Mobile 
# Vehicles
#
# Authors: Jörg Fliege and Walton Coutinho
# Date: 04/03/2024
# Version: 0.0
# Known bugs: Not tested
###############################################################

###############################################################
########################## Model sets #########################
###############################################################

# Set of agents in fleet 1
set FLEET1;
# Set of agents in fleet 2
set FLEET2;
# Set of all agents
set ALL_AGENTS := FLEET1 union FLEET2;

###############################################################
####################### Model parameters ######################
###############################################################

# Scalar parameters
param varrho; # Antenna parameter
param alpha; # Path-loss coefficient
param sigmasqr; # Background noise
param eps; # General purpose tolerance parameter
param reg; # Regularization parameter

# Vector parameters
# initial (x,y) position of agent i in fleet f
param x0{i in ALL_AGENTS};
param y0{i in ALL_AGENTS};
param maxpow{i in ALL_AGENTS}; # Max power of agent i in fleet f
param c{i in ALL_AGENTS}; # conversion factor (W/m)
param delta{i in ALL_AGENTS}; # Maximum distance between agent's i origin and final position

param x_hat{i in FLEET1}; # x position of agent i in the previous iteration
param y_hat{i in FLEET1}; # y position of agent i in the previous iteration

###############################################################
####################### Model variables #######################
###############################################################
# Power transmitted by agent i to agent j in the same fleet
var pow{i in ALL_AGENTS, j in ALL_AGENTS: ((i in FLEET1 and j in FLEET1) or (i in FLEET2 and j in FLEET2)) and i != j} >= 0; 
# Jamming power sent by agent i to agent j in the other fleet
var jam_pow{i in ALL_AGENTS, k in ALL_AGENTS: ((i in FLEET1 and k not in FLEET1) or (i in FLEET2 and k not in FLEET2)) and (k != i)} >= 0; 

# x,y position of agent i
var x{i in ALL_AGENTS}; 
var y{i in ALL_AGENTS};

# Auxiliary variables
# Distance between agents i and j
var distance{i in ALL_AGENTS, j in ALL_AGENTS: i != j} = sqrt((x[i] - x[j])^2 + (y[i] - y[j])^2);
# Distance between agent i and its origin
var dist_orig{i in ALL_AGENTS} = sqrt((x[i] - x0[i])^2 + (y[i] - y0[i])^2);
# Regularization term distance
var reg_dist{i in FLEET1} = sqrt((x[i] - x_hat[i])^2 + (y[i] - y_hat[i])^2);

###############################################################
####################### Model formulation #####################
###############################################################

# Objective function
minimize obj: sum{i in FLEET1, j in FLEET1: i != j}
               (varrho * pow[i, j] * (distance[i, j]^(-alpha)))/
               (sigmasqr + varrho * sum{k in ALL_AGENTS: k not in FLEET1} jam_pow[k, j] * (distance[k, j]^(-alpha)))
               - sum{k in ALL_AGENTS, l in ALL_AGENTS: k not in FLEET1 and l not in FLEET1 and k != l}
                  (varrho * pow[k, l] * (distance[k, l]^(-alpha)))/
                  (sigmasqr + varrho * sum{i in FLEET1} jam_pow[i, l] * (distance[i, l]^(-alpha)))
               + reg * sum{i in FLEET1} reg_dist[i];

# Constraints
subject to

# Power constraint
power_constraint{i in FLEET1}:
      sum{j in FLEET1: j != i} pow[i, j] + sum{k in ALL_AGENTS: k not in FLEET1} jam_pow[i, k] + c[i] * dist_orig[i] <= maxpow[i];

# Max distance constraint
max_distance_constraint{i in FLEET1}:
      dist_orig[i] <= delta[i];






