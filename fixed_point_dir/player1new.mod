###############################################################
# AMPL model for the paper: Jamming Games for Fleets of Mobile 
# Vehicles
#
# Authors: Jörg Fliege and Walton Coutinho
# Date: 13/03/2024
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
# set ALL_AGENTS := FLEET1 union FLEET2;  # We should probably never use this.

###############################################################
####################### Model parameters ######################
###############################################################

# Scalar parameters for physics and numerics:
# Can be read from constant file.
param varrho; # Antenna parameter
param alpha; # Path-loss coefficient
param sigmasqr; # Background noise
param eps; # General purpose tolerance parameter
param reg; # Regularization parameter

# Vector parameters for fleet 1:
# Can be read from constant file, but we need TWO files, one for player 1 and one for player 2.
# initial (x,y) position of agent i in fleet 1
param x0{i in FLEET1};
param y0{i in FLEET1};
param maxpow{i in FLEET1}; # Max power of agent i in fleet 1
param c{i in FLEET1}; # conversion factor (W/m)
param delta{i in FLEET1}; # Maximum distance between agent's i origin and final position

# The following are read from a file we have written before in the last iteration:
# param x_hat{i in FLEET1}; # x position of agent i in the previous iteration 
# param y_hat{i in FLEET1}; # y position of agent i in the previous iteration 

# Parameters that are variables for player 2:
# These area read from some file that we just written.
param x_p2{k in FLEET2}; # x position of agent k in FLEET 2
param y_p2{k in FLEET2}; # y position of agent k in FLEET 2
param pow_p2{k in FLEET2, ell in FLEET2: k != ell} >= 0;
param jam_pow_p2{k in FLEET2, j in FLEET1} >= 0;
param distance_p2{k in FLEET2, l in FLEET2: k != l} >= 0;

###############################################################
####################### Model variables #######################
###############################################################
# Power transmitted by agent i to agent j in the same fleet
var pow{i in FLEET1, j in FLEET1:  i != j} >= 0; 
# Jamming power sent by agent i to agent j in the other fleet
var jam_pow{i in FLEET1, k in FLEET2} >= 0; 

# x,y position of agent i
var x{i in FLEET1}; 
var y{i in FLEET1};

# Auxiliary variables
# Distance between agents i and j, and i and k
var distance1{i in FLEET1, j in FLEET1: i != j} = sqrt((x[i] - x[j])^2 + (y[i] - y[j])^2);
var distance2{i in FLEET1, k in FLEET2} = sqrt((x[i] - x_p2[k])^2 + (y[i] - y_p2[k])^2);
# Distance between agent i and its origin
var dist_orig{i in FLEET1} = sqrt((x[i] - x0[i])^2 + (y[i] - y0[i])^2);
# Regularization term distance
# var reg_dist{i in FLEET1} = sqrt((x[i] - x_hat[i])^2 + (y[i] - y_hat[i])^2);

###############################################################
####################### Model formulation #####################
###############################################################

# Objective function
minimize obj: - sum{i in FLEET1, j in FLEET1: i != j}
               (varrho * pow[i, j] * (distance1[i, j]^(-alpha)))/
               (sigmasqr + varrho * sum{k in FLEET2} jam_pow_p2[k, j] * (distance2[j, k]^(-alpha)))
               + sum{k in FLEET2, l in FLEET2:  k != l}
                  (varrho * pow_p2[k, l] * (distance_p2[k, l]^(-alpha)))/
                  (sigmasqr + varrho * sum{i in FLEET1} jam_pow[i, l] * (distance2[i, l]^(-alpha)));
              # + reg * sum{i in FLEET1} (reg_dist[i])^2;

# Constraints
subject to

# Power constraint
power_constraint{i in FLEET1}:
      sum{j in FLEET1: j != i} pow[i, j] + sum{k in FLEET2} jam_pow[i, k] + c[i] * dist_orig[i] <= maxpow[i];

# Max distance constraint
max_distance_constraint{i in FLEET1}:
      dist_orig[i] <= delta[i];





