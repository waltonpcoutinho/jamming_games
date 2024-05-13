###############################################################
# AMPL model for the paper: Jamming Games for Fleets of Mobile 
# Vehicles
#
# Authors: Jörg Fliege and Walton Coutinho
# Date: 06/05/2024
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
# Temporary set to swap the values of FLEET1 and FLEET2
set TEMP_FLEET;
# Set of time indices
set Time;

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
# Big M parameter
param M := 100;

# Vector parameters for fleet 1:
# Can be read from constant file, but we need TWO files, one for player 1 and one for player 2.
# Initial (x,y) position of agent i in fleet 1
param x0{i in FLEET1};
param y0{i in FLEET1};
param z0{i in FLEET1};
param maxpow{i in FLEET1}; # Max power of agent i in fleet 1
param c{i in FLEET1}; # Conversion factor (W/m)
param delta{i in FLEET1}; # Maximum distance between agent's i origin and final position

# Communication network of Fleet 1
param comms_matrix_p1{i in FLEET1, j in FLEET1} >= 0 <= 1;
param jam_matrix_p1{i in FLEET1, l in FLEET2} >= 0 <= 1;

# The following are read from a file we have written in the previous iteration:
param x_hat{i in FLEET1, t in Time}; # x position of agent i in the previous iteration 
param y_hat{i in FLEET1, t in Time}; # y position of agent i in the previous iteration
param z_hat{i in FLEET1, t in Time}; # z position of agent i in the previous iteration

# Communication network of Fleet 2
param comms_matrix_p2{k in FLEET2, l in FLEET2} >= 0 <= 1;
param jam_matrix_p2{k in FLEET2, j in FLEET1} >= 0 <= 1;

# Parameters that are variables for player 2:
# These area read from some file that we just written.
param x_p2{k in FLEET2, t in Time}; # x position of agent k in FLEET 2
param y_p2{k in FLEET2, t in Time}; # y position of agent k in FLEET 2
param z_p2{k in FLEET2, t in Time}; # z position of agent k in FLEET 2
param pow_p2{k in FLEET2, ell in FLEET2, t in Time: k != ell} >= 0 <= M * comms_matrix_p2[k, ell];
param jam_pow_p2{k in FLEET2, i in FLEET1, j in FLEET1, t in Time: i != j} >= 0 <= M * jam_matrix_p2[k, j];
param distance_p2{k in FLEET2, l in FLEET2, t in Time: k != l} >= 0;

###############################################################
####################### Model variables #######################
###############################################################
# Power transmitted by agent i to agent j in the same fleet
var pow{i in FLEET1, j in FLEET1, t in Time:  i != j} >= 0 <= M * comms_matrix_p1[i, j];
# Jamming power sent by agent i to agent j in the other fleet
var jam_pow{i in FLEET1, k in FLEET2, l in FLEET2, t in Time: k != l} >= 0 <= M * jam_matrix_p1[i, l];

# x,y position of agent i
var x{i in FLEET1, t in Time}; 
var y{i in FLEET1, t in Time};
var z{i in FLEET1, t in Time};

# Auxiliary variables
# Distance between agents i and j, and i and k
var distance1{i in FLEET1, j in FLEET1, t in Time: i != j} = sqrt((x[i,t] - x[j,t])^2 + (y[i,t] - y[j,t])^2 + (z[i,t] - z[j,t])^2);
var distance2{i in FLEET1, k in FLEET2, t in Time} = sqrt((x[i,t] - x_p2[k,t])^2 + (y[i,t] - y_p2[k,t])^2 + (z[i,t] - z_p2[k,t])^2);
# Distance between agent i and its origin
var dist_travel{i in FLEET1, t in Time: t > 0} = sqrt((x[i,t] - x[i,t-1])^2 + (y[i,t] - y[i,t-1])^2 + (z[i,t] - z[i,t-1])^2);
# Regularization term distance
var reg_dist{i in FLEET1, t in Time} = sqrt((x[i,t] - x_hat[i,t])^2 + (y[i,t] - y_hat[i,t])^2 + (z[i,t] - z_hat[i,t])^2);

###############################################################
####################### Model formulation #####################
###############################################################

# Objective function
minimize obj: - sum{t in Time, i in FLEET1, j in FLEET1: i != j}
               (varrho * comms_matrix_p1[i, j] * pow[i, j, t] * (distance1[i, j, t]^(-alpha)))/
               (sigmasqr + varrho * sum{k in FLEET2} jam_matrix_p2[k, j] * jam_pow_p2[k, i, j, t] * (distance2[j, k, t]^(-alpha)))
               + sum{t in Time, k in FLEET2, l in FLEET2:  k != l}
                  (varrho * comms_matrix_p2[k, l] * pow_p2[k, l, t] * (distance_p2[k, l, t]^(-alpha)))/
                  (sigmasqr + varrho * sum{i in FLEET1} jam_matrix_p1[i, l] * jam_pow[i, k, l, t] * (distance2[i, l, t]^(-alpha)))
              + reg * sum{t in Time, i in FLEET1} (reg_dist[i, t])^2;

# Constraints
subject to

# Power constraint
power_constraint{i in FLEET1}:
      sum{t in Time, j in FLEET1: j != i} comms_matrix_p1[i, j] * pow[i, j, t] + sum{t in Time, k in FLEET2, l in FLEET2: k != l} jam_matrix_p1[i, l] * jam_pow[i, k, l, t] + sum{t in Time: t > 0} c[i] * dist_travel[i, t] <= maxpow[i];

# Max vel constraint
max_vel_1{i in FLEET1, t in Time: t > 0}:
      eps <= dist_travel[i, t] <= delta[i];
# Velocity at 0
max_vel_0{i in FLEET1}: 
      eps <= sqrt((x[i, 1] - x[i, 0])^2 + (y[i, 1] - y[i, 0])^2 + (z[i, 1] - z[i, 0])^2) <= delta[i];

# Initial position constraint
initial_x{i in FLEET1}:
      x[i, 0] - x0[i] == 0;
initial_y{i in FLEET1}:
      y[i, 0] - y0[i] == 0;
initial_z{i in FLEET1}:
      z[i, 0] - z0[i] == 0;

# A few anti-degeneracy constraints (necessary for most of the tested solvers)
anti1{i in FLEET1, j in FLEET1, t in Time: i != j}:
      eps <= distance1[i, j, t];
anti2{i in FLEET1, k in FLEET2, t in Time}:
      eps <= distance2[i, k, t];
anti3{i in FLEET1, t in Time}:
      eps <= reg_dist[i, t];

