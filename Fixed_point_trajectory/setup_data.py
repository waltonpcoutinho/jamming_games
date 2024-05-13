#!/usr/bin/env python3

'''
Python sript to write the required data for running the
AMPL's jamming.run script.

The input file (.dat) has the following format:

      Fleet #<fleet number>
      Agents <number of agents>
      //for each agent in Fleet
      agent_indx x0 y0 z0 maxpow c delta
      Fleet #<fleet number>
      Agents <number of agents>
      //for each agent in fleet
      agent_indx x0 y0 z0 maxpow c delta

Each output file has the following format:

Consider the example below where fleet 1 has to agents A and B.

The file is named "data_team1.dat" and describes the following:

- x0, y0, z0: initial position of the agents
- maxpow: maximum power available for each agent
- c: conversion factor from power to distance (cost to move)
- delta: distance between the agents

###############################################################
####################### Set parameters ########################
###############################################################

param: x0 y0 z0 maxpow c delta :=
A 1 1 1 1 0.1 1
B -1 1 1 1 0.1 1
;

The second type of output file writes the AMPL sets that describe 
each fleet. The file is named "fleet_team<i>.dat" and has the 
format:

Take the example below where fleet 1 has two agents A and B.

###############################################################
###################### Set fleet members ######################
###############################################################

# Define fleet 1
set FLEET1 := 'A' 'B' ;

The same is done for fleet 2.

The third type of output file writes the time index set for the
players. The file is named "time_set.dat" and has the format:

Example: consider N = 5 time steps.

###############################################################
####################### Set time index ########################
###############################################################
# Time set
set Time := 0 1 2 3 4;
'''

import sys
import math
import random
import re

# Class to write the data for the AMPL script
class DataWriter:
   
      def __init__(self, file_name):
         self.file_name = "Inst_Data/" + file_name
         self.fyle = open(self.file_name, "w")
   
      def write_header(self):
         self.fyle.write("###############################################################\n")
         self.fyle.write("####################### Set parameters ########################\n")
         self.fyle.write("###############################################################\n")
         self.fyle.write("param: x0 y0 z0 maxpow c delta :=\n")
   
      def write_data(self, agent, x0, y0, z0, maxpow, c, delta):
         self.fyle.write(agent + " " + str(x0) + " " + str(y0) + " " + str(z0) + " " + str(maxpow) + " " + str(c) + " " + str(delta) + "\n")

      # def write_matrix(self, matrix, agents, fleet):
      #    # Flatten the matrix
      #    vals = [el for row in matrix for el in row]
      #    # Write the matrix to file
      #    count = 0
      #    self.fyle.write("param comms_matrix_p%s := \n" % fleet)
      #    for ag1 in agents:
      #       for ag2 in agents:
      #          self.fyle.write(ag1 + " " + ag2 + " " + str(vals[count]) + "\n")
      #          count += 1

      def write_matrix(self, matrix, fleet):
         for row in matrix:
            self.fyle.write(" ".join([str(val) for val in row]) + "\n")
   
      def write_footer(self):
         self.fyle.write(";\n\n")
   
      def close(self):
         self.fyle.close()

# Class to write the fleet data for the AMPL script
class FleetWriter:
   def __init__(self, file_name):
      self.file_name = "Inst_Data/" + file_name
      self.fyle = open(self.file_name, "w")

   def write_fleet(self, fleet, agents):
      self.fyle.write("# Define fleet " + fleet + "\n")
      self.fyle.write("set FLEET" + fleet + " := ")
      for agent in agents:
         self.fyle.write("'" + agent + "' ")
      self.fyle.write(";\n")

   def close(self):
      self.fyle.close()

# Class for writing feasible initial guesses for the agents of
# fleets 1 and 2
class InitialGuessWriter:
   def __init__(self, file_name):
      self.file_name = "Manual_Data/" + file_name
      self.fyle = open(self.file_name, "w")

   def compute_init_guesses(self, fleet, agent, N, all_fleets, rand_walk, comms_arr, jams_arr):
      # Parameters
      fleet_size = len(all_fleets[fleet])
      maxpow = all_fleets[fleet][agent][4]
      c = all_fleets[fleet][agent][5]
      delta = all_fleets[fleet][agent][6]

      # Compute tavelled distances
      # travel_dist = \sum_{t in T\{0}} ||x_{it} - x_{i,t-1}||
      x = rand_walk[fleet][agent][0]
      y = rand_walk[fleet][agent][1]
      z = rand_walk[fleet][agent][2]
      travel_dist = sum([math.sqrt((x[t] - x[t-1])**2 + (y[t] - y[t-1])**2 + (z[t] - z[t-1])**2) for t in range(1, N)])

      # Power allocation
      # \sum_{j in fleet} pow_{ij} + \sum_{k not in fleet} jam_{ik} + c * dist <= maxpow
      remaining_pow = maxpow - c * travel_dist
      # Part of power goes for communication and part for jamming 
      #12/05/2024: Checked the numbers, it seems to work for >= 2 agents
      
      # Comms allocation
      pow = {}
      for j in all_fleets[fleet]:
         if j != agent:
            # Get index of agents from the dictionary
            ag_idx = list(all_fleets[fleet]).index(agent)
            j_idx = list(all_fleets[fleet]).index(j)
            # Comms matrix val
            comms_net = comms_arr[ag_idx][j_idx]
            if comms_net >= 1 - 1e-6:
               pow[j] = [(0.5*remaining_pow)/(fleet_size*N) for t in range(N)]
            else:
               pow[j] = [0 for t in range(N)]
      
      # Jamming allocation
      jam = {}
      for f2 in all_fleets:
         if f2 != fleet:
            f2_size = len(all_fleets[f2])
            for enmy1 in all_fleets[f2]:
               jam[enmy1] = {}
               for enmy2 in all_fleets[f2]:
                  if enmy1 != enmy2:
                     # Get index of agents
                     ag_idx = list(all_fleets[fleet]).index(agent)
                     enmy2_idx = list(all_fleets[f2]).index(enmy2)
                     # Jamming matrix val
                     jam_net = jams_arr[ag_idx][enmy2_idx]
                     if jam_net >= 1 - 1e-6:
                        jam[enmy1][enmy2] = [(0.5*remaining_pow)/((f2_size**2 - f2_size)*N) for t in range(N)]
                     else:
                        jam[enmy1][enmy2] = [0 for t in range(N)]

      # Compute the sum of power and jamming
      sum_pow = sum([sum(pow[i]) for i in pow])
      sum_jam = sum([sum(sum(jam[k][l]) for l in jam[k]) for k in jam])
      power_expenditure = sum_pow + sum_jam + c * travel_dist

      # Print some log
      print(agent, ": Comms= ", sum_pow, ", Jams= ", sum_jam, ", Move= ", c*travel_dist, ", Total= ", power_expenditure, ", Max= ", maxpow)

      #Assert that the power allocation is feasible
      assert power_expenditure <= maxpow, "Power expenditure exceeds the maximum power available for agent " + agent + " in fleet " + fleet + ".\nPower expenditure: " + str(power_expenditure) + " Max power: " + str(maxpow) + "\n"

      # Compute distances between agents
      dist = {}
      for j in all_fleets[fleet]:
         if j != agent:
            dist[j] = [math.sqrt((x[t] - rand_walk[fleet][j][0][t])**2 + (y[t] - rand_walk[fleet][j][1][t])**2 + (z[t] - rand_walk[fleet][j][2][t])**2) for t in range(N)]
      
      # Return the initial guesses
      return x, y, z, pow, jam, dist

   def write_initial_guess(self, fleet, agents, N, all_fleets, rand_walk, comms_arr, jams_arr):
      for ag in agents:
         # Compute the initial guesses for agent i
         x_guess, y_guess, z_guess, pow_guess, jam_guess, dist_guess = self.compute_init_guesses(fleet, ag, N, all_fleets, rand_walk, comms_arr, jams_arr)
         for t in range(N):
            self.fyle.write(str(x_guess[t]))
            self.fyle.write('\n')
            print(ag, t, x_guess[t])
         for t in range(N):
            self.fyle.write(str(y_guess[t]))
            self.fyle.write('\n')
            print(ag, t, y_guess[t])
         for t in range(N):
            self.fyle.write(str(z_guess[t]))
            self.fyle.write('\n')
            print(ag, t, z_guess[t])
         for j in agents:
            if j != ag:
               for t in range(N):
                  self.fyle.write(str(pow_guess[j][t]))
                  self.fyle.write('\n')
                  print(ag, j, t, pow_guess[j][t])
         for f2 in all_fleets:
            if f2 != fleet:
               for k in all_fleets[f2]:
                  for l in all_fleets[f2]:
                     if k != l:
                        for t in range(N):
                           self.fyle.write(str(jam_guess[k][l][t]))
                           self.fyle.write('\n')
                           print(ag, k, l, t, jam_guess[k][l][t])
         for j in agents:
            if ag != j:
               for t in range(N):
                  self.fyle.write(str(dist_guess[j][t]))
                  self.fyle.write('\n')
                  print(ag, j, t, dist_guess[j][t])

   def close(self):
      self.fyle.close()

# Function to extract the data from the input file
def extract_data(lines, N, user_maxpow, user_c, user_delta):
   # Each fleet has a dictonary with the agents and their data
   fleets = {}
   comms = {}
   jams = {}
   for line in lines:      
      if "Fleet" in line:
         fleet = line.split()[1][1]
         fleets[fleet] = {}
      elif "Agents" in line:
         agents =int(line.split()[1])
      elif "Comms" in line:
         # Extract index after "Comms_"
         fleet_idx = re.search(r'Comms_(\d+)', line).group(1)
         # Extract data inside square brackets
         data_str = re.search(r'\[(.*?)\]', line).group(1)
         # Parse data string into an array
         comms_arr = [[int(num) for num in row.split()] for row in data_str.split(';')]
         # Store the communication matrix
         comms[fleet_idx] = comms_arr
      elif "Jams" in line:
         # Extract index after "Jams_"
         fleet_idx = re.search(r'Jams_(\d+)', line).group(1)
         # Extract data inside square brackets
         data_str = re.search(r'\[(.*?)\]', line).group(1)
         # Parse data string into an array
         jams_arr = [[int(num) for num in row.split()] for row in data_str.split(';')]
         # Store the jamming matrix
         jams[fleet_idx] = jams_arr
      else:
         data = line.split()
         agent = data[0]
         x0 = float(data[1])
         y0 = float(data[2])
         z0 = float(data[3])
         maxpow = float(data[4])
         c = float(data[5])
         delta = float(data[6])
         # Check if the user has provided the maxpow, c, and delta
         if user_maxpow:
            maxpow = user_maxpow
         if user_c:
            c = user_c
         if user_delta:
            delta = user_delta
         # Add the agent to the fleet
         fleets[fleet][agent] = (agent, x0, y0, z0, maxpow, c, delta)
   
   # Return data structures
   return fleets, comms, jams
   
def main(file_name, N=5, *args, **kwargs):

   # First check for optional keyword arguments
   user_maxpow = kwargs.get('maxpow', None)
   user_c = kwargs.get('c', None)
   user_delta = kwargs.get('delta', None)

   # Read the input file
   with open(file_name, "r") as fyle:
      lines = fyle.readlines()

   # Extract the data from the input file
   fleets, comms, jams = extract_data(lines, N, user_maxpow, user_c, user_delta)

   # Write the fleets for the AMPL script
   for fleet in fleets:
      # Get the agents in the fleet
      agents = [agent for agent in fleets[fleet]]

      fleet_file = "fleet_team" + fleet + ".dat"
      fleet_writer = FleetWriter(fleet_file)
      fleet_writer.write_fleet(fleet, agents)
      fleet_writer.close()

      # Write the data for the AMPL script
      data_file = "data_team" + fleet + ".dat"
      data_writer = DataWriter(data_file)
      data_writer.write_header()
      for agent in fleets[fleet]:
         data_writer.write_data(*fleets[fleet][agent])
      data_writer.write_footer()
      data_writer.close()

      # Write the communication and jamming matrices
      comms_file = "comms_matrix_p" + fleet + ".dat"
      mtrx_writer = DataWriter(comms_file)
      mtrx_writer.write_matrix(comms[fleet], fleet)

      jams_file = "jams_matrix_p" + fleet + ".dat"
      mtrx_writer = DataWriter(jams_file)
      mtrx_writer.write_matrix(jams[fleet], fleet)

   # Create a 3D random walk for the initial guesses
   rand_walk = {}
   for fleet in fleets:
      agents = [agent for agent in fleets[fleet]]
      rand_walk[fleet] = {}
      for agent in agents:
         walking_power = 0
         ag_delta = fleets[fleet][agent][6]
         if ag_delta > 0 + 1e-6:
            walking_power = fleets[fleet][agent][4]/(2*N)
         x = [fleets[fleet][agent][1]]*N
         y = [fleets[fleet][agent][2]]*N
         z = [fleets[fleet][agent][3]]*N
         for t in range(1, N):
            theta = random.vonmisesvariate(math.pi, 0)
            x[t] = x[t-1] + walking_power*math.cos(theta)
            y[t] = y[t-1] + walking_power*math.sin(theta)
            z[t] = z[t-1] + walking_power*math.sin(theta)
         rand_walk[fleet][agent] = (x, y, z)

   # Write the initial guesses for the agents using
   # the random walk to compute distances
   print("Computing initial guesses for the agents")
   for fleet in fleets:
      agents = [agent for agent in fleets[fleet]]
      # Write the initial guesses for the agents
      init_guess_file = "init_guess_p" + fleet + ".dat"
      init_guess_writer = InitialGuessWriter(init_guess_file)
      init_guess_writer.write_initial_guess(fleet, agents, N, fleets, rand_walk, comms[fleet], jams[fleet])
      init_guess_writer.close()
   
   # Write the time index set to file
   time_file = "time_set.dat"
   with open("Inst_Data/" + time_file, "w") as fyle:
      fyle.write("###############################################################\n")
      fyle.write("####################### Set time index ########################\n")
      fyle.write("###############################################################\n")
      fyle.write("# Time set\n")
      fyle.write("set Time := ")
      for i in range(N):
         fyle.write(str(i) + " ")
      fyle.write(";\n")

if __name__ == "__main__":
   # Print statement to show the script is running
   print("\n\nRunning setup_data.py")

   # Check if the user has provided the input file
   if len(sys.argv) < 3:
      print("Usage: python setup_data.py <input_file> <# time_steps>")
      sys.exit(1)
   
   # Get the input file
   instance = sys.argv[1]
   discretisation_size = int(sys.argv[2])

   # Call the main function
   main(instance, N=discretisation_size)

   # Print statement to show the script has finished
   print("Finished running setup_data.py\n\n")



