'''
This script is used to plot the solution of the jamming game.

Instance data is also found in the "Inst_Data" folder. The script
needs to read the files "data_team<i>.dat", "fleet_team<i>.dat" and 
"time_set.dat" to plot the solution.

- data_team<f>.dat: contains the initial data for each agent in fleet
f in the following format:
###############################################################
####################### Set parameters ########################
###############################################################
param: x0 y0 z0 maxpow c delta :=
A 1.0 1.0 1.0 1.0 0.2 0.5
B -1.0 1.0 1.0 1.0 0.2 0.49
;

- fleet_team<f>.dat: contains the agents in fleet f in the following
format (example with two agents in fleet 1):
# Define fleet 1
set FLEET1 := 'A' 'B' ;

- time_set.dat: contains the time index set (example with N=5)
###############################################################
####################### Set time index ########################
###############################################################
# Time set
set Time := 0 1 2 3 4 ;

Solution files have been saved in the "Solution" folder. There
are many types of solution files. You'll find a brief description
of each type below:

- obj_convergence_p<i>.dat: contains the objective function value 
at each step of the block-descent algorithm

- power_p<f>.dat: contains the power spent by each agent in fleet f
to communicate with each other agent in fleet f and to jam each agent 
not in the same fleet. The file has the following format:
   for i in fleet:
      for j in fleet:
         for t in time:
            power[i,j,t]
      for k not in fleet:
         for l not in fleet:
            if l != k:
               for t in time:
                  power_J[i,k,l,t]

- reg_value_p<f>.dat: contains the regularization value at each step
for fleet f. The format is as follows:
   for i in fleet:
      for t in time:
         x_hat[i,t]
      for t in time:
         y_hat[i,t]
      for t in time:
         z_hat[i,t]

- solution_p<f>.dat: contains the solution of the jamming game at each
time step for the fleet f. The format is as follows:
   for i in fleet:
      for t in time:
         x[i,t]
      for t in time:
         y[i,t]
      for t in time:
         z[i,t]
      for j in fleet:
         if j != i:
            for t in time:
               power[i,j,t]
      for k not in fleet:
         for l not in fleet:
            if l != k:
               for t in time:
                  power_J[i,k,l,t]
      for j in fleet:
         for t in time:
            distance[i,j,t]

The scalar values can be read from the file "scalars.dat". The format
is the following:

###############################################################
####################### Set parameters ########################
###############################################################

data;
param varrho := 1;
param alpha := 1;
param sigmasqr := 1;
param eps := 1e-6;
param reg := 1e-3;
'''

import matplotlib.pyplot as plt
import matplotlib.animation as anime
import numpy as np

# Declare global variables
rho, alpha, sigmasqr, eps, reg = 0, 0, 0, 0, 0

# Read constants
def read_constants():
   file = open("Inst_Data/scalars.dat", "r")
   # Skip header
   for i in range(5):
      next(file)

   #Constants are global variables
   global rho, alpha, sigmasqr, eps, reg
   
   # Read the constants
   rho = float(file.readline().split(":=")[1].split(';')[0])
   alpha = float(file.readline().split(":=")[1].split(';')[0])
   sigmasqr = float(file.readline().split(":=")[1].split(';')[0])
   eps = float(file.readline().split(":=")[1].split(';')[0])
   reg = float(file.readline().split(":=")[1].split(';')[0])

# Function to read the fleet data
def read_fleets():
   fleets = []
   for i in range(1, 3):
      file = open("Inst_Data/fleet_team" + str(i) + ".dat", "r")
      # Skip first line
      next(file)
      # Read fleet and extract agents (each agent is in single quotes)
      fleet_i = file.read().split("\n")[0].split(" := ")[1].split(" ;")[0].split(" ")
      # Remove the single quotes and save the data
      fleets.append([agent[1] for agent in fleet_i])

   # Return the list of fleets
   return fleets

# Function to read the time set
def read_time_set():
   file = open("Inst_Data/time_set.dat", "r")
   # Skip header
   for i in range(4):
      next(file)
   timeSet = file.read().split("\n")[0].split(" := ")[1].split(" ;")[0].split(" ")
   # Convert to integers
   timeSet = [int(t) for t in timeSet]

   # Return the time set
   return timeSet

# Function to read the parameters
def read_fleet_params(fleet):
   file = open("Inst_Data/data_team" + str(fleet) + ".dat", "r")
   # Skip header
   for i in range(4):
      next(file)
   # Read the data
   fleet_params = {}
   for line in file:
      if line == ";\n":
         break
      data = line.split()
      agent = data[0]
      x0 = float(data[1])
      y0 = float(data[2])
      z0 = float(data[3])
      maxpow = float(data[4])
      c = float(data[5])
      delta = float(data[6])
      fleet_params[agent] = (x0, y0, z0, maxpow, c, delta)

   # Return the fleet parameters
   return fleet_params

def read_solution(fleet_idx, fleet, timeSet, all_fleets):
   # Open the solution file
   fyle = open("Solution/solution_p" + str(fleet_idx) + ".dat", "r")

   # Read the solution from file ordered as in the documentation
   # That is: fleet_sol = {} is a dictionary of size len(fleet) with
   # keys being the agents in the fleet and values being numpy arrays.
   # For each agent, the dimensions of each array is as follows:
   #   x <- (len(timeSet), 1)
   #   y <- (len(timeSet), 1)
   #   z <- (len(timeSet), 1)
   #   power <- (len(fleet) - 1, len(timeSet))
   #   power_J <- (len(fleet2), len(fleet2), len(timeSet))
   #   dist <- (len(fleet) - 1, len(timeSet))
   fleet_sol = {}
   for agent in fleet:
      x = [] # x coordinates
      for t in timeSet:
         x.append(float(fyle.readline()))
      y = [] # y coordinates
      for t in timeSet:
         y.append(float(fyle.readline()))
      z = [] # z coordinates
      for t in timeSet:
         z.append(float(fyle.readline()))
      power = [] # communication power
      for ag2 in fleet:
         if ag2 != agent:
            aux = []
            for t in timeSet:
               aux.append(float(fyle.readline()))
            power.append(aux)
      power_J = [] # jamming power
      for f2 in all_fleets:
         if f2 != fleet:
            for enmy1 in f2:
               for enmy2 in f2:
                  if enmy1 != enmy2:
                     aux = []
                     for t in timeSet:
                        aux.append(float(fyle.readline()))
                     power_J.append(aux)
      dist = [] # distance between agents
      for ag2 in fleet:
         if ag2 != agent:
            aux = []
            for t in timeSet:
               aux.append(float(fyle.readline()))
            dist.append(aux)
      # Store the solution for the agent
      fleet_sol[agent] = [x, y, z, power, power_J, dist]

   # Return the solution
   return fleet_sol

def checkSolutions(fleet_params, fleet_sol):
   for fleet in fleet_sol:
      fleet_index = fleet_sol.index(fleet)
      for agent in fleet:
         x0, y0, z0, maxpow, c, delta = fleet_params[fleet_index][agent]
         x, y, z, pow, pow_J, dist = fleet[agent]
         # Compute energy spent by moving
         moving_pow = 0
         for t in range(1, len(x)):
            dist = np.sqrt((x[t] - x[t-1])**2 + (y[t] - y[t-1])**2 + (z[t] - z[t-1])**2)
            moving_pow += c * dist
         # Compute energy spent by communicating and jamming
         comms_pow = sum([sum(p) for p in pow])
         jam_pow = sum([sum(p) for p in pow_J])

         # Compute total power
         total_pow = moving_pow + comms_pow + jam_pow
            
         print("Agent ", agent, " moving power =", moving_pow)
         print("Agent ", agent, " comms power =", comms_pow)
         print("Agent ", agent, " jamming power =", jam_pow)
         RED = '\033[91m'
         END = '\033[0m'
         print(RED + "Agent ", agent, " total power =", total_pow, END)
         print("\n")

         assert total_pow < maxpow, "Agent " + agent + " exceed max pow: " + str(total_pow) + " > " + str(maxpow)

# Function to animate the trajectories of the agents by
# plotting the x and y coordinates of each agent at each
# time step using arrows
def animate_trajectories(frame, ax, fleet_sol):
   prev_fr = 0
   if frame > 1:
      prev_fr = frame - 2
   for fleet in fleet_sol:
      for agent in fleet:
         x, y, z, pow, pow_J, dist = fleet[agent]
         xp = x[prev_fr:frame]
         yp = y[prev_fr:frame]
         zp = z[prev_fr:frame]
         # Plot path
         ax.plot(xp, yp, zp, c = 'k', alpha=0.7)

# Main function that triggers the plotting
def main():
   # Read constants
   read_constants()

   # Read instance files
   fleets = read_fleets()
   timeSet = read_time_set()
   fleet_params = [None for i in range(len(fleets))]
   for f in range(len(fleets)):
      fleet_params[f] = read_fleet_params(f+1)

   # Read solution files
   fleet_sol = [None for i in range(len(fleets))]
   for f in range(len(fleets)):
      fleet_sol[f] = read_solution(f+1, fleets[f], timeSet, fleets)

   # Check solutions
   checkSolutions(fleet_params, fleet_sol)

   # Plot each agent's initial positions in 3D
   figure = plt.figure()
   ax = figure.add_subplot(111, projection='3d')

   # Plot agents in the same team with similar colours
   for f in range(len(fleets)):
      # Set the colours for the agents in the fleet
      shade = 0
      agent_colour = (0, 0, 0)
      for agent in fleets[f]:
         # Color shade
         shade = np.random.choice(range(256))
         # Update agent_colour
         if f == 0: # 0 for blue, 1 for red
            agent_colour = (0, shade/256, 0)
         else:
            agent_colour = (shade/256, 0, 0)
         x0, y0, z0, _, _, _ = fleet_params[f][agent]
         ax.plot(x0, y0, z0, "o", c = agent_colour, label=agent)
   
   # Common settings
   global box_s 
   box_s = 2
   ax.set_xlim(-box_s, box_s)
   ax.set_ylim(-box_s, box_s)
   ax.set_zlim(-box_s, box_s)
   ax.set_xlabel("x")
   ax.set_ylabel("y")
   ax.set_zlabel("z")
   ax.set_aspect('equal', 'box')

   # Animation constants
   speed = 800 # miliseconds
   num_frames = len(timeSet) + 1

   # Plot the agents' paths using arrows
   animation1 = anime.FuncAnimation(fig=figure, func=animate_trajectories, fargs=(ax, fleet_sol,), frames=num_frames, interval=speed, repeat=False)
   # Save animation
   animation1.save("Solution/trajectories.gif")

   # Show plot
   plt.show()


if __name__ == "__main__":
   main()
