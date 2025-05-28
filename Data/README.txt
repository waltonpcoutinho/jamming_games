Fleet #<fleet number>
Agents <number of agents>
//for each agent in Fleet
agent_indx x0 y0 z0 maxpow c delta
A 1 1 1 2 0.5 0.1
B 1 -1 1 2 0.5 0.1
Fleet #<fleet number>
Agents <number of agents>
//for each agent in fleet
agent_indx x0 y0 z0 maxpow c delta
a -1 1 1 2 0.5 0.1
b -1 -1 1 2 0.5 0.1
//communications adjacency matrix #<fleet number>
//not necessarily symmetric
Comms_1 [0 1; 1 0]
  A B
A 0 1
B 1 0
//jamming adjacency matrix #<fleet number>
//not necessarily symmetric
Jams_1 [1 1; 1 1]
//Represents:
  a b
A 1 1
B 1 1
// !!!!!!!!!!!!!!!!
// We assume a compact form of the jamming matrix were
// entry i,j means that i jams all arcs comming into j
// ex: "A" jamming "a" means that A jams (b,a), (c,a), ...
// !!!!!!!!!!!!!!!!
//communications adjacency matrix #<fleet number>
//not necessarily symmetric
Comms_2 [0 1; 1 0]
//Represents:
  a b
a 0 1
b 1 0
//jamming adjacency matrix #<fleet number>
//not necessarily symmetric
Jams_1 [1 1; 1 1]
//Respresents
  A B
a 1 1
b 1 1
// !!!!!!!!!!!!!!!!
// We assume a compact form of the jamming matrix were
// entry i,j means that i jams all arcs comming into j
// ex: "A" jamming "a" means that A jams (b,a), (c,a), ...
// !!!!!!!!!!!!!!!!

// !!!!!!!!!!!!!!!!
Jams_i(i,j) = 1 <=> agent i in fleet F is able to jam
all incoming communications of agent j in G
// !!!!!!!!!!!!!!!!
