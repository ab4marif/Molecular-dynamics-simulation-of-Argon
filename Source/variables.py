import numpy as np


# natural constanst
eps = 1 #119.8*1.38*10**-23
sigma = 1 #3.405E-10
mass = 6.6E-26
kb = 1.38*10**-23

# primitive unit cell of an FCC lattice
unitCell = np.array([[[0, 0, 0], [0.5, 0.5, 0],
                           [0.5, 0, 0.5], [0, 0.5, 0.5] ]])

multiplier = 5                           # amount of multiplication of the unitcell in all directions

D = unitCell.shape[2]                   # dimensions in accordence with the unit cell
N = multiplier**D*unitCell.shape[1]     # number of particles consistent witht the unit cell
rho = 0.88                              # density in units of sigma^3
L = (N/rho)**(1/D)                      # Length of the box


dt = 3e-3                               # timestep for evolving the simulation
simulation_time = 20
time = np.arange(0, simulation_time, dt)
Timesteps = time.shape[0]               # number of timesteps

# initial temperature
T_initial = 1.095

print(N)
