import numpy as np
import  Source.variables as vr
import Source.functions as fnc

def initilaze():
    """ This is to initilaze the simulation. It creates matrices for variables and calculates the
        first elements of the matices.


    Returns:
    --------
    v:                  velocity obtained from gaussian distribution

    r:                  position of each atom structed in FFC lattice

    F:                  force acting on each atom due to Lennard jones potential

    kineticE:           kinetic energy of the System

    potentialE:         potential energy of the System

    Temperature:        temperature according to the equipartion theory

    vir:                the virial

    hist:               histogram of the distances between atoms

    bins:               bins of the histogram

    """
    # initilazing the coordinates and velocities
    meanVelocity = 0
    devVelocity = np.sqrt(vr.T_initial)
    v = np.random.normal(meanVelocity, devVelocity, (vr.N,vr.D))          # create intitial velocity from gauss distribution
    r = np.zeros((vr.N, vr.D, vr.Timesteps))                                 # create matrix to store positions
    r[:,:,0] = fnc.createLattice(vr.multiplier, vr.unitCell, vr.L)           # initilaze position in FFC lattice
    F = fnc.calculateForce(r[:,:,0])

    # initilazing the Energies
    kineticE = np.zeros(vr.Timesteps)
    potentialE = np.zeros(vr.Timesteps)
    kineticE[0] = fnc.kineticEnergy(v)
    potentialE[0] = fnc.potentialEnergy(r[:,:,0])


    # Restabilization of the energy
    Temperature = np.zeros(vr.Timesteps)
    Temperature[0] = kineticE[0]/((vr.N-1)*3/2)    # initilaze the Temperature with equipartion


    # initilazing the viral  to calcute the compressibility
    vir = np.zeros(vr.Timesteps)
    hist = np.zeros((vr.Timesteps, 200))
    normDis = fnc.normDistance(r[:,:,0])
    hist[0], bins = np.histogram(normDis, 200)

    return (v, r, F, kineticE, potentialE,
        Temperature, vir, hist, bins)
