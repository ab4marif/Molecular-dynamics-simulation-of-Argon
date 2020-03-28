import numpy as np
from scipy.optimize import curve_fit
import Source.variables as vr

""" This file contains the functions used in the maincode file to simulate argon gas.
Please refer to the maincode to run the file.
"""


def calculateDistance(r):
    """ Calculates the distance between each particle to any other particle, it takes as input the
    current position and gives as output a matrix with all the distances with the minimal image implementation.

    Parameters
    ----------

    r: (N,D) matrix with current position in x,y,z direction
    L : length of the box parameter

    Returns:
    --------
    deltaR: (N,N,D) matrix with the distances of each particle, where the diagonal of this matrix is 0

    """

    N = r.shape[0]  # Number of particles
    D = r.shape[1]  # Dimension

    bigR = np.broadcast_to(r, (N, N, D))    # broadcast to (N,N,D) matrix to avoid for loops
    bigRT = np.transpose(bigR, (1, 0, 2))

    deltaR = bigR - bigRT                   # subtracting them gives all possible distances

    # minimal image implementation
    deltaR = (deltaR + vr.L/2)%vr.L - vr.L/2

    return deltaR

def normDistance(r):
    """ Calculates the normed between each particle to any other particle by sqrt(x^2 + y^2 + z^2), it takes as input the
    current position and gives as output a matrix with all the normed distances with the minimal image implementation.

    Parameters
    ----------

    r: (N,D) matrix with current position in x,y,z direction
    L : length of the box parameter

    Returns:
    --------
    norm_distance: (N,D) matrix with the normed distances of each particle, where the diagonal of this matrix is 0

    """

    min_distance = calculateDistance(r)
    norm_distance = np.linalg.norm(min_distance, axis = 2)
    return norm_distance

def calculateForce(r):
    """ Calculates the force on each particle due to the Lennard-Jones potential

    Parameters
    ----------

    r: (N,D) matrix with current position in x,y,z direction
    L : length of the box parameter

    Returns:
    --------
    F: (N,D) matrix with the total force on each particle

    """
    N = r.shape[0]  # number of particles

    deltaR = calculateDistance(r) # gives the distances needed for the potential

    # deltaR has zeros on the diagonal so we fill them with ones
    # because we multiply with deltaR again this is valid
    normR = np.linalg.norm(deltaR, axis = 2, keepdims = True)
    mask = np.eye(N, dtype = bool)
    normR[mask,:] = 1


    deltaU = 4*(6/normR**8 - 12/normR**14)*deltaR
    F = np.sum(deltaU, axis = 1 )

    return F


def potentialEnergy(r):
    """This function gives us the total potential energy from the Lennard Jones potential

    Parameters
    ----------

    r: (N,D) matrix with current position in x,y,z direction
    L : length of the box parameter

    Returns:
    --------
    Utot: scalar, The total potential Energy of the System

    """

    N = r.shape[0]                              # number of particles
    deltaR = calculateDistance(r)             # distances of the particles
    normR = np.linalg.norm(deltaR, axis = 2)    # converting the x,y,z coordinates into a r

    mask = np.eye(N, dtype = bool)              # filling the diagonal with ones to avoid division by zero
    normR[mask] = 1

    U = 4*(1/normR**12 - 1/normR**6 )           # calculating potential for each particles
    U[mask] = 0                                 # setting diagonal back to 0

    Utot = 0.5*np.sum(U)                        # calculating total potential of system
    return Utot


def kineticEnergy(velocity):
    return np.sum(0.5*np.linalg.norm(velocity)**2)


def totalEnergy(kineticEnergy, potentialEnergy):
    return kineticEnergy + potentialEnergy


def createLattice(multiplier, unitCell, L):
    """ This function will create a Lattice from a unit cell, by viewing the unitcells in
    lattice as matrix elements.

    Parameters
    ----------

    multiplier: number of times one wants to move a unit cell in D dimensions

    unitCell:   A unit cell in the primitive vector notation

    L :         length of the box parameter

    Returns:
    --------
    lattice:    Gives the positions of the particles in in the lattice within a box of size L

    """

    latticeStructure = np.ones((multiplier,multiplier,multiplier))          # create the structe of the lattice matrix
    allUnitCells = np.argwhere(latticeStructure)
    allUnitCells = np.broadcast_to(allUnitCells, (1,multiplier**3,3))

    movingUnitCells = np.swapaxes(allUnitCells, 0, 1) + unitCell            # add the unitcell to translate it in the directions
    positions = movingUnitCells.transpose(2,0,1).reshape(3,-1).transpose()  # gives the positions of the particles in the lattice
    lattice = L/multiplier*positions

    return lattice
