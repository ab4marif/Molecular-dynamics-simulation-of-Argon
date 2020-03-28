import numpy as np
np.seterr(divide='ignore', invalid='ignore')

from Source.functions import *
import Source.variables as vr

def virial(r):
    """ This function will calculate the virial: sum(r_ij * dU/dr_ij), where r_ij is the distance
        between pair of particles. and U is the Lennard Jones Potential U = 4(r^-12 - r^-6) in natural units

    Parameters
    ----------
    r:          (N,D) matrix with current position in x,y,z direction

    L:          length of the box parameter

    Returns:
    --------
    vir:    the virial term sum(r_ij * dU/dr_ij)

    """
    N = r.shape[0]
    normR = normDistance(r)                   # calculate distances
    mask = np.eye(N, dtype = bool)
    normR[mask] = 1                             # make diagonal 1 to avoid division by zero

    deltaU = (24/normR**7 - 48/normR**13)       # calculate dU/dr
    deltaU[mask] = 0                            # make diagonal 0 again
    vir = np.sum(normR*deltaU)
    return vir


def calculateCompressibility(virial):
    """ Uses the viral to calculate the compressibility

    Parameters
    ----------
    virial:          The virial
    N:              Number of particles
    T:              Temperature of the system

    Returns:
    --------
    compressibility

    """

    compressibility = 1 - 0.5*virial/(3*vr.N*vr.T_initial)
    return compressibility

def calculateSpecificHeat(kineticE):
    """ Calculates the Specific heat of the system from the kinetic energy, with the formula as shown in
        lecture notes week 4.

    Parameters
    ----------
    kineticE:           array of size Timesteps with the total kinetic energy of the system at each timestep


    Returns:
    --------
    cv:                 scalar the Specific heat of the system

    """
    average = np.sum(kineticE)/vr.N
    var = np.sum(kineticE**2)/vr.N - average**2
    A = var/(average**2)

    cv = (6*vr.N)/(2-3*vr.N*A)
    #(2/(3*vr.N) -A)**(-1)
    return cv

def specificHeatInteraction(specificHeat, potentialEnergy):
    """ Calculates the Interaction part of the specific heat C^i according formula (3.9)
        of the Lebowitz paper: J. L. Lebowitz, J.K. Percus, and L. Verlet, ‘Ensemble dependence of fluctuations
                                    with application to machine calculations,’ Phys. Rev., 253, 250–254, 1967.

    Parameters
    ----------
    specificHeat:           The total specific heat
    potentialEnergy:        Total potential energy at each time step


    Returns:
    --------
    c_i:                 Interaction part of the specific heat

    """
    average = np.mean(potentialEnergy)
    var = np.mean(potentialEnergy**2) - average**2

    c_i = specificHeat/(vr.N*vr.T_initial**2)*2/3*var

    return c_i


def errorSpecificHeat(N_error, kineticE):
    """ Calculates the error of the specific heat with the bootstrap method. It picks a
        random set from the kinetic energy array, with this subset the specific is calculated.
        This is done N_error times, then the standard deviation of this gives the error of the
        specific heat.

    Parameters
    ----------
    kineticE:           array of size Timesteps with the total kinetic energy of the system at each timestep
    N_error             Number of times to calculate the specific from the subset

    Returns:
    --------
    error_cv:          error of specific heat

    """
    N_subset = 10000
    Cv = np.zeros(N_error)

    for i in range(N_error):
        kineticE_subset = np.random.choice(kineticE,N_subset)
        Cv[i] = calculateSpecificHeat(kineticE_subset)

    average_cv = 1/N_error*np.sum(Cv)
    error_cv = np.sqrt((1/N_error)*np.sum(Cv**2) - average_cv**2)

    return error_cv

def bootstrapError(data, N_error):
    """ Calculates error of a data set with the bootstrap Method. It picks a
        random set from the data, with this subset the mean is calculated.
        This is done N_error times, then the standard deviation of this gives the error of the
        data.

    Parameters
    ----------
    data:           Data set of N points
    N_error         Number of times to calculate the mean data from the subset


    Returns:
    --------
    data_error:      error of the data

    """
    N_subset = 10000
    data_new = np.zeros(N_error)

    for i in range(N_error):
        data_subset = np.random.choice(data,N_subset)       # draw random N_subset elements from data
        data_new[i] = np.mean(data_subset)

    average_data = 1/N_error*np.sum(data_new)
    error_data = np.sqrt((1/N_error)*np.sum(data_new**2) - average_data**2)

    return error_data

def pairCorrelationFunction(hist, bins):
    """ Calculates the pair Correlation function g(r), from the distances between the atoms.
        From the histogram is calculates the mean and uses that for the g(r) from the following
        formula: g(r) = 2V/N(N-1) * <n(r)>/4*pi*r^2* delta_r, where n(r) is the histogram.


    Parameters
    ----------
    hist:           histogram of the distances
    bins:           bins of the histogram


    Returns:
    --------
    pair_correlation_function:      g(r)

    """
    avg_hist = np.mean(hist, axis = 0)              # mean of the histogram
    pair_correlation_function = np.zeros(200)
    delta_r = bins[1]-bins[0]

    A = avg_hist/(4*np.pi*delta_r*bins[0:-1]**2)
    pair_correlation_function = 2*vr.L**3/(vr.N*(vr.N-1))*A

    return pair_correlation_function
