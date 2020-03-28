import numpy as np
from matplotlib import pyplot as plt

from Source.initialization import *
from Source.timeEvolution import *

""" This file will calculate the Interaction part of the specific heat for different densities and temperatures.
    The values will then be plotted against the values from Lebowitz paper :
     L. Lebowitz, J.K. Percus, and L. Verlet, ‘Ensemble dependence of fluctuations
    with application to machine calculations,’ Phys. Rev., 253, 250–254, 1967.
"""


def specficHeatComparison(density, temperature):
    """ Calculates the Interaction specific heat and the error, given a density and temperature.
        Using the initilaze and timeEvolution functions.
    """
    vr.rho = density
    vr.T_initial = temperature
    vr.L = (vr.N/vr.rho)**(1/vr.D)

    (v, r, F, kineticE, potentialE,
        Temperature, vir, hist, bins) = initilaze()

    (kineticE, potentialE,
        vir, hist, bins) = timeEvolution(v, r, F, kineticE, potentialE,
                                                Temperature, vir, hist, bins)

    cv = dp.calculateSpecificHeat(kineticE)
    error_cv = dp.errorSpecificHeat(100, kineticE)
    cv_interaction = dp.specificHeatInteraction(cv, potentialE)
    return cv_interaction, error_cv

def specificHeatPlot(cv_interaction, lebowitsCv):
    number = np.arange(cv_interaction.shape[0], dtype = "int")

    plt.plot(number, cv_interaction, linestyle='--', marker='o', label = "Simulation")
    plt.plot(number, lebowitsCv, linestyle='--', marker='o', label = "Lebowitz paper", color = "r")
    plt.xlabel("simulation")
    plt.ylabel(r" $C^i$")
    plt.title(r"Interaction part Specific Heat $C^i$ Comparison Lebowitz Paper")
    plt.legend()
    plt.savefig("plots/Specific_Heat_comparison.png")
    plt.show()


# density, temp and Interaction specific heat according to the Lebowitz paper
density = np.array([0.85, 0.85, 0.85, 0.75, 0.75, 0.45, 0.45, 0.45])
temperature = np.array([2.89, 2.20, 1.21, 2.84, 0.827, 4.62, 2.93, 1.71])
lebowitsCv = np.array([0.59, 0.78, 0.84, 0.47, 0.89, 0.25, 0.22, 0.46])

# calculate cv_interaction for all densities and temperatures
cv_interaction = np.zeros(density.shape[0])
error_cv = np.zeros(density.shape[0])
for i in range(temperature.shape[0]):
    cv_interaction[i], error_cv[i] = specficHeatComparison(density[i], temperature[i])




specificHeatPlot(cv_interaction, lebowitsCv)
print(cv_interaction)
print(error_cv)
