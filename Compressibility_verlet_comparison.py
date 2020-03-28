import numpy as np
from matplotlib import pyplot as plt

from Source.initialization import *
from Source.timeEvolution import *

""" This file will calculate the compressibility for different densities and temperatures.
    These values are then plotted witht the values from the Verlet Paper:
     Verlet, L. (1967, Jul). Computer "experiments" on classical fluids.
      i. thermodynamical properties of lennard-jones molecules. Phys. Rev., 159 , 98–103.
"""


def verletComparisonPressure(density, temperature):
    """ Calculates the average compressibility and the error, given a density and temperature.
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

    compressibility = dp.calculateCompressibility(vir)
    average_p = np.mean(compressibility)
    error_p = dp.bootstrapError(compressibility, 100)

    return (average_p, error_p)


# densities, temperatures and compressibility as in the verlet paper
density = np.array([0.88, 0.88, 0.85, 0.75, 0.65, 0.55, 0.45, 0.4, 0.35])
temperature = np.array([1.095, 0.940, 2.889, 2.849, 2.557, 2.645, 4.625, 1.462, 1.620])
verletcompressibilty = np.array([3.48, 2.72, 4.20, 3.10, 2.14, 1.63, 1.68, 0.41, 0.58])


Compressibility = np.zeros(temperature.shape[0])
error_p = np.zeros(temperature.shape[0])

for i in range(temperature.shape[0]):
    Compressibility[i], error_p[i] = verletComparisonPressure(density[i], temperature[i])



def compressibilityPlot(Compressibility, verletcompressibilty):
    number = np.arange(Compressibility.shape[0], dtype = "int")

    plt.plot(number, Compressibility,linestyle='--', marker='o', label = "Simulation")
    plt.plot(number, verletcompressibilty,linestyle='--', marker='o', label = "Verlet paper")
    plt.xlabel("simulation")
    plt.ylabel(r"Compressibility $\beta P/\rho$")
    plt.title("Compressibility Comparison Verlet Paper")
    plt.legend()
    #plt.savefig("plots/Compressibility_comparison.png")
    plt.show()



compressibilityPlot(Compressibility, verletcompressibilty)
print(Compressibility)
print(error_p)
