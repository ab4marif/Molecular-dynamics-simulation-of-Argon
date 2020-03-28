import numpy as np
from matplotlib import pyplot as plt

import Source.functions as fnc
import Source.dataProcessing as dp
from Source.initialization import *
from Source.timeEvolution import *
import Source.variables as vr


def calculateCorrelation(density, temp):

    vr.rho = density
    vr.T_initial = temp
    vr.L = (vr.N/vr.rho)**(1/vr.D)

    (v, r, F, kineticE, potentialE,
        Temperature, vir, hist, bins) = initilaze()

    (kineticE, potentialE,
        vir, hist, bins) = timeEvolution(v, r, F, kineticE, potentialE,
                                                Temperature, vir, hist, bins)
    gr = dp.pairCorrelationFunction(hist, bins)

    return (gr, bins)

def corr_plot(gr, bins):

    plt.figure()
    plt.plot(bins[1:201], gr[0], label=fr"$\rho = {density[0]}$, T = {temp[0]}")
    plt.plot(bins[1:201], gr[1], label=fr"$\rho = {density[1]}$, T = {temp[1]}")
    plt.plot(bins[1:201], gr[2], label=fr"$\rho = {density[2]}$, T = {temp[2]}")
    plt.xlabel('r')
    plt.ylabel('g(r)')
    plt.title('Pair Correlation Function of Different Phases')
    plt.legend()
    plt.savefig("plots/pair_correlation_plot.png")
    plt.show()

temp = np.array([1.095, 0.8, 3])
density = np.array([0.88, 1.5, 0.5])

gr = np.zeros((3,200))

for i in range(len(temp)):
    gr[i], bins = calculateCorrelation(density[i], temp[i])



corr_plot(gr, bins)
