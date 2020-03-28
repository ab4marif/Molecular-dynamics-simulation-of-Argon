import numpy as np
from matplotlib import pyplot as plt

import Source.functions as fnc
import Source.dataProcessing as dp
import Source.variables as vr
from Source.timeEvolution import *


def energyPlot(kineticEnergy, potentialEnergy, totalEnergy):
    plt.plot(vr.time, kineticEnergy, label = "Kinetic" )
    plt.plot(vr.time, potentialEnergy, label = "Potential")
    plt.plot(vr.time, totalEnergy, label =  "Total")
    plt.ylabel(r"Energy [$\epsilon$]")
    plt.xlabel("time [s]")
    plt.title("Energy of the system")
    plt.legend()
    plt.savefig("plots/Energy_over_time.png")
    plt.show()


def pairCorrelationPlot(gr, bins, density, temperature):
    plt.plot(bins[1:201], gr, label = fr"T = {temperature}, $\rho$ = {density}")
    plt.xlabel("r")
    plt.ylabel("g(r)")
    plt.title("Pair Correlation Function")
    plt.legend()
    plt.savefig("plots/pair_correlation.png")
    plt.show()
