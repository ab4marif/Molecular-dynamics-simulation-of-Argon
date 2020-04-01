import numpy as np
from matplotlib import pyplot as plt
import time

import Source.functions as fnc
import Source.dataProcessing as dp
from Source.initialization import *
from Source.timeEvolution import *
from plots import *

""" This is the main file of the simulation. Running this file will initialize the system and let it evolve over time,
    with parameters rho and T find at the end of this file. it will output the compressibility, specific heat and if set to
    True also the plots of the energies and pari correlation function.
"""

start_time = time.time()

def run_simulation(density, temperature, plots ):

    vr.rho = density
    vr.T_initial = temperature
    vr.L = (vr.N/vr.rho)**(1/vr.D)

    # initilazing the simulation
    (v, r, F, kineticE, potentialE,
        Temperature, vir, hist, bins) = initilaze()

    # letting the simulation evolve in time and obtaining the variables
    (kineticE, potentialE,
        vir, hist, bins) = timeEvolution(v, r, F, kineticE, potentialE,
                                                Temperature, vir, hist, bins)

    totalEnergy = fnc.totalEnergy(kineticE, potentialE)

    # Calculating Observables
    gr = dp.pairCorrelationFunction(hist, bins)

    # obtaining the compressibility of the system and the error
    compressibility = dp.calculateCompressibility(vir)
    average_p = np.mean(compressibility)
    error_p = dp.bootstrapError(compressibility, 100)
    print(f"compressibility = {average_p} with error of {error_p} ")

  
    # Calculating Specific Heat and error
    cv = dp.calculateSpecificHeat(kineticE)
    N_error = 100
    error_cv = dp.errorSpecificHeat(N_error, kineticE)
    cv_SI = cv*vr.kb/vr.mass*10**-3                            # converting cv into SI units J/kg*K
    print(f"Specific Heat = {cv} with error of {error_cv}")
    print(fr"specific heat SI units = {cv_SI} [J/kgK]")
    
    #Plots
    if plots == True:
        energyPlot(kineticE, potentialE, totalEnergy)
        pairCorrelationPlot(gr, bins, density, temperature)



# simulation parameters
density = 1
temperature = 1
plots = True                   # set to False if don't want to show plots

run_simulation(density, temperature, plots)

elapsed_time = time.time() - start_time
print(f"time taken for simulation :{elapsed_time}")
