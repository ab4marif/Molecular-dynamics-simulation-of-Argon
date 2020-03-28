import numpy as np
import Source.functions as fnc
import Source.variables as vr
import Source.dataProcessing as dp


def timeEvolution(v, r, F, kineticE, potentialE,
    Temperature, vir, hist, bins):
    """ This lets the simulation evolve over time, at ech timestep calculating the variables.
    """

    checkTime = np.arange(0,vr.Timesteps, 10)      # create times to check for the Temperature
    reachedTemperature = False                  # set a boolean for the disired Temperature

    for i in range(1,vr.Timesteps):

        # calcute the energies and Temperature at each timestep
        kineticE[i] = fnc.kineticEnergy(v)
        potentialE[i] = fnc.potentialEnergy(r[:,:,i-1])
        Temperature[i] = kineticE[i]/((vr.N-1)*3/2)

        # the virial at each timestep
        vir[i] = dp.virial(r[:,:,i-1])

        # verlet algorithm
        r[:,:,i] = r[:,:,i-1] + vr.dt*v + 0.5*fnc.calculateForce(r[:,:,i-1])*vr.dt**2
        r[:,:,i] = r[:,:,i]%vr.L

        normDis = fnc.normDistance(r[:,:,i])
        hist[i], _ = np.histogram(np.sort(normDis.flatten()), bins)

        # at some times check of the desired Temperature has been reached, if not multiply the
        # velocities with lambda
        if i in checkTime and not reachedTemperature:
            if np.abs(((vr.T_initial- Temperature[i]))/vr.T_initial) < 0.01:
                reachedTemperature = True
                v = v + 0.5*vr.dt*(fnc.calculateForce(r[:,:,i])+ fnc.calculateForce(r[:,:,i-1]))
            else:
                v = np.sqrt(((vr.N-1)*3/2*vr.T_initial)/kineticE[i])*v
        else:
            v = v + 0.5*vr.dt*(fnc.calculateForce(r[:,:,i])+ fnc.calculateForce(r[:,:,i-1]))

    return (kineticE, potentialE, vir, hist, bins)
