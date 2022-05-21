#%% Solver for Kepler's equation

#%% Import libraries
import numpy as np
from numpy import sin
from numpy import cos

#%% Functions
def kepler(M, e):
    # A solver to compute the eccentric anomaly at any given time via the newton raphson iteration
    E = M
    error = 1
    
    while error >= 0.0001:
        f = E - (e * sin(E)) - M
        df = 1 - (e * cos(E))
        
        # newton raphson formula
        E_updated = E - (f / df)
        
        # update error and eccentric anomaly
        error = E_updated - E
        E = E_updated
        
    return E

#%% Test program
if __name__ == "__main__" :
    
    #%%% Import libraries
    import matplotlib.pyplot as plt

    #%%% Define program constants
    ''' Mission profile : initial altitude 800km, circular orbit'''
    altitude = 800
    r_earth = 6379
    G = 6.67e-11
    M = 5.972e24
    e = 0.9

    #%%% Program
    mu = G * M * 10 ** -9
    a = altitude + r_earth
    period = 2 * np.pi * np.sqrt(a**3 / mu)
    [tmin, tmax] = [0, period]
    timesteps = np.linspace(tmin, tmax, 1000)
    
    # COEs
    n = 2 * np.pi / period
    M = n * timesteps    
    E = np.zeros(len(M))
    x = np.zeros(len(M))
    y = np.zeros(len(M))
    
    # Calculate eccentric anomaly
    for i in range(len(M)):
        E[i] = kepler(M[i], e)
        x[i] = a * (cos(E[i]) - e)
        y[i] = a * sin(E[i]) * np.sqrt(1 - e**2)
    
    plt.plot(x,y)
    plt.axis('square')
    plt.show()
    
    
    
    