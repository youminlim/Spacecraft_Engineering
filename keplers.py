#%% Solver for Kepler's equation

#%% Import libraries
import numpy as np
from numpy import sin
from numpy import cos

#%% Functions
def solver(M, e):
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
    
    #%%% Define functions
    def plotter(x, y):
        
        # Define earth to plot
        r = r_earth
        theta = np.linspace(0, 2*np.pi, 1000)
        x_e = r * cos(theta)
        y_e = r * sin(theta)
        
        fig, ax = plt.subplots()
        earth_plot, = plt.plot(x_e, y_e, label = 'Earth')
        trajectory, = plt.plot(x, y, label = 'Trajectory')
        ax.legend(handles = [earth_plot, trajectory], loc = 'upper right')
        plt.axis('square')
        plt.show()

    #%%% Define program constants
    ''' Mission profile : initial altitude 800km, circular orbit'''
    altitude = 800
    r_earth = 6379
    G = 6.67e-11
    M = 5.972e24
    [e, a, i, RAAN, omega, M0] = [0, altitude+r_earth, 0, 0, 0, 60]
    M0 = np.deg2rad(M0)

    #%%% Program
    mu = G * M * 10 ** -9
    period = 2 * np.pi * np.sqrt(a**3 / mu)
    [tmin, tmax] = [0, period]
    dt = np.linspace(tmin, tmax, 1000)
    
    # COEs
    E = np.zeros(len(dt))
    x = np.zeros(len(dt))
    y = np.zeros(len(dt))
    n = np.sqrt(mu / a**3)
    
    # Calculate eccentric anomaly
    for i in range(len(dt)):
        M = M0 + (n * dt[i])
        E[i] = solver(M, e)
        x[i] = a * (cos(E[i]) - e)
        y[i] = a * sin(E[i]) * np.sqrt(1 - e**2)
    
    plotter(x, y)
    
    
    