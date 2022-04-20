#%% Program for satellite simulation

#%%% Importing libraries
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp as ode

#%%% Functions
def twobody_orbit(t, y, mu):
    
    # Unpack state vector
    rx, ry, rz, vx, vy, vz = y
    r = np.array([rx, ry, rz])
    
    # Norm of position vector
    r_norm = np.linalg.norm(r)
    
    # Calculate acceleration
    ax, ay, az = -r*mu / r_norm**3
    
    # Return complete differentiated state vector
    return [vx, vy, vz, ax, ay, az]

def planet_radius(planet_name):
    
    if planet_name == "Earth":
        R = 6371
    else:
        R = 0
        
    return R

def planet_mu(planet_name):
    
    # Initialise gravitational constant
    G = 6.67 * 10 ** -11
    
    if planet_name == "Earth":
        M = 5.972 * 10 ** 24
        mu = G * M
    else: 
        mu = 0
        
    return mu

def plotter(state_matrix):
    
    # Unpack state matrix 
    rx = state_matrix[:,0]
    ry = state_matrix[:,1]
    rz = state_matrix[:,2]
    
    ax = plt.figure().add_subplot(projection='3d')
    ax.plot(rx, ry, rz)

#%%% Main program
if __name__ == "__main__":
    
    # Obtain planet properties
    planet = "Earth"
    R = planet_radius(planet)
    mu = planet_mu(planet)
    
    # Initial orbit parameters
    altitude = R + 500 # km
    v_initial = np.sqrt(mu/altitude) # km / s
    
    # Initial state vector
    r0 = [altitude, 0, 0]
    v0 = [0, v_initial, 0]
    y0 = r0 + v0
    
    # Time step
    # (2 * np.pi * np.linalg.norm(r0)) / np.sqrt(mu)
    tspan = np.array([0, 100 * 60]) # Period of orbit
    dt = 100 # Timestep
    n = int(np.ceil(tspan[1] / dt))
    
    # Initialise state arrays to store values
    y = np.zeros((n,6))
    t = np.zeros((n,1))
    y[0] = np.array(y0)
    
    # Initialise ODE solver
    solver = ode(twobody_orbit, tspan, y0, args=[mu])
    t = solver.t
    y = solver.y
    y = np.transpose(y)
    
    plotter(y)