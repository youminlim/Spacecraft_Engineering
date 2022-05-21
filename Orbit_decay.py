#%% Import external libraries
import numpy as np
from numpy import sin
from numpy import cos
import matplotlib.pyplot as plt 
import keplers

#%% Define program constants
# chosen astral body : Earth
G = 6.67e-11          # universal gravitational constant - m³/kg*s²
M = 5.972e24          # earth mass - kg
r_earth = 6379        # earth radius - km
mu = G * M * 10 ** -9 # earth gravitaional parameter - km³/s²

# initialise spacecraft parameters
BC = 0.00000229         # ballistic coefficient
m = 1000                # spacecraft mass - kg
altitude = 100          # spacecraft initial altitude - km

#%% Functions
def atmos(r):
    # A function to find the density of air at various altitudes
    
    # re-scaling
    altitude = (r - r_earth) * 10 ** 3
    
    # define constants for atmosphere model
    p0 = 101325 # sea level atmospheric pressure Pa
    T0 = 288.15 # sea level atmospheric temperature K
    g = 9.80665 # earth surface gravitational acceleration m/s²
    R = 8.31446 # ideal gas constant J/mol*K
    M = 0.02897 # molar mass of dry air kg/mol
    
    # altitude upper limit levels
    troposphere = 11 * 10 ** 3
    tropopause = 20 * 10 ** 3
    stratosphere1 = 32 * 10 ** 3
    stratosphere2 = 47 * 10 ** 3
    stratopause = 51 * 10 ** 3
    
    if altitude <= troposphere:
        # constants within the troposphere
        L = 0.0065  # temperature lapse rate K/m
    
        # governing equations
        T = T0 - (L * altitude)
        p = p0 * (1 - ((L*altitude)/T0)) ** (g*M/R*L)
        rho = (p * M)  / (R * T)
        
    elif troposphere < altitude <= tropopause:
        rho = (((altitude - troposphere)/(tropopause-troposphere))*(0.088 - 0.3639)) + 0.3639

    elif tropopause < altitude <= stratosphere1:
        rho = (((altitude - tropopause)/(stratosphere1-tropopause))*(0.0132 - 0.088)) + 0.088
        
    elif stratosphere1 < altitude <= stratosphere2:
        rho = (((altitude - stratosphere1)/(stratosphere2-stratosphere1))*(0.002 - 0.0132)) + 0.0132
        
    elif stratosphere2 < altitude <= stratopause:
        rho = (((altitude - stratosphere2)/(stratopause-stratosphere2))*(0 - 0.002)) + 0.002
        
    elif stratopause < altitude :
        rho = 0
        
    return rho * 10 ** 9

def period(r):
    # computes the keplerian period of an orbit at a specific orbit radius
    return 2 * np.pi ** np.sqrt(r**3 / mu)

def velocity(r, a):
    # computes the orbital velocity at a specific orbit radius
    return np.sqrt(mu * ((2/r) - (1/a)))

def accel_drag(r, a):
    # computes the deceleration due to atmospheric drag at a specific orbit radius
    rho = atmos(r)
    v = velocity(r, a)
    return -BC * 0.5 * rho * v ** 2

def drdt(r, a):
    # calculates the change in orbital radius due to atmospheric drag
    return accel_drag(r, a) * period(r) / np.pi

    
#%% Main program

if __name__ == "__main__":
    
    # initialise solver parameters
    r = altitude + r_earth
    r_propagate = [r]
    
    while altitude >= 0 :
        
        # orbit decay propagation
        
        # update orbit radius and altitude
        r = r - dr
        altitude = r - r_earth
        r_propagate.append(r)
    
    