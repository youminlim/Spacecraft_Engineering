#%% Variation of thrust against time to deorbit

#%%% Load libraries
import numpy as np
from numpy import sqrt
import matplotlib.pyplot as plt

#%%% Define global constants

# Earth parameteres
M = 5.972 * 10 ** 24
G = 6.67 * 10 ** -11
g = 9.81
r_earth = 6379
mu = G * M * 10 ** -9

# Satellite parameters
m = 2000
Isp = 220

#%%% Functions
def hohmann(r_i, r_f):
    # Calculates delta v required to perform a hohmann transfer
    a = (r_i + r_f) / 2
    v_i = sqrt(mu / r_i)
    v_1 = sqrt(mu * ((2/r_i) - (1/a)))
    v_2 = sqrt(mu * ((2/r_f) - (1/a)))
    v_f = sqrt(mu / r_f)

    return (v_i - v_1) + (v_2 - v_f)

def period(a):
    # Calculates the period of orbit
    return 2 * np.pi * sqrt(a**3 / mu)

def propellant(m0, delta_v, Isp):
    # Calculates the propellant mass required to perform a manouvre
    return m0 * (1 - np.exp(-delta_v / (Isp * g)))

def t_deorbit(T, Isp, m_prop):
    # Calculates the time to deorbit for a specific thrust and propellant mass
    return m_prop * Isp * g / T

#%%% Main program

if __name__ == "__main__" :
    
    r_i = 800 + r_earth
    r_f = 100 + r_earth
    a = (r_i + r_f) / 2
    delta_v = hohmann(r_i, r_f)
    
    m_prop = propellant(m, delta_v, Isp)
    
    T = np.linspace(100, 450, 450)
    #time = t_deorbit(T, Isp, m_prop)
    time = np.zeros(len(T))
    for i in range(len(T)):
        time[i] = t_deorbit(T[i], Isp, m_prop)
        
    fig, ax = plt.subplots()
    plt.plot(T, time)
    ax.set_ylabel('Time / s')
    ax.set_xlabel('Thrust / N')
    plt.show
    
    ''' Assuming hohmann transfer by thrusting at opposite ends of orbit starting from initial burn'''
    t_transfer = period(a) / 2
    #t_burn = time
    print(t_transfer/60) # minutes
    