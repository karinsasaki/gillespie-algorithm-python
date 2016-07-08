#!/usr/bin/env python2.7

#-----------------------------------------------------------------
# Gillespie algorithm for a simple decay reaction A -> 0 of rate k
#-----------------------------------------------------------------

import sys
import random as rdm              # pseudo-random numbers generator, random() generates random number between 0 and 1
from math import log, exp, floor  # basic math functions: exp(), log(), floor()
import numpy as np
import scipy.integrate as integrate # to solve ODEs

from matplotlib import pyplot as plt
from matplotlib import interactive
interactive(True)

# INITIATION
Na = 100            # Initial number of molecules
k  = 1.             # Reaction rate
FinalTime = 10.     # Simulation time

# SIMULATION RUN
def run():
    t = 0   # current time
    n = Na  # current number of molecules
    times = [t]
    nbMol = [n]
    while (t < FinalTime and n > 0):
        tau = - log(rdm.random()) / (n * k) 
        n = n-1
        t += tau
        times.append(t)
        nbMol.append(n)
    return [times,nbMol]
    
# DETERMINISTIC RUN: solves dn/dt = -kn (could write the analytical result, but embryonic basis to generalize later on)
def ode(n,t):
    return (- k*n)
    
def deterministic(ode,time):
    n = integrate.odeint(ode,Na,time)
    return n
    
# PLOT
def plot(data):
        plt.plot(data[0],data[1],linewidth=1.)        
        plt.show()

# KEY INPUT (Python 2.7 only...)
def newRunInput():
    newRun = False
    unknownKey = True
    if __name__ == '__main__':
        while (unknownKey == True):
            key = raw_input("New run ? (y/n)")
            if (key == "y"):
                newRun = True
                unknownKey = False
            elif (key == "n"):
                newRun = False
                unknownKey = False
    return newRun

# MAIN
newRun = True
nRun = 0
time = np.linspace(0, FinalTime, num=100)
nDeterministic = deterministic(ode,time)

if __name__ == '__main__':
    fig = plt.figure()  
    fontsize = 14
    plt.xlabel("Time", fontsize=fontsize)
    plt.ylabel("Molecules", fontsize=fontsize)
    plt.title("Deterministic solution and Gillespie simulations A->0", fontsize=fontsize)
    plt.plot(time, nDeterministic, "k--", linewidth=4.)
    plt.show()
    
while newRun:
      nRun += 1        
      result = run()
      plot(result)
      plt.plot(time, nDeterministic, "k--", linewidth=4.)
      newRun = newRunInput()              