#!/usr/bin/env python2.7

#-----------------------------------------------------------------
# Gillespie algorithm for a binary equilibrium reaction:
# A -> B rate k_AB
# B -> A rate k_BA
#-----------------------------------------------------------------

import sys
import random as rdm       # pseudo-random numbers generator, random() generates random number between 0 and 1
from math import log, exp  # basic C math functions: exp(), log()
import numpy as np
import scipy.integrate as integrate # to solve ODEs


from matplotlib import pyplot as plt
from matplotlib import interactive
interactive(True)

# STOCHASTIC RUN
def run(params):
    
    # would be better with a general input format & parser...
    nA = params[0]  
    nB = params[1]
    kAB = params[2]
    kBA = params[3]

    t = 0   # current time
    times = [t]
    nbMolA = [nA]
    nbMolB = [nB]
    
    while (t < FinalTime):
        # random numbers generation
        r1 = rdm.random()
        r2 = rdm.random()
        # calculate next reaction time t+\tau
        k1 = nA * kAB
        k2 = nB * kBA       
        K  = k1 + k2
        tau = - log(r1) / K
        t += tau
        # find and perform next reaction        
        if r2 * K < k1:
            nA -= 1
            nB += 1
        else:
            nA += 1
            nB -= 1
        
        times.append(t)
        nbMolA.append(nA)
        nbMolB.append(nB)
        
    return [times,nbMolA, nbMolB]


# DETERMINISTIC RUN: redundant and unpractical for now, but to generalize into dn/dt = An with n vector of solutes and A reaction matrix
def deterministic(params,time):
    ni = np.array([params[0],params[1]])
    kAB = params[2]
    kBA = params[3]
    def ode(n,t):
        dn0 = -kAB*n[0] + kBA*n[1]
        print dn0
        dn1 =  kAB*n[0] - kBA*n[1]
        print dn1
        return [dn0, dn1]
        
    n = integrate.odeint(ode,ni,time)
    return n

# PLOT FUNCTION
def plotf(resultAB):
        plt.plot(resultAB[0],resultAB[1],'b-')        
        plt.plot(resultAB[0],resultAB[2],'r-')        
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
# Initiation
 # Initial number of molecules
Na = 100
Nb = 0            
 # Reaction rates
kAB  = 1              
kBA  = 1
 # Simulation time
FinalTime = 10     
params = [Na,Nb,kAB,kBA]
time = np.linspace(0, FinalTime, num=100)
nDeterministic = deterministic(params,time)

# Plotting
if __name__ == '__main__':
    fig = plt.figure()  
    fontsize = 14
    plt.xlabel("Time", fontsize=fontsize)
    plt.ylabel("Molecules", fontsize=fontsize)
    plt.title("Gillespie simulations A->0", fontsize=fontsize)
    plt.plot(time,nDeterministic[:,0], 'darkblue', dashes=[3,0,3], linewidth=3.)
    plt.plot(time,nDeterministic[:,1], 'darkred', dashes=[3,0,3], linewidth=3.)
    plt.show()

# Main loop
newRun = True
nRun = 0    
while newRun:
    nRun += 1        
    result = run(params)
    plotf(result)
    newRun = newRunInput()
    plt.plot(time,nDeterministic[:,0], 'darkblue', dashes=[3,0,3], linewidth=3.)
    plt.plot(time,nDeterministic[:,1], 'darkred', dashes=[3,0,3], linewidth=3.)              