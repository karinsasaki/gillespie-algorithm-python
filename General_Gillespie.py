#!/usr/bin/env python2.7
# 
# GILLESPIE
# Herve Turlier - March 2016
#
#-----------------------------------------------------------------
# GILLESPIE ALGORITHM FOR A GENERAL SYSTEM OF REACTIONS:
#
# The reaction i of rate ki, involves species j with signed stoechiometric coefficients v_ij
#
#      reaction 1:  v_A1*A1 + v_B1*B1 + ... + v_F1*F1 --k1--> 0 
#      reaction 2:  v_A2*A2 + v_B2*B2 + ... + v_F2*F2 --k2--> 0
#      ...
#      reaction n:  v_An*An + v_Bn*Bn + ... + v_Fn*Fn --kn--> 0
#
# Reaction matrix R = (Stoechiometric matrix S, rates vector K)  
# with Stoechiometric matrix S = transpose of classical definition 
# -> size nx(m-1), with m-1=nb metabolites & n=nb reactions
#
#            metabolites          
#      ( v_A1 v_B1 ... v_F1 ) r
#      ( v_A2 v_B2 ... v_F2 ) e
# S =  ( ...  ...  ...  ... ) a
#      ( ...  ...  ...  ... ) c
#      ( v_An v_Bn ... v_Fn ) t
#
# Reaction rates (column) vector R (size n):
#   
#     ( k_1 )
#     ( k_2 )
# K = ( ... )
#     ( ... )
#     ( k_n )
#
# Reaction matrix (size nxm ):
#
#            metabolites          
#      r ( v_A1 v_B1 ... v_F1  k1  )
#      e ( v_A2 v_B2 ... v_F2  k2  )
# R =  a ( ...  ...  ...  ...  ... )
#      c ( ...  ...  ...  ...  ... )
#      t ( v_An v_Bn ... v_Fn  kn  )

# Metabolites number (line) vector N (size m):
#
# N = ( N_A, N_B, ..., N_F)
#
#-----------------------------------------------------------------
# TO DO:
#
# - Parser analyzing a list of reactions like: 
#       A + B -> C, k1
#       C -> 0, k2
#       2B -> D, k3 
# 
# - Rewrite the main function to get an input reactions file to parse.
#
# - Modify plot function to control colors & get a legend/metabolite
#
#-----------------------------------------------------------------
# LIBRARIES IMPORT
import sys, math, random, operator

# for matrix manipulation
try:
    import numpy as np
except ImportError:
    print("  Error: could not load numpy in python " + sys.version)
    sys.exit()

# for solving deterministic ODEs
try:
    import scipy.integrate as integrate
except ImportError:
    print("  Error: could not load scipy in python " + sys.version)
    sys.exit()

# for plotting
try:
    import matplotlib
except ImportError:
    print("  Error: could not load matplotlib in python " + sys.version)
    sys.exit()
    
import matplotlib.pyplot as plt
from matplotlib import interactive
interactive(True)
#-----------------------------------------------------------------


#-----------------------------------------------------------------
# USEFUL add-on FUNCTIONS

# function to send an error
def error(msg): # msg should be a string
    sys.stderr.write("  Error:" + msg)
    sys.exit()

# Home-made product function (return 1 if list is empty)
from functools import reduce
def prod(factors):
    return reduce(operator.mul, factors, 1)

# Function to add vector Nt (at time t) to matrix of metabolites number evolution Nmatrix
def add_line_to_Nmatrix(Nmatrix,Nt):
    if (len(Nmatrix) > 0):
        assert( len(Nt) == len(Nmatrix) )
        for i in range(len(Nt)):
            Nmatrix[i].append(Nt[i])
    else:
        for i in range(len(Nt)):
            Nmatrix.append([Nt[i]])

# Key input (Python 2.7 only...)
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

# Function to plot results of a run
def plot_run(run_results):
        time = run_results[0]
        Nmatrix = run_results[1]
        assert(len(time) == len(Nmatrix[0]))
        
        for i in range(len(Nmatrix)):
            plt.plot(time,Nmatrix[:][i])
        plt.show()

#-----------------------------------------------------------------
# PARSER of parameters file
 


#-----------------------------------------------------------------
# STOCHASTIC STEP from reaction matrix R and current metabolites vector N
# return new reaction time (immutable type) and update the metabolites number vector N (mutable type)
def time_step(R,N,t):
    
    n = len(R)     # Nb of lines = nb of reactiosn 
    m = len(R[0])  # Nb of columns = nb of metabolies + 1
    
    # Compute the current vector of propensities alphas
    #   For a reaction i v_iA*A + v_iB*B + ... + v_iF*F --ki--> 0 propensity is defined as 
    #       alpha_i = k_i * product_j ( N_j! / (N_j-v_ij')! )
    #       where j=A,B,...,F and v_ij' = abs( min(0,v_ij) ) >= 0
    alphas = []
    # loop on the reactions
    for i in range(n): 
        alpha_i = R[i][m-1] # rate of reaction i
        # loop on the metabolite to fill NList_i with N_j x (N_j-1) x ... x (N_j-v_ij) for all j such that v_ij < 0
        NList_i = []
        for j in range(m-1):
            v_ij = abs( min(0,R[i][j]) )
            if (v_ij > 0.):
                # check that stoechiometric coefficient is integer type
                if not( isinstance(v_ij, (int, long)) ):
                    error("Reaction matrix contains non-integer stoechiometric coefficients")
                for k in range(v_ij):
                    NList_i.append( N[j]-k )
        # the propensity for reaction i is obtained as product of NList_i elements with rate of reaction i
        alpha_i *= prod(NList_i)
        alphas.append( alpha_i ) 
     
    # Calculate the sum of all propensities
    alpha0 = sum(alphas) 
    
    # Set reaction time t+\tau
    r1 = random.random() # random number between 0 and 1
    tau = - math.log(r1) / alpha0
    t += tau

    # Set next reaction (= column index i in the reaction matrix R)
    r2 = random.random() # random number between 0 and 1
    x = alpha0 * r2
    i = -1
    while x > 0:
        x -= alphas[i+1]
        i += 1
        
    # Update the number of molecules with the reaction j
    for j in range(m-1):
        N[j] += R[i][j]
            
    # Return time (float is an immutable type)
    return t

#-----------------------------------------------------------------
# STOCHASTIC RUN from stoechiometric matrix S, rates vector R, metabolites vector N
# run the Gillespie algorithm for time from 0 to tmax, returns tuple with time vector + matrix of metabolites number evolution

def stochastic_run(R,Ninit,tmax):
    
    n = len(R)      # Nb of lines = nb of metabolites 
    m = len(R[0])   # Nb of columns = nb of reactions
    if ( len(Ninit) != m-1 ):
        error("Number of metabolites in the reaction matrix and initial vector do not match")
    if (tmax < 0):
        error("Final simulation time shall be positive")    
    
    # initiate time
    t = 0    
    times = [t]
    
    # current list of metabolites number 
    N = Ninit[:]
    
    # initiate a matrix of metabolite number evolution
    Nmatrix = []
    add_line_to_Nmatrix(Nmatrix,N)
    
    while (t < tmax):
        t = time_step(R,N,t)
        times.append(t)
        add_line_to_Nmatrix(Nmatrix,N)
        
    run_results = (times,Nmatrix) # tuple
    return run_results

#-----------------------------------------------------------------
# DETERMINISTIC RUN:

def deterministic_run(R,Ninit,times):

    N0 = Ninit[:]
    n = len(R)     # Nb of lines = nb of reactions
    m = len(R[0])  # Nb of columns = nb of metabolies + 1
    assert( n > 0 )
    assert( m > 1 )
    assert( len(N0) == m-1 )

    # system of ode(s)
    def ode(N,t):
        dN = []
        # loop on the metabolite to write dN_j/dt = ... for each metabolite j
        for j in range(m-1):
            dN_j = 0
            # loop on the reactions to find those involving metabolite j
            for i in range(n):
                dN_j_i = R[i][j] * R[i][m-1] # np.sign(x) is 0 if x=0, -1 if x negative, +1 if x positive
                # loop on the metabolites k to find which ones react in reaction i
                for k in range(m-1):
                    if ( np.sign(R[i][k]) < 0 ):
                        dN_j_i *=  math.pow( N[k], abs(R[i][k]) )
                dN_j += dN_j_i
            dN.append(dN_j)
        return dN
    
    ## solves the system of equation ode, with initial conditions Ninit
    Nsol = integrate.odeint(ode,N0,times)
   
    # convert solution into same format as for a stochastic run 
    Nmatrix = []
    for i in range(len(times)):
        add_line_to_Nmatrix(Nmatrix,Nsol[i])
            
    run_results = (times,Nmatrix) # tuple
    return run_results

#-----------------------------------------------------------------

# MAIN
# Initiation
 # Initial number of molecules
Ninit = [100,0]
 # Reaction matrix
R = [[-1,2,0.5],[0,-1,1.]]
 # Simulation time
tmax = 10
 # Time range from 0 to tmax (for deterministic run)
times = np.linspace(0, tmax, num=1000)
det_run_results = deterministic_run(R,Ninit,times)

# Plotting
if __name__ == '__main__':
    fig = plt.figure()  
    fontsize = 14
    plt.xlabel("Time", fontsize=fontsize)
    plt.ylabel("Molecules", fontsize=fontsize)
    plt.title("Gillespie simulations", fontsize=fontsize)
    plt.show()

# Main loop
newRun = True
nRun = 0 
while newRun:
    N = Ninit[:] 
    stoch_run_results = stochastic_run(R,N,tmax)
    plot_run(stoch_run_results)
    plot_run(det_run_results)
    nRun += 1
    newRun = newRunInput()