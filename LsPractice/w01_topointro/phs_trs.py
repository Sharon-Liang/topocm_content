#------This is the 1st task of week-1 assignment-------------------------------
import numpy as np
#----- This is the example code of generating a spinful TRS Hamiltonian--------

np.random.seed(1)   # make the random numbers predictable, i.e the same random
                    # number sequence of each specific seed

# Define a function make random Hamiltonian
randn = np.random.randn    # generate a ramdon tensor

def make_random_ham(N):
    H = randn(N, N) + 1j * randn(N, N)
    H += H.T.conj()  # arr.T : transpose of an array
    return H / 2
    
