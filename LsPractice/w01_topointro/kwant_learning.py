'''This is a note of learning kwant'''

'''Kwant is a free (open source) Python package for numerical calculations
   on tight-binding models with a strong focus on quantum transport. '''

'''The typical workflow with Kwant is as follows:
    1. Create an “empty” tight binding system.
    2. Set its matrix elements and hoppings.
    3. Attach leads (tight binding systems with translational symmetry).
       Pass the finalized system to a solver.'''
import kwant
#from types import SimpleNamespace
import matplotlib.pyplot as plt

''' Example-1 : Transport through a quantum wire '''
# Build a general tight-binding system
sys = kwant.Builder()
# Set the lattice constant and the lattice shape
a=1
lat = kwant.lattice.square(a) #defaut : a=1
# Set a set of parameters about the system
t = 1.0
W = 10
L = 30

# Define onsite interactions and hoopings
for i in range(L): #range() starts from 0
    for j in range(W):
        sys[lat(i,j)] = 4 * t

        if j>0:
            sys[lat(i,j),lat(i,j-1)] = -t
            # There is no need to define its Hermitian conjugate
        if i>0:
            sys[lat(i,j),lat(i-1,j)] = -t

#kwant.plot(sys)
'''Now we have a finite squre lattice system with size L * W '''

'''Now buil a 1D system'''
# Build a general tight-binding system
def chian(L):
    sys_chain = kwant.Builder()
    # Set the lattice constant and the lattice shape
    #a=1
    lat = kwant.lattice.chain() #defaut : a=1
    # Set a set of parameters about the system

    t = 1.0

    # Define onsite interactions and hoopings
    for j in range(L):
        sys_chain[lat(j)] = 4 * t
        if j>0:
            sys_chain[lat(j),lat(j-1)] = -t
    return sys_chain

syst = chian(10)
#kwant.plot(syst)


kwant.plot(syst)
