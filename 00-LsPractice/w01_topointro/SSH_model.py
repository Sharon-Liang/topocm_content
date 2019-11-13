#------This is the 2nd task of week-1 assignment-------------------------------
'''Mind the spaces in python !!! '''
import sys
sys.path.append('./code')  #import modules in folder : code
import pfaffian as pf

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

import kwant
'''Kwant is a free (open source) Python package for numerical calculations
   on tight-binding models with a strong focus on quantum transport. '''

'''Now write SSH model Hamiltonain'''
def ssh_chain(L=None, periodic=False):
    lat = kwant.lattice.chain()
    # with defaut a=1

    if L is None:
        sys = kwant.Builder(kwant.TranslationalSymmetry((-1,)))
        L = 1
    else:
        sys = kwant.Builder()


    # Define onsite energy
    for x in range(L):
        if x%2 == 0:
            sys[lat(x)] = ua
        else:
            sys[lat(x)] = ub

    # Define hooping terms
    if L%2 == 0:
        for x in range(L//2):
            sys[lat(2*x),lat(2*x + 1)] = -t1
        for x in range(L//2 -1):
            sys[lat(2*x + 1),lat(2*x +2)] = -t2
    else:
        for x in range(L//2):
            sys[lat(2*x),lat(2*x + 1)] = -t1
            sys[lat(2*x + 1),lat(2*x +2)] = -t2
            # Builder ensures hermiticity

    if periodic:
        if L%2 == 0:
            t_last = -t2
        else:
            t_last = -t1
        sys[lat(0), lat(L - 1)] = t_last
    return sys


'''First calculate the energy spectrum in k-space.'''
'''[t1=0,t2=1] and [t1=1,t2 = 0] should be degenerate, since they are simply
   several pairs of dimmers'''
# Now define a infinite SSH chain, and check it's spectrum
# Ek = -2 * [t1 + t2 cos(2ka)] , defaut a=1
def SSH_inf_chain(L):
    lat = kwant.lattice.chain()
    sys= kwant.Builder(kwant.TranslationalSymmetry((2,)))
    # Define onsite energy
    for x in range(L):
        sys[lat(x)] = 1

    # Define hooping terms
    if L%2 == 0:
        for x in range(L//2):
            sys[lat(2*x),lat(2*x + 1)] = -t1
        for x in range(L//2 -1):
            sys[lat(2*x + 1),lat(2*x +2)] = -t2
    else:
        for x in range(L//2):
            sys[lat(2*x),lat(2*x + 1)] = -t1
            sys[lat(2*x + 1),lat(2*x +2)] = -t2
    return sys

'''Now calculte the energy spectrum as a function of t2/t1 '''
t1 = 1 # Set as default
ua = ub = 0
N = 30 # Number of sites
       # If N=odd, there always be a zero mode
ts = np.linspace(0,2.5,1000)
energies = np.zeros((1000,N))

i=0
for t2 in ts:
    sys = ssh_chain(L=N).finalized()
    ham = sys.hamiltonian_submatrix()
    energies[i,:] = np.linalg.eigvalsh(ham)
    i = i + 1

fig1 = plt.figure()
ax1 = fig1.add_subplot(121)
plt.title('ua = ub =0')
ax1.plot(ts, energies[:,0:N-1], color = 'black')
#plt.ylim(ymin = -10, ymax = 10)
plt.xlim(xmin = 0, xmax = 2.5)
plt.xlabel('t2 / t1')
plt.ylabel('Energy')


'''Now add a sublattice potential '''
t1 = 1 # Set as default
ua = 0.2
ub = -0.2
N = 30 # Number of sites
       # If N=odd, there always be a zero mode
ts = np.linspace(0,2.5,1000)
energies = np.zeros((1000,N))

i=0
for t2 in ts:
    sys = ssh_chain(L=N).finalized()
    ham = sys.hamiltonian_submatrix()
    energies[i,:] = np.linalg.eigvalsh(ham)
    i = i + 1

ax2 = fig1.add_subplot(122)
plt.title('ua=0.2, ub=-0.2')
ax2.plot(ts, energies[:,0:N-1], color = 'black')
#plt.ylim(ymin = -10, ymax = 10)
plt.xlim(xmin = 0, xmax = 2.5)
plt.xlabel('t2 / t1')
#plt.ylabel('Energy')
plt.show()


'''Calculate the wave function of the zero modes'''
# zero modes: the L/2-1 and L/2 column of the eigen vectors
# 1st excited state: L/2 + 1 column of the eigen vectors
t1 = 1
t2 = 2
ua = ub = 0
N = 30
sys = ssh_chain(L=N).finalized()
ham = sys.hamiltonian_submatrix()
ev, evec = np.linalg.eigh(ham)
zero1 = evec[:,N//2].copy()
zero1 = pow(abs(zero1),2)
plt.bar(x=np.arange(1,31),height = zero1[:])
plt.show()


'''Connect two chains with t1>>t2 and t2>>t1 respectively'''


'''Observations:
   1. t2/t1 ~ 1.5, the zero mode emerges
   2. add sublattice chemical potential can destory the
      topological non-trivial phase.
   '''
