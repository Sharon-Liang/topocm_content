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
        sys[lat(x)] = 0

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



'''First calculate the energy spectrum change with respect to if t2/t1.'''
'''[t1=0,t2=1] and [t1=1,t2 = 0] should be degenerate, since they are simply
   several pairs of dimmers'''
# Now define a infinite SSH chain, and check it's spectrum
# Ek = -2 * [t1 + t2 cos(2ka)] , defaut a=1
def SSH_inf_chain(L):
    lat = kwant.lattice.chain()
    sys= kwant.Builder(kwant.TranslationalSymmetry((2,)))
    # Define onsite energy
    for x in range(L):
        sys[lat(x)] = 0

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

t1 = 1
t2 = 0.5
sys_inf = SSH_inf_chain(L=10).finalized()
kwant.plotter.bands(sys_inf)









'''Connect two chains with t1>>t2 and t2>>t1 respectively'''
