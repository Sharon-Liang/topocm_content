#------This is the 2nd task of week-1 assignment-------------------------------
#import sys
#sys.path.append('./code')  #import modules in folder : code
#import pfaffian as pf
#import functions as fun

#import numpy as np
#import matplotlib as mpl
#import matplotlib.pyplot as plt

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
            sys[lat(2*x),lat(2*x + 1)] = t1
        for x in range(L//2 -1):
            sys[lat(2*x + 1),lat(2*x +2)] = t2
    else:
        for x in range(L//2):
            sys[lat(2*x),lat(2*x + 1)] = t1
            sys[lat(2*x + 1),lat(2*x +2)] = t2
    # Builder ensures hermiticity
    return sys
t1=1
t2=1
sys = ssh_chain(L=6)
kwant.plot(sys)
