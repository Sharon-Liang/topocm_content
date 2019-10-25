#------This is the 2nd task of week-1 assignment-------------------------------
import sys
sys.path.append('./code')  #import modules in folder : code
import pfaffian as pf
import functions as fun

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

import kwant
'''Kwant is a free (open source) Python package for numerical calculations
   on tight-binding models with a strong focus on quantum transport. '''

'''Example: Build 1D Kitaev chain'''
def kitaev_chain(L=None, periodic=False):
    lat = kwant.lattice.chain()
         # make 1-D chain, defaut parameters: a=1, name='', norbs=None
         # a is the lattice parameter
         # the 1D chain only have 1 parameter to label lattice sites, lat(i)

    if L is None:
        syst = kwant.Builder(kwant.TranslationalSymmetry((-1,)))
        # kwant.Builder([symmetry]) : A tight binding system defined on a graph
        # kwant.TranslationalSymmetry((-1,)) : build a system with period -1
        L = 1
        # if L is not assigned, then L = 1
    else:
        syst = kwant.Builder()
        # syst defines a tight binding model, no specific symmetry assigned.

    # transformation to antisymmetric basis
    U = np.array([[1.0, 1.0], [1.j, -1.j]]) / np.sqrt(2)

    def onsite(onsite, p):  #onsite interactions
        return - p.mu * U @ pauli.sz @ U.T.conj()
        ''' p : SimpleNamespace object
                A container used to store Hamiltonian parameters.
                The parameters that are sequences are used as plot axes.'''

    for x in range(L):
        syst[lat(x)] = onsite
        # sys[lat(i)] sets the onsite energy for point i

    def hop(site1, site2, p):
        return U @ (-p.t * pauli.sz - 1j * p.delta * pauli.sy) @ U.T.conj()

    syst[kwant.HoppingKind((1,), lat)] = hop
    '''kwant.builder.HoppingKind(delta, family_a, family_b=None)
       no family_b means family_b is the same kind as family_a'''

    if periodic:
        def last_hop(site1, site2, p):
            return hop(site1, site2, p) * (1 - 2 * p.lambda_)

        syst[lat(0), lat(L - 1)] = last_hop
    return syst

'''Now write SSH model Hamiltonain'''
def ssh_chain(L=None, periodic=False):
    lat = kwant.lattice.chain()

    if L is None:
        syst = kwant.Builder(kwant.TranslationalSymmetry((-1,)))
        L = 1
    else:
        syst = kwant.Builder()

    for x in range(L//2):
        syst[lat(2 * x - 1),lat(2 * x)] = p.t1
        syst[lat(2 * x),lat(2 * x + 1)] = p.t2
    # Builder ensures hermiticity
    return syst

    sys = ssh_chain(L=25)
    kwant.plotter.plot(sys)
