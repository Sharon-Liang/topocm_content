#-----Test the function of the codes-------------------------------------------
import sys
sys.path.append('./code')  #import modules in folder : code
import pfaffian as pf
import functions as fun

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from types import SimpleNamespace

import kwant

def ssh_chain(L=None, periodic=False):
    lat = kwant.lattice.chain()

    if L is None:
        syst = kwant.Builder(kwant.TranslationalSymmetry((-1,)))
        L = 1
    else:
        syst = kwant.Builder()

    for x in range(1, L//2):
        syst[lat(2 * x - 1),lat(2 * x)] = t1
        syst[lat(2 * x),lat(2 * x + 1)] = t2
    # Builder ensures hermiticity
    return syst

    sys = ssh_chain(L=25)
    p = SimpleNamespace(t1=1, t2=1)
    fun.spectrum(sys,p)
