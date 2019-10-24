#------This is the 1st task of week-1 assignment-------------------------------
import sys
sys.path.append('./code')  #import modules in folder : code
import pfaffian as pf

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
#import pfaffian as pf     # Pfaffin package

#----- This is the example code of generating a spinful TRS Hamiltonian--------


# Define a function make random Hamiltonian
randn = np.random.randn    # generate a ramdon tensor
alphas = np.linspace(0, 1, 1000)    # evenly select 1000 points between [0,1]

def make_random_ham(N):
    H = randn(N, N) + 1j * randn(N, N)
    H += H.T.conj()  # arr.T : transpose of an array
    return H / 2

# Define a function make a Halmiltonian with spin-1/2 TRS
def make_random_symplectic_ham(N):
    if N % 2:
        raise ValueError('Matrix dimension should be a multiple of 2')
    # sy equals the direct production of a identity matrix(dim N//2) and the
    # Pauli matrix \sigma_y
    sy = np.kron(np.eye(N // 2), np.array([[0, -1j], [1j, 0]]))
    h = randn(N, N) + 1j * randn(N, N)
    h += h.T.conj()    # define a random Hamiltonian
    Th = sy @ h.conj() @ sy    # @ is inner product, * is the product of
                               # corresponding elements. A(11)*B(11)
    return (h + Th) / 4


'''Generate a Hamiltonian with PHS : H = - sx @ H* @ sx
    Note that this matrix is not a skrew matrix'''

def make_random_phs_ham(N):
    if N % 2:
        raise ValueError('Matrix dimension should be a multiple of 2')
    # sy equals the direct production of a identity matrix(dim N//2) and the
    # Pauli matrix \sigma_y
    ##sx = np.kron(np.eye(N // 2), np.array([[0, 1], [1, 0]]))
    sx = np.kron(np.array([[0, 1], [1, 0]]), np.eye(N // 2))
    h = randn(N, N) + 1j * randn(N, N)
    h += h.T.conj()    # define a random Hamiltonian
    Th = - sx @ h.conj() @ sx    # @ is inner product, * is the product of
                               # corresponding elements. A(11)*B(11)
    return (h + Th) / 4

def make_H_skrew(h):
    N = h.shape[0]
    if N % 2:
        raise ValueError('Matrix dimension should be a multiple of 2')
    p = np.kron(np.array([[1, 1], [1j, -1j]]), np.eye(N // 2))
    h = p @ h @ p.T.conj()
    return h/2


def make_BdG_ham(N):
    # This is antisymmetric basis
    H = 1j * randn(2*N, 2*N)
    H += H.T.conj()
    return H / 2

# Find spectrum
def energies(alpha, h0, h1):
    h = (1 - alpha) * h0 + alpha * h1
    return np.linalg.eigvalsh(h)

def find_spectrum(alphas,h0,h1):
    spectrum = [energies(a,h0,h1) for a in alphas]
    return np.array(spectrum)

# Phaffin
def find_pfaffian(alphas, h0, h1):
    """Function caculates the Pfaffian for a Hamiltonian.

    Parameters:
    -----------
    alphas : numpy array
        Range of alphas for which the energies are calculated.
    H0 : numpy array
        Hamiltonian, same size as H1.
    H1 : numpy array
        Hamiltonian, same size as H0.

    Returns:
    --------
    pfaffians : numpy array
        Pfaffians for each alpha.
    """
    def H(alpha):
        return (1 - alpha) * h0 + alpha * h1
    pfaffians = [np.sign(np.real(pf.pfaffian(1j*H(a)))) for a in alphas]
    return np.array(pfaffians)

#------Task 1 : add TRS to a Halmitonian with p-h symmetry
'''
Guess: A Hamiltonian with TRS has Kramers degenercy, each En is doublly
       degenerate.
       In this case, there are two states at k=0 and k=π, so the energy
       level crossing won't change the systems topology
'''

np.random.seed(102)   # make the random numbers predictable, i.e the same random
                    # number sequence of each specific seed

# Now start with a PHS Hamiltonian
H0 = make_random_phs_ham(8)
H1 = make_random_phs_ham(8)

H0 = make_H_skrew(H0)
H1 = make_H_skrew(H1)

pfaffian = find_pfaffian(alphas, H0, H1)
print(pfaffian.shape)

# Energy spectrum with the change of α
spectrum = find_spectrum(alphas, H0, H1)

fig1 = plt.figure()
ax1 = fig1.add_subplot(211)
plt.title('Energy spectrum as a function of α')
ax1.plot(alphas, spectrum[:,0:len(spectrum[0])], color = 'black')
plt.ylim(ymin = -1.5, ymax = 1.5)
plt.xlim(xmin = 0, xmax = 1)
plt.xlabel('α')
plt.ylabel('Energy')

ax2 = fig1.add_subplot(212)
ax2.plot(alphas, pfaffian, color = 'black')
plt.ylim(ymin = -1.5, ymax = 1.5)
plt.xlim(xmin = 0, xmax = 1)
plt.xlabel('α')
plt.ylabel('Pfaffin')

plt.show()
