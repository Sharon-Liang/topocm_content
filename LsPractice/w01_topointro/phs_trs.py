#------This is the 1st task of week-1 assignment-------------------------------
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

#----- This is the example code of generating a spinful TRS Hamiltonian--------


# Define a function make random Hamiltonian
randn = np.random.randn    # generate a ramdon tensor

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

#------Task 1 : add TRS to a Halmitonian with p-h symmetry
'''
Guess: A Hamiltonian with TRS has Kramers degenercy, each En is doublly
       degenerate.
       In this case, there are two states at k=0 and k=π, so the energy
       level crossing won't change the systems topology
'''

# Generate a Hamiltonian with PHS : H = - sx @ H* @ sx
def make_random_phs_ham(N):
    if N % 2:
        raise ValueError('Matrix dimension should be a multiple of 2')
    # sy equals the direct production of a identity matrix(dim N//2) and the
    # Pauli matrix \sigma_y
    sx = np.kron(np.eye(N // 2), np.array([[0, 1], [1, 0]]))
    h = randn(N, N) + 1j * randn(N, N)
    h += h.T.conj()    # define a random Hamiltonian
    Th = - sx @ h.conj() @ sx    # @ is inner product, * is the product of
                               # corresponding elements. A(11)*B(11)
    return (h + Th) / 4


alphas = np.linspace(0, 1, 1000)    # evenly select 1000 points between [0,1]

# Find spectrum
def energies(alpha, h0, h1):
    h = (1 - alpha) * h0 + alpha * h1
    return np.linalg.eigvalsh(h)

def find_spectrum(alphas,h0,h1):
    spectrum = [energies(a,h0,h1) for a in alphas]
    return np.array(spectrum)


np.random.seed(102)   # make the random numbers predictable, i.e the same random
                    # number sequence of each specific seed

# Now start with a PHS Hamiltonian
H0 = make_random_phs_ham(10)
H1 = make_random_phs_ham(10)

# Energy spectrum with the change of α
spectrum = find_spectrum(alphas, H0, H1)

plt.title('Energy spectrum with the change of α')
plt.plot(alphas, spectrum[:,0:len(spectrum[0])], color = 'black')
plt.ylim(ymin = -1.5, ymax = 1.5)
plt.xlim(xmin = 0, xmax = 1)
plt.show()
