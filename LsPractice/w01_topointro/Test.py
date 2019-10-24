import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

randn = np.random.randn
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

H1 = make_random_phs_ham(4)
H2 = make_H_skrew(H1)
H2 += H2.T
print(H2)
