# --------------------------------------------------------------------
# ******************  Importing Packages *****************************
# --------------------------------------------------------------------
from __future__ import print_function  # Py 2.6+; In Py 3k not needed

import numpy as np
from scipy.linalg import block_diag
from scipy.optimize import minimize
import timeit
import pandas as pd
import matplotlib.pyplot as plt

# --------------------------------------------------------------------
# ************************  Quantum Gate *****************************
# --------------------------------------------------------------------
start = timeit.default_timer()

np.set_printoptions(precision=4, suppress=True)

# Pauli matrices

I = np.array([[1, 0],
              [0, 1]])

Sx = np.array([[0, 1],
               [1, 0]])

Sy = np.array([[0, -1j],
               [1j, 0]])

Sz = np.array([[1, 0],
               [0, -1]])

# Hadamard matrix

H = (1 / np.sqrt(2)) * np.array([[1, 1],
                                 [1, -1]])

# Phase matrix
S = np.array([[1, 0],
              [0, 1j]])

# single qubit basis states |0> and |1>

q0 = np.array([[1],
               [0]])

q1 = np.array([[0],
               [1]])

# Projection matrices |0><0| and |1><1|
P0 = np.dot(q0, q0.conj().T)
P1 = np.dot(q1, q1.conj().T)

# Rotation matrices as a function of theta, e.g. Rx(theta), etc.
Rx = lambda theta: np.array([[np.cos(theta / 2), -1j * np.sin(theta / 2)],
                             [-1j * np.sin(theta / 2), np.cos(theta / 2)]])

Ry = lambda theta: np.array([[np.cos(theta / 2), -np.sin(theta / 2)],
                             [np.sin(theta / 2), np.cos(theta / 2)]])

Rz = lambda theta: np.array([[np.exp(-1j * theta / 2), 0.0],
                             [0.0, np.exp(1j * theta / 2)]])

# CNOTij, where i is the control qubit and j is target qubit

CNOT10 = np.kron(P0, I) + np.kron(P1, Sx)  # control -> q1, target -> q0
CNOT01 = np.kron(I, P0) + np.kron(Sx, P1)  # control -> q0, target -> q1

SWAP = block_diag(1, Sx, 1)

# --------------------------------------------------------------------
# ***********************  initial basis ************************
# --------------------------------------------------------------------
# initial basis, put in |01> state with Sx operator on q0

psi0 = np.zeros((4, 1))
psi0[0] = 1
psi0 = np.dot(np.kron(Sx, I), psi0)

# --------------------------------------------------------------------
# ***********************  Ansatz Preparation ************************
# --------------------------------------------------------------------

# read right-to-left (bottom-to-top?)

ansatz = lambda theta: (np.dot(np.dot(np.kron(-Ry(np.pi / 2), Rx(np.pi / 2)),
                                      np.dot(CNOT10,
                                             np.dot(np.kron(I, Rz(theta)),
                                                    CNOT10))),
                               np.kron(Ry(np.pi / 2), -Rx(np.pi / 2))))

# --------------------------------------------------------------------
# *****************  Parallel Energy Calculation *********************
# --------------------------------------------------------------------

def projective_expected(theta, ansatz, psi0):
    # this will depend on the hard-coded Hamiltonian + coefficients
    circuit = ansatz(theta[0])
    psi = np.dot(circuit, psi0)

    # for 2 qubits, assume we can only take Pauli Sz measurements (Sz \otimes I)
    # we just apply the right unitary for the desired Pauli measurement
    measureZ = lambda U: np.dot(np.conj(U).T, np.dot(np.kron(Sz, I), U))

    energy = 0.0

    # although the paper indexes the hamiltonian left-to-right (0-to-1)
    # qubit-1 is always the top qubit for us, so the tensor pdt changes
    # e.g. compare with the "exact Hamiltonian" we explicitly diagonalized

    # <I1 I0>
    energy += g0  # it is a constant

    # <I1 Sz0>
    U = SWAP
    energy += g1 * np.dot(psi.conj().T, np.dot(measureZ(U), psi))

    # <Sz1 I0>
    U = np.kron(I, I)
    energy += g2 * np.dot(psi.conj().T, np.dot(measureZ(U), psi))

    # <Sz1 Sz0>
    U = CNOT01
    energy += g3 * np.dot(psi.conj().T, np.dot(measureZ(U), psi))

    # <Sx1 Sx0>
    U = np.dot(CNOT01, np.kron(H, H))
    energy += g4 * np.dot(psi.conj().T, np.dot(measureZ(U), psi))

    return np.real(energy)[0, 0]


# --------------------------------------------------------------------
# ************  Energy Calculation by Python Manualy *************
# --------------------------------------------------------------------

db = pd.read_excel('/Users/kourosh/Desktop/Demo/Demo_H2/Demo_Full.xlsx')

energy_py_ucc = []
dist = []

g0_list = db['g0']
g1_list = db['g1']
g2_list = db['g2']
g3_list = db['g3']
g4_list = db['g4']
d_list = db['d']
nrp_list = db['NR']

for i in range(len(g0_list)):
    g0 = float(g0_list[i])
    g1 = float(g1_list[i])
    g2 = float(g2_list[i])
    g3 = float(g3_list[i])
    g4 = float(g4_list[i])
    nrp = float(nrp_list[i])
    theta = [0.0]
    result = minimize(projective_expected, theta, method='SLSQP', args=(ansatz, psi0))
    theta = result.x[0]
    val = result.fun
    final_res = val + nrp
    energy_py_ucc += [final_res]
    dist.append(d_list[i])



stop = timeit.default_timer()
runtime = stop - start
print('Run Time: ', runtime, 'sec', ' or ', runtime/60, 'min','\n')
# --------------------------------------------------------------------
# ********************* Saving The Data locally **********************
# --------------------------------------------------------------------
df = pd.DataFrame(list(zip(dist, energy_py_ucc)),
                  columns = ['distance', 'Python Manually'])
#df.to_excel("Diag.xlsx")
print('-'*8, ' Final Results ', '-'*8)
print(df)


# --------------------------------------------------------------------
# ******************************  Plot *******************************
# --------------------------------------------------------------------

plt.plot(dist, energy_py_ucc, '|', color='Red', label='Py_UCC_Parity')
plt.grid(True, linestyle='-.', linewidth=0.5, which='major')
plt.title("Ground State Energy Curve of hydrogen Molecule")
plt.xlabel("H-H Distance ($\AA$) ")
plt.ylabel("Energy (Hartree)")
plt.legend()
plt.show()

