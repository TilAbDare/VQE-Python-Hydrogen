# VQE-Python-Hydrogen
simulation of Hydrogen molecule to calculate ground state energy by python manually over bonding 
distance in the range of 0.2 to 2.8 angstrom.

It was implemented without any quantum computing libraries and just it was done by python
programming packages.

Only The hamiltonian form generated in terms of pauli matrices by qiskit library and afterwards 
the pauli matrices' coefficients were saved in a datasheet (I used excel) on local machine and by Pandas libraries
they were read by the program


mapping to qubits by: Parity transformation with two qubit reduction
number of pauli terms: 5
No. of Qubits: 2
ansatz and refrence wavefunction: UCCSD and HF, respectively. 
algorithm: VQE
Optimizer: SLSQP









it is worth mentioning that the preqiuistess to prepare this project was used via
a project by "jjgoings", https://github.com/jjgoings/jjgoings/commits?author=jjgoings
