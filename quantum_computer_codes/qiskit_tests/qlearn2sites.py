#import pauli as p
import numpy as np
import matplotlib.pyplot as plt
from qiskit.primitives import Estimator as pEstimator
from qiskit_nature.second_q.mappers import JordanWignerMapper,QubitConverter
from qiskit_nature.second_q.operators import FermionicOp
from qiskit import QuantumCircuit
from qiskit import transpile
from qiskit_nature.second_q.hamiltonians import FermiHubbardModel
from qiskit_nature.second_q.hamiltonians.lattices import (
    BoundaryCondition,
    HyperCubicLattice,
    Lattice,
    LatticeDrawStyle,
    LineLattice,
    SquareLattice,
    TriangularLattice,
)
import sys

np.set_printoptions(threshold=sys.maxsize, linewidth=800)

# Building the Hamiltonian of Hubbard's model

N = 2 # 2 sites
t = -1
U = 4
mu = 2 #U/2

    # Initializing latice
    
boundary_condition = BoundaryCondition.OPEN
lattice = LineLattice(num_nodes = N, boundary_condition = boundary_condition)
#lattice.draw()
#plt.show()

    # Initializing the Hubbard Hamiltonian

Hamiltonian = FermiHubbardModel(
    lattice.uniform_parameters(
        uniform_interaction = t,
        uniform_onsite_potential = -mu,
    ),
    onsite_interaction = U,
).second_q_op()


# Exact result

Hamiltonian_matrix = Hamiltonian.to_matrix().toarray() # Transform H into matrix.
print(np.real(Hamiltonian_matrix))
print(Hamiltonian_matrix.shape)


'''Fetches eigen vectors and values. 
E are ground_state energies,S are Ground_State vectors.'''
E,S = np.linalg.eigh(Hamiltonian_matrix)

print()
print('Exact ground state energy:',E[0])
print('Exact ground state:',S[:,0])
print()


# Circuit

'''We *know* that the ground state is in the |5> |6> |9> and |10> block. Hence, we can derive from the bloc what the circuit should look like, and we obtain the following circuit:'''

'''The text before might not be true after all. The actual circuits that works
is in test.py.'''

'''
theta = 0.7854074074074073

circuit = QuantumCircuit(2*N) # Creating circuit with 4 qubits.
circuit.ry(theta, 2)
circuit.h(0)
circuit.cx(2,3)
circuit.cx(0,1)
circuit.x(2)
circuit.x(0)
circuit.cx(1,3)
circuit.cx(1,2)
circuit.swap(0,1)
circuit.cz(1,2)
circuit.swap(1,2)
'''

theta = 2*1.178097245

circuit = QuantumCircuit(2*N) # Creating circuit with 4 qubits.
circuit.ry(theta, 2)
circuit.h(0)
circuit.cx(2,3)
circuit.cx(0,1)
circuit.x(2)
circuit.x(0)
circuit.cx(1,3)
circuit.cx(1,2)
circuit.cz(1,2)
circuit.swap(1,2)

circuit.draw('mpl')
print(circuit)
print()
#plt.show()

qubit_converter = QubitConverter(JordanWignerMapper()) # Defines transformations as the Jordan-Wigner transformation.
HHam = qubit_converter.convert(Hamiltonian) # Converts the Hamiltonian "operator".
'''print(dir(HHam))
print(type(HHam))
print(HHam)
exit()'''

estimator = pEstimator() # Evaluates mean value of observable
job = estimator.run(circuit,HHam) # Runs the simulation of the quantum computer
result = job.result()
values = result.values # Transforms this object into a list that we can later index.

print('Result of simulation:', values[0])
