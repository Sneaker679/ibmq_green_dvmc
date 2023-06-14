#Packages----------------------------------------------------

import numpy as np
import matplotlib.pyplot as plt
from parameters import *

from qiskit.primitives import Estimator as pEstimator
from qiskit_nature.second_q.mappers import JordanWignerMapper
from qiskit_nature.second_q.operators import FermionicOp
from qiskit import QuantumCircuit,QuantumRegister
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
import copy
import time

#------------------------------------------------------------


# Building the Hamiltonian of Hubbard's model using qiskit

    # Initializing latice
    
boundary_condition = BoundaryCondition.OPEN
lattice = SquareLattice(rows=2, cols=2, boundary_condition = boundary_condition)
#lattice.draw()
#plt.show()


    # Initializing the Hubbard Hamiltonian using the lattice

Hamiltonian = FermiHubbardModel(
    lattice.uniform_parameters(
        uniform_interaction = t,
        uniform_onsite_potential = -mu,
    ),
    onsite_interaction = U,
).second_q_op() # The last line transforms the newly created Hamiltonian into a FermionicOp object


# Circuit
vec = np.linalg.eigh(Hamiltonian.to_matrix().toarray())[1][:,0].real.tolist()

q = QuantumRegister(8)

circuit = QuantumCircuit(q)

circuit.initialize(vec, [q[0],q[1],q[2],q[3],q[4],q[5],q[6],q[7]])
circuit.draw()


# Ground state
qubit_hamiltonian = JordanWignerMapper.mode_based_mapping(Hamiltonian)

job = pEstimator().run(circuit,qubit_hamiltonian)
result = job.result()
values = result.values
print('Ground State:',values[0])

omega = values[0]
