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
lattice = LineLattice(num_nodes=N, boundary_condition = boundary_condition)
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

"""The following circuit must be derived from the GS by hand."""

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
#print(circuit)
#print()
#plt.show()


# Ground state
qubit_hamiltonian = JordanWignerMapper.mode_based_mapping(Hamiltonian)

job = pEstimator().run(circuit,qubit_hamiltonian)
result = job.result()
values = result.values
print('Ground State:',values[0])

omega = values[0]
