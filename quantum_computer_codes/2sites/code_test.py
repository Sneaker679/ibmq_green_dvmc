import sys
import numpy as np
import focktest_classes as fc
import excitation_matrices2 as em
from qiskit.quantum_info import Pauli,Operator
from qiskit.primitives import Estimator as pEstimator
from qiskit_nature.second_q.mappers import JordanWignerMapper
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

# Which Test?
answer = input('Tests to be ran?  ')


# Parameters
types = ['H+','H-','S+','S-']
N_values = [2]
N_min_values = [1]
spins_left = ['+','-']
spins_right = ['+','-']
t_values = [1]
U_values = range(9)
mu_values = range(5)
excit_doc = 'excitation.def'

#Test 1 - H/S matrix with fock and qc
if '1' in answer:
    t_q3 = -1
    t3 = 1
    U3 = 4
    mu3 = 2
    N3 = 2
    N_min3 = 1
    boundary_condition = BoundaryCondition.OPEN
    lattice = LineLattice(num_nodes = N3, boundary_condition = boundary_condition)

    Hamiltonian = FermiHubbardModel(
        lattice.uniform_parameters(
            uniform_interaction = t_q3,
            uniform_onsite_potential = -mu3,
        ),
        onsite_interaction = U3,
    ).second_q_op()
    
    theta = 2*1.178097245
    circuit = QuantumCircuit(2*N3) # Creating circuit with 4 qubits.
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
    
    for type in types:
        for spin_left in spins_left:
            for spin_right in spins_left:
                result1 = fc.matrix(type,excit_doc,N3,N_min3,spin_left,spin_right,t3,U3,mu3)
                result2 = em.matrix(type,excit_doc,N3,N_min3,spin_left,spin_right,Hamiltonian,circuit)
                                
                print(['type','spin_left','spin_right','N','N_min','t','U','mu'])
                print([type,spin_left,spin_right,N3,N_min3,t3,U3,mu3])
                bool = np.allclose(result1, result2, rtol=1e-05, atol=1e-05, equal_nan=False)
                print('Matrices are equal:',bool)
                if not bool:
                    print('Test failed!')
                    sys.exit()
                print()

