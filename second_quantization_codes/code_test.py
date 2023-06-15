import sys
import numpy as np
import hubbard as hs
import hubbard_classes as hc
import focktest as fg
import focktest_classes as fc
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


# Test 1 - gs calculation
if '1' in answer:
    for N in N_values:
        for t in t_values:
            for U in U_values:
                for mu in mu_values:
                    GS,E,GS_bloc_matrix = hs.hubbard(N,t,U,mu)
                    gs_energy,gs_numerical_state,gs_block_matrix = hc.hubbard(N,t,U,mu)
                    print(['N','t','U','mu'])
                    print([N,t,U,mu])
                    print('GS_energies: ',E == gs_energy)
                    if not E == gs_energy:
                        print('Test failed!')
                        sys.exit()
                    print()
                    

# Test 2 - H/S matrix with fock
if '2' in answer:
    for type in types:
        for spin_left in spins_left:
            for spin_right in spins_left:
                for N in N_values:
                    for N_min in N_min_values:
                        for t in t_values:
                            for U in U_values:
                                for mu in mu_values:
                                    result1 = fc.matrix(type,excit_doc,N,N_min,spin_left,spin_right,t,U,mu)
                                    result2 = fg.matrix(type,excit_doc,N,N_min,spin_left,spin_right,t,U,mu)
                                    
                                    print(['type','spin_left','spin_right','N','N_min','t','U','mu'])
                                    print([type,spin_left,spin_right,N,N_min,t,U,mu])
                                    bool = np.allclose(result1, result2, rtol=1e-05, atol=1e-05, equal_nan=False)
                                    print('Matrices are equal:',bool)
                                    if not bool:
                                        print('Test failed!')
                                        sys.exit()
                                    print()
