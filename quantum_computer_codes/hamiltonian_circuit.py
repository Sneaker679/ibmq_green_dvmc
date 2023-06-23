#################### LATTICE SETUP ####################
import sys,os
import matplotlib.pyplot as plt
import numpy as np
if len(sys.argv) == 2:
    sys.path.insert(0,'./examples/'+sys.argv[1]+'sites')

from parameters import *

from qiskit.quantum_info import Pauli,Operator
from qiskit.primitives import Estimator as pEstimator
from qiskit_nature.second_q.mappers import JordanWignerMapper
from qiskit_nature.second_q.operators import FermionicOp
from qiskit_nature.second_q.problems import EigenstateResult,LatticeModelProblem
from qiskit_nature.second_q.algorithms import GroundStateEigensolver
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
if use_qcm == 'Y':
    import pyqcm

# Find number of rows and columns of lattice
if not len(factors) == 2 and force_custom_lattice == 'N':
    print('Using automatic lattice...')

    if len(factors) % 2 == 0:
        num_columns = factors[int(len(factors)/2 - 1)]
        num_rows = factors[int(len(factors)/2)]
    else:
        num_columns = factors[int((len(factors)-1)/2)]
        num_rows = num_columns

    # Initializing latice
    boundary_condition = BoundaryCondition.OPEN
    lattice = SquareLattice(rows=num_rows, cols=num_columns, boundary_condition = boundary_condition)
    #lattice.draw()
    #plt.show()

    # Initializing the Hubbard Hamiltonian using the qiskit lattice.
    Hamiltonian = FermiHubbardModel(
        lattice.uniform_parameters(
            uniform_interaction = t,
            uniform_onsite_potential = -mu,
        ),
        onsite_interaction = U,
    ).second_q_op()

    if use_qcm == 'Y':
        
        # QCM lattice.
        CM = pyqcm.cluster_model(N)

        clus_coordinates = []
        for row in range(num_rows):
            for column in range(num_columns):
                clus_coordinates.append((row,column,0))
        clus_coordinates[:] = (x for x in clus_coordinates)

        clus = pyqcm.cluster(CM,clus_coordinates)
        model = pyqcm.lattice_model('custom', clus, ((1000,0,0),))
        model.interaction_operator('U')
        model.hopping_operator('t', (1,0,0), -1)  # NN hopping
        model.hopping_operator('t', (0,1,0), -1)  # NN hopping
    
    # Automatic Fock hopping matrix
    lattice_fock = []
    for row in range(num_rows):
        for column in range(num_columns):
            lattice_fock.append(np.matrix([row,column]))
    t_fock = np.zeros((N,N))
    for index1,site1 in enumerate(lattice_fock):
        for index2,site2 in enumerate(lattice_fock):
            hop = np.subtract(site1,site2)
            if (np.array_equal(hop,np.matrix([1,0]))
                or np.array_equal(hop,np.matrix([-1,0]))
                or np.array_equal(hop,np.matrix([0,1]))
                or np.array_equal(hop,np.matrix([0,-1]))):

                t_fock[index1,index2] = t

else:
    if len(factors) == 2:
        print('Number of sites is a prime number.')
    print('Using custom lattice...')
    
    # The lattice here has been initialized in parameters.py
    Hamiltonian = FermiHubbardModel(
        lattice.uniform_parameters(
            uniform_interaction = t,
            uniform_onsite_potential = -mu,
        ),
        onsite_interaction = U,
    ).second_q_op()
    

#################### CIRCUIT SETUP ####################

if force_custom_circuit.upper() == 'N': 
    print('Using exact diagonalisation state...')

    vec = np.linalg.eigh(JordanWignerMapper().map(Hamiltonian).to_matrix())[1][:,0].real.tolist()
    
    q = QuantumRegister(2*N)
    circuit = QuantumCircuit(q)
    circuit.initialize(vec,q)
    circuit.draw()
    #print(circuit)

else:
    print('Using custom circuit...')
print()
#######################################################


########### CALCULATING GROUND STATE ENERGY ###########
qubit_hamiltonian = JordanWignerMapper.mode_based_mapping(Hamiltonian)

job = pEstimator().run(circuit,qubit_hamiltonian)
result = job.result()
values = result.values
omega = values[0]
#######################################################

