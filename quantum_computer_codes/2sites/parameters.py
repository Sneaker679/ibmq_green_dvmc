# Parameters

use_qcm = 'Y' # 'Y' or 'N'.
force_custom_lattice = 'N'
force_custom_circuit = 'N'

N = 2 # Number of sites.
N_min = 1 # Number of sites the site with the least neighbors can interact with.

t = -1
U = 4
mu = U/2

spin_left = '+' # Either '+' or '-'.
spin_right = '+'

excit_document = f'excitation{N}sites.def'
generate_npy = 'Y' 
generate_matrix = 'ALL' # Should be left on 'ALL'.


#################### LATTICE SETUP ####################
import sys
import numpy as np
import matplotlib.pyplot as plt
from qiskit.quantum_info import Pauli,Operator
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
if use_qcm == 'Y':
    import pyqcm

# Find number of rows and columns of lattice
factors = []
for i in range(1, N + 1):
       if N % i == 0:
           factors.append(i)

if not len(factors) == 2 and force_custom_lattice == 'N':
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

else:
    if len(factors) == 2:
        print('Number of sites is a prime number.')
    print('Using custom lattice...')
    print()

    ### Custom qiskit lattice ### 
    # Remove the (#)s to show the lattice using matplotlib.
    boundary_condition = BoundaryCondition.OPEN
    lattice = LineLattice(num_nodes = N, boundary_condition = boundary_condition)
    #lattice.draw()
    #plt.show()
    ############################

    # Don't change. Initializing the Hubbard Hamiltonian using the lattice.
    Hamiltonian = FermiHubbardModel(
        lattice.uniform_parameters(
            uniform_interaction = t,
            uniform_onsite_potential = -mu,
        ),
        onsite_interaction = U,
    ).second_q_op()
    
    if use_qcm == 'Y':
        
        ### Custom QCM lattice.###
        #Input the coordinates of the sites as tuples in "clus_coordinates" (3 dimensional space).
        # You may also need to add "hopping_operator"s depending on your structure.
        CM = pyqcm.cluster_model(N)

        clus_coordinates = ((0,0,0),(1,0,0))

        clus = pyqcm.cluster(CM,clus_coordinates)
        model = pyqcm.lattice_model('custom', clus, ((1000,0,0),))
        model.interaction_operator('U')
        model.hopping_operator('t', (1,0,0), -1)  # NN hopping
        model.hopping_operator('t', (0,1,0), -1)  # NN hopping
        #########################


#################### CIRCUIT SETUP ####################

if force_custom_circuit.upper() == 'N':
    vec = np.linalg.eigh(JordanWignerMapper().map(Hamiltonian).to_matrix())[1][:,0].real.tolist()
    
    q = QuantumRegister(2*N)
    circuit = QuantumCircuit(q)
    circuit.initialize(vec,q)
    circuit.draw()

else:
    circuit = QuantumCircuit(2*N)
    
    ### INPUT CUSTOM CIRCUIT HERE ###
    theta = 2*1.178097245
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
    #################################

    circuit.draw('mpl')
    #print(circuit)
    #print()
    #plt.show()
