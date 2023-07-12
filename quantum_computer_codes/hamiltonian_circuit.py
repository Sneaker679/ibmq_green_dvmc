### Packages ################################################
import matplotlib.pyplot as plt
import numpy as np
import sys,os
from qiskit.primitives import Estimator
from qiskit_nature.second_q.mappers import JordanWignerMapper
from qiskit_nature.second_q.operators import FermionicOp
from qiskit import QuantumCircuit,QuantumRegister
from qiskit.quantum_info import Statevector
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


### Fetching parameters.py and importing qcm ################
module_directory = os.path.dirname(__file__)
sys.path.insert(0,module_directory)

from parameters import *

sys.path.insert(0,os.path.join(module_directory,'..','second_quantization_codes'))
from fock_class import Fock

if use_qcm == 'Y':
    import pyqcm
#################### LATTICE SETUP ##########################

"""This loops checks if the number of site can be made into a grid of sites.
This is checked by calculating (in parameters.py) if the number of sites is
a prime number or not. If said number is NOT a prime number (in other words,
if the inputted number of site can be made into a grid of sites) and if we
are not forcing the the code to use a custom lattice (configured in parameters.py,
then with the following lines, we generate the automatic lattice for all 3 
implementations of the green function calculation (qcm,fock and ibmq)."""

if not len(factors) == 2 and force_custom_lattice == 'N':
    print('Using automatic lattice...')
    
    # Find number of rows and columns the automatic lattice would have.
    if len(factors) % 2 == 0:
        num_columns = factors[int(len(factors)/2 - 1)]
        num_rows = factors[int(len(factors)/2)]
    else:
        num_columns = factors[int((len(factors)-1)/2)]
        num_rows = num_columns

    ### Initializing qiskit lattice
    boundary_condition = BoundaryCondition.OPEN
    lattice = SquareLattice(rows=num_rows, cols=num_columns, boundary_condition = boundary_condition)
    #lattice.draw()
    #plt.show()

    # Initializing the Hubbard qiskit Hamiltonian using the qiskit lattice.
    Hamiltonian = FermiHubbardModel(
        lattice.uniform_parameters(
            uniform_interaction = t,
            uniform_onsite_potential = -mu,
        ),
        onsite_interaction = U,
    ).second_q_op()


    ### QCM lattice.
    if use_qcm == 'Y':
        # Initializing cluster model
        CM = pyqcm.cluster_model(N)

        # Generating coordinates of every site for that cluster automatically.
        clus_coordinates = []
        for row in range(num_rows):
            for column in range(num_columns):
                clus_coordinates.append((row,column,0))
        clus_coordinates[:] = (x for x in clus_coordinates)

        # Creating the actual lattice and definning possible hoppings.
        clus = pyqcm.cluster(CM,clus_coordinates)
        model = pyqcm.lattice_model('custom', clus, ((1000,0,0),))
        model.interaction_operator('U')
        model.hopping_operator('t', (1,0,0), -1)  # NN hopping
        model.hopping_operator('t', (0,1,0), -1)  # NN hopping
    

    ### Automatic Fock hopping matrix
    
    # Generating lattice's site coordinates
    lattice_fock = []
    for row in range(num_rows):
        for column in range(num_columns):
            lattice_fock.append(np.matrix([row,column]))

    # Generating hopping matrix using the coordinates defined before
    t_fock = np.zeros((N,N))
    for index1,site1 in enumerate(lattice_fock):
        for index2,site2 in enumerate(lattice_fock):
            
            '''The condition for checking if a hopping is possible is here. We are substracting the coordinates
            and checking if the result is a possible hopping in a grid lattice. If it is, the appropriate term
            is added to the hopping matrix.'''
            hop = np.subtract(site1,site2)
            if (np.array_equal(hop,np.matrix([1,0]))
                or np.array_equal(hop,np.matrix([-1,0]))
                or np.array_equal(hop,np.matrix([0,1]))
                or np.array_equal(hop,np.matrix([0,-1]))):

                t_fock[index1,index2] = t

else:
    '''If the number of sites is a prime number, we use the custom lattice from qiskit in parameters.py to
    generate the associated hamiltonian. The custom lattice we have defined for qcm and fock in parameter.py will be handled by 
    qcm_benchmark.py and fock_benchmark.py.'''
    
    if len(factors) == 2:
        print('Number of sites is a prime number. Using custom lattice...')
    else:
        print('Forcing custom lattice...')
    
    # The lattice here has been initialized in parameters.py
    # We make the hamiltonian out of it.
    Hamiltonian = FermiHubbardModel(
        lattice.uniform_parameters(
            uniform_interaction = t,
            uniform_onsite_potential = -mu,
        ),
        onsite_interaction = U,
    ).second_q_op()
    

#################### CIRCUIT SETUP ####################
'''Qiskit needs a circuit so that it can calculate the desired element of the H and S matrices.
This circuit is either created automatically here, or customized in parameters.py.'''

if force_custom_circuit.upper() == 'N': 
    """We calculate the exact ground state and use that vector to initialize the qubits."""
    print('Using exact diagonalisation for the ground state circuit...')

    eigen_states = np.linalg.eigh(JordanWignerMapper().map(Hamiltonian).to_matrix().real)[1]

    def find_vector_spin(vec):
        tol = 0.00001
        for vec_index,component in enumerate(vec):
            if np.abs(component) > tol:
                vec_spin = Fock(N,vec_index,qiskit_notation='Y').total_spin
                break
        return vec_spin

    gs_index = 0
    while True:
        if spin_gs not in ['+','-']:
            raise Exception("spin_gs must be either '+' or '-'.")
 
        vec = eigen_states[:,gs_index].real.tolist()
        if spin_gs == '+':
            vec_spin = find_vector_spin(vec) 
            if vec_spin >= 0:
                break
            
        if spin_gs == '-':
            vec_spin = find_vector_spin(vec) 
            if vec_spin < 0:
                break

        gs_index += 1

    q = QuantumRegister(2*N)
    qc = QuantumCircuit(q)
    qc.initialize(vec,q)
    if decompose_and_print_circuit == 'Y':
        circuit = qc.decompose(reps=4*N)
        print(circuit)
    else:
        circuit = qc
        del qc
    circuit.draw()

else:
    """We don't define any circuits here because is was already done in parameters.py. The previous if
    overrides the custom circuit for the automatic one."""
    print('Using custom ground state circuit...')
print()

#sv = Statevector.from_instruction(circuit)
#print(sv.to_dict())
#######################################################


########### CALCULATING GROUND STATE ENERGY ###########
qubit_hamiltonian = JordanWignerMapper.mode_based_mapping(Hamiltonian)

job = Estimator().run(circuit,qubit_hamiltonian)
result = job.result()
values = result.values
omega = values[0]
#######################################################

