### Packages ################################################
from math import isclose
import matplotlib.pyplot as plt
import numpy as np
import sys,os
from qiskit.primitives import Estimator
from qiskit_nature.second_q.mappers import JordanWignerMapper
from qiskit_nature.second_q.operators import FermionicOp
from qiskit import QuantumCircuit,QuantumRegister
from qiskit.quantum_info import Statevector
from qiskit_nature.second_q.hamiltonians import FermiHubbardModel
from qiskit_nature.second_q.hamiltonians.lattices import Lattice

### Fetching parameters.py and importing qcm ################
module_directory = os.path.dirname(__file__)
sys.path.insert(0,module_directory)

from parameters import *

sys.path.insert(0,os.path.join(module_directory,'..','second_quantization_codes'))
from fock_class import Fock

if use_qcm is True:
    import pyqcm

#################### LATTICE SETUP ##########################

"""This loops checks if the number of site can be made into a grid of sites.
This is checked by calculating (in parameters.py) if the number of sites is
a prime number or not. If said number is NOT a prime number (in other words,
if the inputted number of site can be made into a grid of sites) and if we
are not forcing the the code to use a custom lattice (configured in parameters.py,
then with the following lines, we generate the automatic lattice for all 3 
implementations of the green function calculation (qcm,fock and ibmq)."""

if not len(factors) == 2 and force_custom_lattice is False:
    print('Using automatic lattice...')
    
    # Find number of rows and columns the automatic lattice would have.
    if len(factors) % 2 == 0:
        num_columns = factors[int(len(factors)/2 - 1)]
        num_rows = factors[int(len(factors)/2)]
    else:
        num_columns = factors[int((len(factors)-1)/2)]
        num_rows = num_columns

    ### QCM lattice.
    if use_qcm is True:
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
    

    ### Automatic hopping matrix
    # Generating lattice's site coordinates
    lattice_fock = []
    for row in range(num_rows):
        for column in range(num_columns):
            lattice_fock.append(np.matrix([row,column]))

    # Generating hopping matrix using the coordinates defined before
    hopping_matrix = np.zeros((N,N))
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

                hopping_matrix[index1,index2] = 1
    
    ### Automatic Qiskit lattice
    #lattice = Lattice.from_adjacency_matrix(hopping_matrix)
    boundary_condition = BoundaryCondition.OPEN
    lattice = SquareLattice(rows=num_columns, cols=num_rows, boundary_condition = boundary_condition)


else:
    '''If the number of sites is a prime number, we use the custom lattice from qiskit in parameters.py to
    generate the associated hamiltonian. The custom lattice we have defined for qcm and fock in parameter.py will be handled by 
    qcm_benchmark.py and fock_benchmark.py.'''
    
    if len(factors) == 2:
        print('Number of sites is a prime number. Using custom lattice...')
    else:
        print('Forcing custom lattice...')


# Initializing the Hubbard qiskit Hamiltonian using the qiskit lattice.
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

if force_custom_circuit is False: 

    # Choosing circuit if N == 2. These circuits are hardcoded.
    continue_with_diag = False
    if N == 2:
        q = QuantumRegister(2*N)
        circuit = QuantumCircuit(q)
        if mu >= 5:
            circuit.x(0)
            circuit.x(1)
            circuit.x(2)
            circuit.x(3)
        elif mu < 5 and mu >= 3.8286 and spin_gs == '+':
            circuit.x(0)
            circuit.h(0)
            circuit.cx(0,2)
            circuit.x(2)
            circuit.x(1)
            circuit.x(3)
        elif mu < 5 and mu >= 3.8286 and spin_gs == '-':
            circuit.x(3)
            circuit.h(3)
            circuit.cx(3,1)
            circuit.x(1)
            circuit.x(0)
            circuit.x(2)
        elif mu < 3.8286 and mu > 0.1714:
            eigen = np.linalg.eigh(JordanWignerMapper().map(Hamiltonian).to_matrix().real)
            gs_vec = eigen[1][:,0].tolist()
            angle = np.arccos(-gs_vec[3]*2**(1/2))
            theta = 2*angle

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
        elif mu <= 0.1714 and mu > -1 and spin_gs == '+':
            circuit.h(3)
            circuit.cx(3,1)
            circuit.x(3)
        elif mu <= 0.1714 and mu > -1 and spin_gs == '-':
            circuit.h(0)
            circuit.cx(0,2)
            circuit.x(0)
        elif mu <= -1:
            ...# State is |0000>
        else:
            print('No hardcoded circuits available for that configuration.\nProceeding with exact diagonalisation.')
            continue_with_diag = True
    

    if not N == 2 or continue_with_diag is True:
        """We calculate the exact ground state and use that vector to initialize the qubits."""
        print('Using exact diagonalisation for the ground state circuit...')

        eigen = np.linalg.eigh(JordanWignerMapper().map(Hamiltonian).to_matrix().real)
        eigen_states = eigen[1] 
        eigen_energies = eigen[0]
        
        """In order to prioritize the ground state with the spin closest to 0, we need the following code."""
        def find_vector_spin(vec):
            """
            vec: Vector as a list.

            Returns: Vector state's total spin as an integer
            """

            tol = 0.00001
            for vec_index,component in enumerate(vec):
                if np.abs(component) > tol:
                    vec_spin = Fock(N,vec_index,qiskit_notation=True).total_spin
                    break
            return vec_spin
        
        # Fetching the number of states that have the same ground state energy
        gs_eigen_energies = []
        gs_energy = eigen_energies[0]
        tol = 1e-5
        for energy in eigen_energies:
            if isclose(energy, gs_energy, abs_tol=tol):
                gs_eigen_energies.append(energy)
            else:
                break

        # Fetching the states that have the same ground state energy.
        eigen_states_list = []
        for index in range(len(gs_eigen_energies)):
            eigen_states_list.append(eigen_states[:,index])

        # Looping trough the states to find the one witht the spin closest to 0.
        for vec in eigen_states_list:
            if spin_gs not in ['+','-','0']:
                raise Exception("spin_gs must be either '+', '-' or '0'.")
            
            vec_spin = find_vector_spin(vec)
            if ((spin_gs == '+' and vec_spin >= 0)
            or (spin_gs == '-' and vec_spin <= 0)):
                gs_vec = vec 

        # Initializing quantum circuit
        q = QuantumRegister(2*N)
        qc = QuantumCircuit(q)
        qc.initialize(gs_vec,q)
        if decompose_and_print_circuit is True:
            circuit = qc.decompose(reps=4*N)
            print(circuit)
        else:
            circuit = qc
            del qc
        circuit.draw()


if force_custom_circuit is True:
    """We don't define any circuits here because is was already done in parameters.py. The previous if
    overrides the custom circuit for the automatic one."""
    print('Using custom ground state circuit...')
print()

#######################################################


########### CALCULATING GROUND STATE ENERGY ###########
qubit_hamiltonian = JordanWignerMapper.mode_based_mapping(Hamiltonian)

job = Estimator().run(circuit,qubit_hamiltonian)
result = job.result()
values = result.values
omega = values[0]
#######################################################

