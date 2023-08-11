############## PARAMETERS ###############
# For yes: True
# For no: False

### IBM Credentials and parameters ###
run_on_quantum_computer = False
max_circuit_per_job = 30

channel = "ibm_quantum"
token = "token"
backend_device = "ibm_sherbrooke"
    # List of backends here: https://quantum-computing.ibm.com/services/resources?tab=yours

# Backend Options
custom_qubits = []
optimization_level = 3              # int
resilience_level = 1                # int
max_execution_time = None           # int or None
execution = {
    'shots' : 4000,                 # int
    'init_qubits' : True            # Boolean
    }

# Recover jobs
recover_jobs = False
job_ids = {# Add the ids of the job in the same order they were submitted intially.
    'job0': '',
    'job1': '',
    'job2': '',
    #add more
} 


### Physic parameters ###
N = 2
t = -1
U = 4
mu = 2
spin_green = '+'                                 # Either '+' or '-'.
spin_gs = '+'                                    # Either '+' or '-'.


### Code parameters ###
use_qcm = False

force_custom_lattice = False
force_custom_circuit = False
decompose_and_print_circuit = False
produce_latex_circuit = False

generate_npy = True
generate_matrix = 'ALL'                          # 'H+','H-','S+','S-' or 'ALL'.
excit_document = f'excitation{N}sites.def'


### Noisy simulation ###
noisy_simulation = False

# Noisy Simulation Parameters
"""Some of these options below do not seem to have an effect on the simulation. Qiskit is bugged."""
backend_options = {
    'method': 'automatic',
    'device': 'CPU',
    'precision': 'double',
    'max_job_size': 1,
    'max_shot_size': 8192,
    'max_parallel_threads': 8,
    'max_parallel_experiments': 8,
    'max_parallel_shots': 8
    }
run_options = {
    'shots' : None,                             # int or None
    'seed' : 50                                 # int
    }

"""More info about these options here:
https://qiskit.org/ecosystem/aer/stubs/qiskit_aer.QasmSimulator.html
https://qiskit.org/ecosystem/aer/stubs/qiskit_aer.primitives.Estimator.html
"""


### Graph parameters ###
graph_for_ibmq = True
graph_for_fock = True
graph_for_qcm = True


############### PACKAGES ################
import matplotlib as mpl
from qiskit_ibm_runtime.options import Options
from qiskit import QuantumCircuit,QuantumRegister
from qiskit.compiler import transpile
from qiskit_nature.second_q.hamiltonians.lattices import (
    Lattice,
    BoundaryCondition,
    HyperCubicLattice,
    LatticeDrawStyle,
    LineLattice,
    SquareLattice,
    TriangularLattice,
)
import numpy as np
if use_qcm is True:
    import pyqcm
#########################################


#########################################
####### INPUT CUSTOM LATTICE HERE #######
#######   affects qiskit only     #######
#########################################
boundary_condition = BoundaryCondition.OPEN
lattice = LineLattice(num_nodes = N, boundary_condition = boundary_condition)
#########################################


#########################################
####  INPUT CUSTOM t HOPPING MATRIX  ####
####    affects fock's basis only    ####
#########################################

# A possible hopping is 1
# No possible hopping is 0
hopping_matrix = np.matrix([
    [0,1],
    [1,0]
])

#########################################

factors = []
for i in range(1, N + 1):
       if N % i == 0:
           factors.append(i)

if (use_qcm is True and force_custom_lattice is True) or (use_qcm is True and len(factors) == 2):
    #####################################
    ### INPUT CUSTOM QCM LATTICE HERE ###
    #####################################
    # Input the coordinates of the sites as tuples in "clus_coordinates" (3 dimensional space).
    # You may also need to add "hopping_operator"s depending on your structure.
    
    CM = pyqcm.cluster_model(N)

    clus_coordinates = ((0,0,0),(1,0,0))

    clus = pyqcm.cluster(CM,clus_coordinates)
    model = pyqcm.lattice_model('custom', clus, ((1000,0,0),))
    model.interaction_operator('U')
    model.hopping_operator('t', (1,0,0), -1)  # NN hopping
    model.hopping_operator('t', (0,1,0), -1)  # NN hopping
    #####################################


####################################
#### INPUT CUSTOM CIRCUIT HERE #####
####################################

# Don't modify this following line #
qr = QuantumRegister(2*N,'qr')
circuit = QuantumCircuit(qr)
'''
θ = np.arccos(np.sqrt(8)*0.26136036)
ϕ = np.arccos(2*0.30581423/np.sin(θ))
#we build the spin down:
circuit.ry(2*θ,0)
circuit.ch(target_qubit=2,control_qubit=0,ctrl_state=0)
circuit.cry(2*ϕ,target_qubit=3,control_qubit=0)
circuit.ccx(target_qubit=1,control_qubit1=0,control_qubit2=2,ctrl_state=0)
circuit.cx(target_qubit=0,control_qubit=3)

#we build the spin up:
circuit.h(4)
circuit.h(5)
circuit.ccx(target_qubit=6,control_qubit1=5,control_qubit2=4)
circuit.ccx(target_qubit=7,control_qubit1=5,control_qubit2=4,ctrl_state=0)
circuit.cx(target_qubit=5,control_qubit=6)
circuit.cx(target_qubit=4,control_qubit=6)

# we apply the swaps:
circuit.cswap(control_qubit=4, target_qubit1=3, target_qubit2=0)
circuit.cswap(control_qubit=5, target_qubit1=3, target_qubit2=1)
circuit.cswap(control_qubit=5, target_qubit1=2, target_qubit2=0)
circuit.cswap(control_qubit=6, target_qubit1=3, target_qubit2=2)
circuit.cswap(control_qubit=6, target_qubit1=1, target_qubit2=0)

# we apply the phase:
circuit.cz(control_qubit=4,target_qubit=3)
circuit.cz(control_qubit=4,target_qubit=2)
circuit.cz(control_qubit=4,target_qubit=1)

circuit.cz(control_qubit=5,target_qubit=2)
circuit.cz(control_qubit=5,target_qubit=3)

circuit.cz(control_qubit=6,target_qubit=3)


circuit.draw('mpl')
import matplotlib.pyplot as plt
plt.show()


circuit = transpile(circuit, initial_layout=[qr[0],qr[4],qr[1],qr[5],qr[2],qr[6],qr[3],qr[7]])
'''
####################################


############ FILE PATHS ############
import sys,os

# Junk output directory
output_directory = os.path.join(os.getcwd(),'output')
if not os.path.exists(output_directory):
    os.mkdir(output_directory)

# Pdf output directory
pdf_output_directory = os.getcwd()
####################################

### DON'T MODIFY
quantum_computer_options = Options(
    optimization_level = optimization_level,
    resilience_level = resilience_level,
    max_execution_time = max_execution_time,
    execution = execution
    )

aer_estimator_options = {
    'backend_options' : backend_options,
    'run_options' : run_options
    }
