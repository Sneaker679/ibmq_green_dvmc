############## PARAMETERS ###############
# For yes: 'Y'
# For no: 'N'

# Physic parameters
N = 2                                            # Number of sites.
t = -1
U = 4
mu = 2
spin_green = '+'                                 # Either '+' or '-'.
spin_gs = '+'                                    # Either '+' or '-'.

# Code parameters
use_qcm = 'N'
force_custom_lattice = 'N'
force_custom_circuit = 'N'
decompose_and_print_circuit = 'N'
generate_npy = 'Y' 
generate_matrix = 'ALL'                         # 'H+','H-','S+','S-' or 'ALL'.
excit_document = f'excitation{N}sites.def'

# Simulation parameters
noisy_simulation = 'N'
"""These options below do not seem to have an effect on the simulation. Qiskit is bugged."""
estimator_options = {
        'method': 'automatic',
        #'executor':, 
        'device': 'CPU',
        'precision': 'double',
        'max_job_size': 8,
        'max_shot_size': 1,
        'max_parallel_threads': 8,
        'max_parallel_experiments': 8,
        'max_parallel_shots': 8
        }

"""More info about these options here:
https://qiskit.org/ecosystem/aer/stubs/qiskit_aer.QasmSimulator.html
"""

############### PACKAGES ################
import matplotlib as mpl
from qiskit import QuantumCircuit,QuantumRegister
from qiskit_nature.second_q.hamiltonians.lattices import (
    BoundaryCondition,
    HyperCubicLattice,
    Lattice,
    LatticeDrawStyle,
    LineLattice,
    SquareLattice,
    TriangularLattice,
)
import numpy as np
if use_qcm == 'Y':
    import pyqcm
#########################################


#########################################
####### INPUT CUSTOM LATTICE HERE #######
#########################################
boundary_condition = BoundaryCondition.OPEN
lattice = LineLattice(num_nodes = N, boundary_condition = boundary_condition)
#lattice.draw()
#plt.show()
#########################################

factors = []
for i in range(1, N + 1):
       if N % i == 0:
           factors.append(i)

if (use_qcm == 'Y' and force_custom_lattice == 'Y') or (use_qcm == 'Y' and len(factors) == 2):
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


#########################################
##### INPUT CUSTOM t HOPPING MATRIX #####
#########################################
t_fock = np.matrix([
    [0,t],
    [t,0]
])
#########################################



####################################
#### INPUT CUSTOM CIRCUIT HERE #####
####################################

# Don't modify this following line #
circuit = QuantumCircuit(2*N)
'''
#Circuit for 2 sites, mu=2

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
'''
'''
#Circuit for 2 sites, mu=-0.7

circuit.h(3)
circuit.cx(3,1)
circuit.x(3)
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
