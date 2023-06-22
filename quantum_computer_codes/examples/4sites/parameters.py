############## PARAMETERS ###############

use_qcm = 'Y' # 'Y' or 'N'.
force_custom_lattice = 'N'
force_custom_circuit = 'N'

N = 4 # Number of site.

t = -1
U = 4
mu = U/2

spin_left = '+' # Either '+' or '-'.
spin_right = '+'

excit_document = f'excitation{N}sites.def'
generate_npy = 'Y' 
generate_matrix = 'ALL' # Should be left on 'ALL'.


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


####################################
#### INPUT CUSTOM CIRCUIT HERE #####
####################################

# Don't modify this following line #
circuit = QuantumCircuit(2*N)

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
circuit.draw()
####################################
