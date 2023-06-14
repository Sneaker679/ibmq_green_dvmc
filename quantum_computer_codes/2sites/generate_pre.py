import excitation_matrices as em
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

# Building the Hamiltonian of Hubbard's model using qiskit

N = 2 # Number of sites. This entire code only works for 2 sites because of the circuit.
t = -1
U = 4
mu = U/2


    # Initializing latice
    
boundary_condition = BoundaryCondition.OPEN
lattice = LineLattice(num_nodes = N, boundary_condition = boundary_condition)
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


# Circuit.

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

# Asking the user for the calculation

print('Calculation using the quantum computer.')
print('Calculate one or all matrices?  (H+,H-,S+,S- or ALL)  ALL')
print('Generate .npy files? (Y,N) Y')

N = 2
N_min = 1

print('Name of the document containing the excitation types?  excitation.def')
print('Spin of left side?  +')
print('Spin of right side?  +')

for type in ['H+','H-','S+','S-']:
    print(type+':')
    print(em.matrix(type,'excitation.def',N,N_min,'+','+',Hamiltonian,circuit,'Y'))
    print()

