#Packages----------------------------------------------------

import numpy as np
import matplotlib.pyplot as plt

#import functions as f
#import hubbard_simplified as hs

from qiskit.quantum_info import Pauli,Operator
from qiskit.primitives import Estimator as pEstimator
from qiskit_nature.second_q.mappers import JordanWignerMapper
from qiskit_nature.second_q.operators import FermionicOp
from qiskit import QuantumCircuit, QuantumRegister
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
import sys
import copy
import time

#------------------------------------------------------------

# Print options
np.set_printoptions(linewidth= 1000)


# Creation, destruction and count(n) operators using qiskit

"""
N is the number of sites in total.
site is the site for which we want to create the operator. The first site is #0.
spin is the spin to be used for the creation of the operator.
"""
def create(N,site,spin): 
    site += 1
    if not spin == '-' and not spin == '+':
        raise Exception('Input must be "+" or "-".')
    if spin == '+':
        spin = 1
    if spin == '-':
        spin = 0

    """The following lines create an FermionicOp, which is a class in qiskit.
    The equation in the parenthesis is simply to accomodate qiskit's notation for
    the fock space.

    | 1,up 1,down 2,up 2,down >
    """
    operator = FermionicOp(
        {
            '+_' + str(2*site-spin-1): 1.0,
        },
        num_spin_orbitals=2*N,
    )
    return operator

"""Since the annihilation, or as I call it here, the 'destroy' operator, is the
equivalent of create**dagger, we only need to do the transpose() and conjugate()
of 'creation' to obtain 'destroy'."""
def destroy(N,site,spin):
    return create(N,site,spin).transpose().conjugate()

"""The check operator is simply the name of the operator 'n', which checks if an
electron is there or not. As per qiskit's FermionicOp class, multiplication of 
operators must be written using @."""
def check(N,site,spin):
    return create(N,site,spin) @ destroy(N,site,spin)


# Calculation of any excited state

"""This function reads a excitation.def file and converts it into a easily manipulated object.
In this case, said object is a list. Each items of this list are the individual lines. 
These items are also themselves lists, each containing in this same order: t ri ra rb."""
def excitdef_reader(document_name):

    file = open(document_name).read()
    lines_doc = file.split('\n')[5:-1]
    for line_number in range(len(lines_doc)):
        lines_doc[line_number] = lines_doc[line_number].split() 
        lines_doc[line_number] = [eval(number) for number in lines_doc[line_number]]
    
    return lines_doc

"""This function calculates the product of the operators associated with an exited state.
It is important to note that this doesn't calculate the product of those operators with the GS.
This task is for the quantum computer."""
def ex_operators(side,type,i,m,spin,lines_doc):
    
    """The following block merely ensure proper values are inserted in the function."""
    valid_side = ['left','right']
    valid_type = ['H+','H-','S+','S-']
    valid_spin = ['+','-']
    if side not in valid_side:
        raise Exception('Side has to be any one of these',valid_site,'.')
    if type not in valid_type:
        raise Exception('Type has to be any one of these:',valid_type,'.')
    if i < 0:
        raise Exception('Site numbers start at 0.')
    if m < 0:
        raise Exception('m starts at 0.')
    if spin not in valid_spin:
        raise Exception('Spin has to be any one of these:',valid_spin,'.')

    """This code block specifies the 'c' operator to be used in, for example, |e> = cnn|GS>.
    In the last example, the c should be a c_dagger."""
    if type[1] == '+':
        c_operator = create(N,i,spin)
    if type[1] == '-':
        c_operator = destroy(N,i,spin)
    
    """For the spin defined, we also define its opposite."""
    if spin == '+':
        spin_op = '-'
    if spin == '-':
        spin_op = '+'
    
    """For each definition of an excited state, we define the operation to be performed."""
    if m == 0:
        ex_ops = c_operator
    elif m == 1:
        ex_ops = c_operator @ check(N,i,spin_op)
    else:
        """Here, we fetch the values of ra and rb from the lines_doc we inputed in the function."""
        lines = [values for values in lines_doc if values[1] == i]
        ra = lines[m][2]
        rb = lines[m][3]

        ex_ops = c_operator @ check(N,ra,spin_op) @ check(N,rb,spin)

    return ex_ops


"""This is the function that calculates the observable to be used by the quantum
computer. """

"""
type is the type of the matrix we wish to calculate -> H+, H-, S+ or S-.
excit_document is the document containing all the excited states we want to calculate.
N is the number of sites in total.
"""

def Observable(type,excit_document,N,i,j,m,n,spin_left,spin_right,hamiltonian):
    
    """The excitation.def file in its list form is imported."""
    lines_doc = excitdef_reader(excit_document)

    """We use the previously defined functions to calculate, for example, |e_right> and <e_left|.
    Notice how ex_op_left has .transpose() and .conjugate at its end, because we want the bra, not the ket."""
    ex_op_left = ex_operators('left',type,i,m,spin_left,lines_doc).transpose().conjugate()
    ex_op_right = ex_operators('right',type,j,n,spin_right,lines_doc)
    
    """The result is the final observable to be used by the quantum computer."""
    if type[0] == 'H':
        return ex_op_left @ hamiltonian @ ex_op_right
    if type[0] == 'S':
        return ex_op_left @ ex_op_right



# Calculation of the entire matrix using the quantum computer
"""
type is the type of matrix we wish to calculate -> H+, H-, S+ or S-.
excit_document is the document containing all the excited states we want to calculate.
N is the number of sites in total.
N_min is the number of neighbors the site with the less neighbors has.
hamiltonian is the hamiltonian to use for the calculation.
save is for if we wish to save the matrices we calculate as a .npy file.
"""

def matrix(type,excit_document,N,N_min,spin_left,spin_right,hamiltonian,q_circuit,save='N'):
    
    # Timer start
    start = time.time()
    
    # This is the total possible number of excitation for each site, as explained by the paper.
    N_exc = 2 + N_min*(N_min+1)

    # Initializing matrix
    matrix_size = N * N_exc
    excitation_matrix = np.zeros((matrix_size,matrix_size))

    # Caluculation of half the elements 
    for i in range(N):
        for j in range(N):
            for m in range(N_exc):
                for n in range(N_exc): 
                    """The 2 following lines represents the distribution of the calculated values in the
                    matrix. For each i, there are (for 2 sites) 4 values of m. Thus, the first 4 columns
                    would be for i = 1 and for 0 <= m <= 3, then the next 4 for i = 2 and for 0 <= m <= 3.
                    The same principles for the rows, but i becomes j and m becomes n."""
                    column_num = N * m + i
                    row_num = N * n + j
                    if column_num - row_num >= 0: # This is to calculate only half of the matrix.
                        observable = Observable(type,excit_document,N,i,j,m,n,spin_left,spin_right,hamiltonian)
                        qubit_hamiltonian = JordanWignerMapper.mode_based_mapping(observable)

                        job = pEstimator().run(q_circuit,qubit_hamiltonian)
                        result = job.result()
                        values = result.values 

                        excitation_matrix[row_num,column_num] = values[0]

    # Symmetrizing the matrix
    excitation_matrix = np.tril(excitation_matrix.T,-1) + excitation_matrix

    end = time.time()
    print('Time:',end - start,'seconds.')

    if save.upper() == 'Y':
        if type[1] == '+':
            identifier = '_AC'
        if type[1] == '-':
            identifier = '_CA'
        np.save(type[0]+identifier,excitation_matrix)

    return excitation_matrix


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
print(np.linalg.eigh(Hamiltonian.to_matrix().toarray())[1][:,0].real.tolist())
vec = np.linalg.eigh(Hamiltonian.to_matrix().toarray())[1][:,0].real.tolist()
#desired_vector = [-0.270598050,0.653281482,-0.653281482,-0.270598050]

q = QuantumRegister(4)

circuit = QuantumCircuit(q)

circuit.initialize(vec, [q[0],q[1],q[2],q[3]])
circuit.draw()
#print(circuit)
#print()
#plt.show()


# Asking the user for the calculation
'''
print('Calculation using the quantum computer.')
generate_matrix = input('Calculate one or all matrices?  (H+,H-,S+,S- or ALL)  ')
generate_npy = input('Generate .npy files? (Y,N) ')

N = 2
N_min = 1

excit_document= input('Name of the document containing the excitation types?  ')
spin_left = input('Spin of left side?  ')
spin_right = input('Spin of right side?  ')

if generate_matrix.upper() == 'ALL':
    for type in ['H+','H-','S+','S-']:
        print(type+':')
        print(matrix(type,excit_document,N,N_min,spin_left,spin_right,Hamiltonian,circuit,generate_npy))
        print()
else:
    print(generate_matrix+':')
    print(matrix(generate_matrix,excit_document,N,N_min,spin_left,spin_right,Hamiltonian,circuit,generate_npy))
'''
