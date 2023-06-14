#Packages----------------------------------------------------

import numpy as np
import matplotlib.pyplot as plt

import functions as f
import hubbard_simplified as hs

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

"""This is the function that calculates the observable to be used by the quantum
computer. """

def Observable(type,N,i,j,m,n,spin_left,spin_right,hamiltonian):

    """The following lines merely ensure that no wrong inputs are given to the function."""
    if type != 'H-' and type != 'H+' and type != 'S-' and type != 'S+':
        raise Exception('The "type" parameter has to be "H-","H+","S-" or "S+".')
    if m < 0 or n < 0:
        raise Exception('Parameters "m" and "n" have to be positive integers.')
    if (spin_right != '+' and spin_right != '-') or (spin_left != '-' and spin_left != '+'):
        raise Exception('Parameters "spin_left" and "spin_right" must be equal to "-" or "+" (both being strings).')

    """spin_op means opposite spin."""
    if spin_left == '+':
        spin_left_op = '-'
    elif spin_left == '-':
        spin_left_op = '+'
    if spin_right == '+':
        spin_right_op = '-'
    elif spin_right == '-':
        spin_right_op = '+'
    
    """A file detailing the excitations we wish to calculate needs to exist in the
    same folder as this code.

    The following lines transform the .def file into an integer list."""
    file = open('excitation.def').read()
    lines = file.split('\n')[5:-1]
    for line_number in range(len(lines)):
        lines[line_number] = lines[line_number].split() 
        lines[line_number] = [eval(number) for number in lines[line_number]]

    """The following lines calculate the left side of the observable. Depending
    on if we are calculating H+, H-, S+ or S-, the definition of ex_state_left changes.

    ex_state stands for excited state.
    """
    if m == 0:
        if type == 'H+' or type == 'S+':
            ex_state_left = create(N,i,spin_left)
        if type == 'H-' or type == 'S-':
            ex_state_left = destroy(N,i,spin_left)
    elif m == 1:
        if type == 'H+' or type == 'S+':
            ex_state_left = create(N,i,spin_left) @ check(N,i,spin_left_op)
        if type == 'H-' or type == 'S-':
            ex_state_left = destroy(N,i,spin_left) @ check(N,i,spin_left_op)
    else:
        # This line filters the 'lines' list to only include the left's side necessary lines.
        lines_left = [ex_state for ex_state in lines if ex_state[1] == i]
        
        ra = lines_left[m][2]
        rb = lines_left[m][3]
         
        if type == 'H+' or type == 'S+':
            ex_state_left = create(N,i,spin_left) @ check(N,ra,spin_left_op) @ check(N,rb,spin_left)
        if type == 'H-' or type == 'S-':
            ex_state_left = destroy(N,i,spin_left) @ check(N,rb,spin_left) @ check(N,ra,spin_left_op)
    
    ex_state_left = ex_state_left.transpose().conjugate()
    

    """The following lines calculate the right side of the observable."""
    if n == 0:
        if type == 'H+' or type == 'S+':
            ex_state_right = create(N,j,spin_right)
        if type == 'H-' or type == 'S-':
            ex_state_right = destroy(N,j,spin_right)
    elif n == 1:
        if type == 'H+' or type == 'S+':
            ex_state_right = create(N,j,spin_right) @ check(N,j,spin_right_op)
        if type == 'H-' or type == 'S-':
            ex_state_right = destroy(N,j,spin_right) @ check(N,j,spin_right_op)
    else:
        lines_right = [ex_state for ex_state in lines if ex_state[1] == j]
       
        ra = lines_right[n][2]
        rb = lines_right[n][3]
        
        if type == 'H+' or type == 'S+':
            ex_state_right = create(N,j,spin_right) @ check(N,ra,spin_right_op) @ check(N,rb,spin_right)
        if type == 'H-' or type == 'S-':
            ex_state_right = destroy(N,j,spin_right) @ check(N,rb,spin_right) @ check(N,ra,spin_right_op)

    if type == 'H+' or type == 'H-':
        return ex_state_left @ hamiltonian @ ex_state_right
    if type == 'S+' or type == 'S-':
        return ex_state_left @ ex_state_right


# Calculation of the entire matrix using the quantum computer
"""
type is the type of matrix we wish to calculate -> H+, H-, S+ or S-.
N is the number of sites in total.
N_min is the number of neighbors the site with the less neighbors has.
spin_left is the spin input for the left side.
spin_right is the spin input for the right side.
hamiltonian is the hamiltonian to use for the calculation.
"""

def matrix(type,N,N_min,spin_left,spin_right,hamiltonian):
    
    # Timer start
    start = time.time()
    
    # This is the total possible number of excitation for each site, as stated by the paper.
    N_exc = 2 + N_min*(N_min+1)

    # Initializing matrix
    matrix_size = N * N_exc
    excitation_matrix = np.zeros((matrix_size,matrix_size))
    
    # Initializing estimator
    estimator = pEstimator()

    # Caluculation of half the elements 
    for i in range(N):
        for j in range(N):
            for m in range(N_exc):
                for n in range(N_exc): 
                    """The 2 following lines represents the distribution of the calculated values in the
                    matrix. For each i, there are (for 2 sites) 4 values of m. Thus, the first 4 columns
                    would be for i = 1 and for 0 <= m <= 3, then the next 4 for i = 2 and for 0 <= m <= 3.
                    The same principles for the rows, but i becomes j and m becomes n."""
                    column_num = N_exc * i + m
                    row_num = N_exc * j + n
                    if column_num - row_num >= 0: # This is to calculate only half of the matrix.
                        observable = Observable(type,N,i,j,m,n,spin_left,spin_right,hamiltonian)
                        qubit_hamiltonian = JordanWignerMapper.mode_based_mapping(observable)

                        job = estimator.run(circuit,qubit_hamiltonian)
                        result = job.result()
                        values = result.values 

                        excitation_matrix[row_num,column_num] = values[0]

    # Symmetrizing the matrix
    excitation_matrix = np.tril(excitation_matrix.T,-1) + excitation_matrix
    np.save("H_CA",excitation_matrix)

    end = time.time()
    print('Time:',end - start,'seconds.')

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
print(circuit)
print()
#plt.show()


# Asking the user for the calculation

print('Calculation using the quantum computer.')
type = input('Which matrix? (H+, H-, S+ or S-.)  ')
N_exc = 1
spin_left = input('Spin of left side?  ')
spin_right = input('Spin of right side?  ')

print(matrix(type,2,1,spin_left,spin_right,Hamiltonian))





