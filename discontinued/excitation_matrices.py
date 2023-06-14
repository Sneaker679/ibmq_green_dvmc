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
equivalent to create**dagger, we only need to do the transpose() and conjugate()
of 'creation' to obtain 'destroy'."""
def destroy(N,site,spin):
    return create(N,site,spin).transpose().conjugate()

"""The check operator is simply the name of the operator 'n', which checks if an
electron is there or not. As per qiskit's FermionicOp class, multiplication of 
operators must be written using @."""
def check(N,site,spin):
    return create(N,site,spin) @ destroy(N,site,spin)


# Calculation of any excited state

"""This is the function that calculates the observable(!!) to be used by the quantum
computer. Calling this function 'H' is, in reality, not correct and an abuse of
notation. The actual calculation of the element H_ijmn is done later.

Using this code's function names, H_ijmn is thus = <GS|H_(e/h)|GS>.

Furthermore, the comments provided for H_e also apply to H_h, S_e and S_h.
The 'e' stands for when we add an electron and 'h' stands for when we remove one (hole).

Granted, this code could use simplifications.
"""


def H_e(N,i,j,m,n,spin_left,spin_right,hamiltonian):
    if m < 0 or n < 0:
        raise('Parameter "m" has to be a positive integer.')
    if (spin_right != '+' and spin_right != '-') or (spin_left != '-' and spin_left != '+'):
        raise('Parameters "spin" must be equal to "-" or "+" (both being strings).')
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

    """The following lines calculate the left side of the observable."""
    if m == 0:
        ex_state_left = create(N,i,spin_left)
    elif m == 1:
        ex_state_left = create(N,i,spin_left) @ check(N,i,spin_left_op)
    else:
        # This line filters the 'lines' list to only include the left's side necessary lines.
        lines_left = [ex_state for ex_state in lines if ex_state[1] == i]
        
        ra = lines_left[m][2]
        rb = lines_left[m][3]
        
        ex_state_left = create(N,i,spin_left) @ check(N,ra,spin_left_op) @ check(N,rb,spin_left)
    ex_state_left = ex_state_left.transpose().conjugate()
        
    """The following lines calculate the right side of the observable."""
    if n == 0:
        ex_state_right = create(N,j,spin_right)
    elif n == 1:
        ex_state_right = create(N,j,spin_right) @ check(N,j,spin_right_op)
    else:
        lines_right = [ex_state for ex_state in lines if ex_state[1] == j]
       
        ra = lines_right[n][2]
        rb = lines_right[n][3]
        
        ex_state_right = create(N,j,spin_right) @ check(N,ra,spin_right_op) @ check(N,rb,spin_right)

    return ex_state_left @ hamiltonian @ ex_state_right

def H_h(N,i,j,m,n,spin_left,spin_right,hamiltonian):
    if m < 0 or n < 0:
        raise('Parameter "m" has to be a positive integer.')
    if (spin_right != '+' and spin_right != '-') or (spin_left != '-' and spin_left != '+'):
        raise('Parameters "spin" must be equal to "-" or "+" (both being strings).')
    if spin_left == '+':
        spin_left_op = '-'
    elif spin_left == '-':
        spin_left_op = '+'
    if spin_right == '+':
        spin_right_op = '-'
    elif spin_right == '-':
        spin_right_op = '+'
    
    file = open('excitation.def').read()
    lines = file.split('\n')[5:-1]
    for line_number in range(len(lines)):
        lines[line_number] = lines[line_number].split() 
        lines[line_number] = [eval(number) for number in lines[line_number]]

    if m == 0:
        ex_state_left = destroy(N,i,spin_left)
    elif m == 1:
        ex_state_left = destroy(N,i,spin_left) @ check(N,i,spin_left_op)
    else:
        lines_left = [ex_state for ex_state in lines if ex_state[1] == i]
        
        ra = lines_left[m][2]
        rb = lines_left[m][3]
             
        ex_state_left = destroy(N,i,spin_left) @ check(N,rb,spin_left) @ check(N,ra,spin_left_op)
    ex_state_left = ex_state_left.transpose().conjugate()
        
    if n == 0:
        ex_state_right = destroy(N,j,spin_right)
    elif n == 1:
        ex_state_right = destroy(N,j,spin_right) @ check(N,j,spin_right_op)
    else:
        lines_right = [ex_state for ex_state in lines if ex_state[1] == j]
       
        ra = lines_right[n][2]
        rb = lines_right[n][3]
        
        ex_state_right = destroy(N,j,spin_right) @ check(N,rb,spin_right) @ check(N,ra,spin_right_op)

    return ex_state_left @ hamiltonian @ ex_state_right

def S_e(N,i,j,m,n,spin_left,spin_right):
    if m < 0 or n < 0:
        raise('Parameter "m" has to be a positive integer.')
    if (spin_right != '+' and spin_right != '-') or (spin_left != '-' and spin_left != '+'):
        raise('Parameters "spin" must be equal to "-" or "+" (both being strings).')
    if spin_left == '+':
        spin_left_op = '-'
    elif spin_left == '-':
        spin_left_op = '+'
    if spin_right == '+':
        spin_right_op = '-'
    elif spin_right == '-':
        spin_right_op = '+'

    file = open('excitation.def').read()
    lines = file.split('\n')[5:-1]
    for line_number in range(len(lines)):
        lines[line_number] = lines[line_number].split() 
        lines[line_number] = [eval(number) for number in lines[line_number]]

    if m == 0:
        ex_state_left = create(N,i,spin_left)
    elif m == 1:
        ex_state_left = create(N,i,spin_left) @ check(N,i,spin_left_op)
    else:
        lines_left = [ex_state for ex_state in lines if ex_state[1] == i]
        
        ra = lines_left[m][2]
        rb = lines_left[m][3]
             
        ex_state_left = create(N,i,spin_left) @ check(N,ra,spin_left_op) @ check(N,rb,spin_left)
    ex_state_left = ex_state_left.transpose().conjugate()
    if n == 0:
        ex_state_right = create(N,j,spin_right)
    elif n == 1:
        ex_state_right = create(N,j,spin_right) @ check(N,j,spin_right_op)
    else:
        lines_right = [ex_state for ex_state in lines if ex_state[1] == j]
       
        ra = lines_right[n][2]
        rb = lines_right[n][3]
        
        ex_state_right = create(N,j,spin_right) @ check(N,ra,spin_right_op) @ check(N,rb,spin_right)
    
    return ex_state_left @ ex_state_right

def S_h(N,i,j,m,n,spin_left,spin_right):
    if m < 0 or n < 0:
        raise('Parameter "m" has to be a positive integer.')
    if (spin_right != '+' and spin_right != '-') or (spin_left != '-' and spin_left != '+'):
        raise('Parameters "spin" must be equal to "-" or "+" (both being strings).')
    if spin_left == '+':
        spin_left_op = '-'
    elif spin_left == '-':
        spin_left_op = '+'
    if spin_right == '+':
        spin_right_op = '-'
    elif spin_right == '-':
        spin_right_op = '+'
    
    file = open('excitation.def').read()
    lines = file.split('\n')[5:-1]
    for line_number in range(len(lines)):
        lines[line_number] = lines[line_number].split() 
        lines[line_number] = [eval(number) for number in lines[line_number]]

    if m == 0:
        ex_state_left = destroy(N,i,spin_left)
    elif m == 1:
        ex_state_left = destroy(N,i,spin_left) @ check(N,i,spin_left_op)
    else:
        lines_left = [ex_state for ex_state in lines if ex_state[1] == i]
        
        ra = lines_left[m][2]
        rb = lines_left[m][3]
             
        ex_state_left = destroy(N,i,spin_left) @ check(N,rb,spin_left) @ check(N,ra,spin_left_op)
    ex_state_left = ex_state_left.transpose().conjugate()
        
    if n == 0:
        ex_state_right = destroy(N,j,spin_right)
    elif n == 1:
        ex_state_right = destroy(N,j,spin_right) @ check(N,j,spin_right_op)
    else:
        lines_right = [ex_state for ex_state in lines if ex_state[1] == j]
        
        ra = lines_right[n][2]
        rb = lines_right[n][3]
        
        ex_state_right = destroy(N,j,spin_right) @ check(N,rb,spin_right) @ check(N,ra,spin_right_op)

    return ex_state_left @ ex_state_right



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


# Calculation of the observable. This is more for testing purposes. Not related to the full matrix.

#result = H_e(2,1,0,0,0,'-','+',Hamiltonian)
result = H_e(2,1,0,3,2,'-','-',Hamiltonian)
#result = H_h(2,1,0,3,2,'-','-',Hamiltonian)


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


# Calculation of one H using the quantum computer. This is for testing purposes.

HHam = JordanWignerMapper.mode_based_mapping(result)

estimator = pEstimator()
job = estimator.run(circuit,HHam)
result = job.result()
values = result.values

print('Result of simulation :', values[0])
print()


# Calculation of the entire matrix using the quantum computer

def H_matrix(N,N_min,excitation_type,spin_left,spin_right,hamiltonian):
    
    # Timer start
    start = time.time()

    if excitation_type != 'e' and excitation_type != 'h':
        raise('Parameter "excitation_type" must be equal to "e" or "h" (both being strings).')
    else:
        func = 'H_' + excitation_type
        H_func = globals()[func]

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
                    column_num = N_exc * i + m
                    row_num = N_exc * j + n
                    if column_num - row_num >= 0:
                        observable = H_func(N,i,j,m,n,spin_left,spin_right,hamiltonian)
                        qubit_hamiltonian = JordanWignerMapper.mode_based_mapping(observable)

                        job = estimator.run(circuit,qubit_hamiltonian)
                        result = job.result()
                        values = result.values 

                        excitation_matrix[row_num,column_num] = values[0]

    excitation_matrix = np.tril(excitation_matrix.T,-1) + excitation_matrix

    end = time.time()
    print('Time:',end - start,'seconds.')

    return excitation_matrix

print(H_matrix(2,1,'h','+','+',Hamiltonian))
print()


def S_matrix(N,N_min,excitation_type,spin_left,spin_right):
    
    # Timer start
    start = time.time()

    if excitation_type != 'e' and excitation_type != 'h':
        raise('Parameter "excitation_type" must be equal to "e" or "h" (both being strings).')
    else:
        func = 'S_' + excitation_type
        S_func = globals()[func]

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
                    column_num = N_exc * i + m
                    row_num = N_exc * j + n
                    if column_num - row_num >= 0:
                        observable = S_func(N,i,j,m,n,spin_left,spin_right)                
                        qubit_hamiltonian = JordanWignerMapper.mode_based_mapping(observable)

                        job = estimator.run(circuit,qubit_hamiltonian)
                        result = job.result()
                        values = result.values 

                        excitation_matrix[row_num,column_num] = values[0]

    excitation_matrix = np.tril(excitation_matrix.T,-1) + excitation_matrix

    end = time.time()
    print('Time:',end - start,'seconds.')

    return excitation_matrix

print(S_matrix(2,1,'h','-','-'))


# Calculation using a normal computer

'''
GS,GSE,bloc =  hs.hubbard(N,t,U,mu)

print(np.matmul(GS.transpose().conjugate(),np.matmul(bloc,GS)))
'''

# Tests
'''
#result = n(2,2,'-') @ n(2,2,'+') @ destroy(2,1,'+') @ Hamiltonian @ create(2,1,'+') @ n(2,2,'+') @ n(2,2,'-')
#result = check(2,0,'-') @ destroy(2,0,'+') @ Hamiltonian @ create(2,1,'-') @ check(2,1,'+')
#result = check(2,0,'-') @ check(2,0,'+') @ destroy(2,1,'-') @ Hamiltonian @ create(2,0,'-') @ check(2,1,'+') @ check(2,1,'-')
#print('result',result.to_matrix().toarray())
'''
