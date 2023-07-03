### Packages ################################################
import numpy as np
import matplotlib.pyplot as plt
import copy
import time
import sys,os
from mpire import WorkerPool
from alive_progress import alive_bar
from qiskit.quantum_info import Pauli,Operator
from qiskit.primitives import Estimator as pEstimator
from qiskit_nature.second_q.mappers import JordanWignerMapper
from qiskit_nature.second_q.operators import FermionicOp
from qiskit import QuantumCircuit,QuantumRegister
from qiskit import transpile
from qiskit_nature.second_q.hamiltonians import FermiHubbardModel
from qiskit.tools.monitor import job_monitor
from qiskit.providers.ibmq.job import IBMQJob
from qiskit_nature.second_q.hamiltonians.lattices import (
    BoundaryCondition,
    HyperCubicLattice,
    Lattice,
    LatticeDrawStyle,
    LineLattice,
    SquareLattice,
    TriangularLattice,
)

### Fetch parameters.py, hamiltonian_circuit.py, excitdef_reader.py, excitation_directory  ##########
module_directory = os.path.dirname(__file__)
sys.path.insert(0,module_directory)

working_directory = os.getcwd()
sys.path.insert(0,working_directory)

from parameters import N,t,U,mu,generate_matrix,excit_document,spin_left,spin_right,generate_npy,output_directory,excit_document
from hamiltonian_circuit import Hamiltonian, circuit

excitation_directory = os.path.join(module_directory,'excitation_files')
sys.path.insert(0,excitation_directory)

from excitdef_reader import excitdef_reader

if os.path.exists(os.path.join(working_directory,excit_document)):
    excitation_directory = os.getcwd()
    print('Using the excitation file in this directory.\n')
else:
    print('Using included excitation.def files.\n')


# Print options
np.set_printoptions(linewidth= 10000)

### FUNCTIONS ###############################################
## Creation, destruction and check(n and n_dag) operators using qiskit

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
operators must be written using @. The check operator has 2 variants for if we are
checking for the presence or absence or a fermion."""
def check(type,N,site,spin):
    if not type == 'presence' and not type == 'absence':
        raise Exception('Type must be either "presence" or "absence".')
    if type == 'presence':
        return create(N,site,spin) @ destroy(N,site,spin)
    if type == 'absence':
        return destroy(N,site,spin) @ create(N,site,spin)


## Calculation of any excited state
"""This function calculates the product of the operators associated with one exited state.
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

    """This code block specifies the 'c' operator and 'n' operator to be used in, for example, |e> = cnn|GS>.
    In the last example, the c should be a c_dagger and check_operator should be a standard 'presence' operator."""
    if type[1] == '+':
        c_operator = create(N,i,spin)
        check_operator = 'presence'
    if type[1] == '-':
        c_operator = destroy(N,i,spin)
        check_operator = 'absence'
    
    """For the spin defined, we also define its opposite."""
    if spin == '+':
        spin_op = '-'
    if spin == '-':
        spin_op = '+'
    
    """For each definition of an excited state, we define the operation to be performed."""
    if m == 0:
        ex_ops = c_operator
    else:
        lines = [values for values in lines_doc if values[1] == i]
        t = lines[m][0]
        
        """Here, we fetch the values of ra and rb from the lines_doc we inputed in the function."""
        ra = lines[m][2]
        rb = lines[m][3]

        if t == 1:
            ex_ops = c_operator @ check(check_operator,N,ra,spin_op)
        elif t == 3:
            ex_ops = c_operator @ check(check_operator,N,ra,spin_op) @ check(check_operator,N,rb,spin_op)
        else:
            ex_ops = c_operator @ check(check_operator,N,ra,spin_op) @ check(check_operator,N,rb,spin)
    return ex_ops


"""This is the function that calculates the observable to be used by the quantum
computer. """

def Observable(type,ex_op_left,ex_op_right,hamiltonian):
    '''Notice how ex_op_left has .transpose() and .conjugate at its end, because we want the bra, not the ket.'''
    ex_op_left = ex_op_left.transpose().conjugate()
    
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

def qubit_Observable(hamiltonian,spin_left,spin_right,lines_doc,type,i,m,j,n):
    from qiskit_nature.second_q.mappers import JordanWignerMapper
    from qiskit_nature.second_q.operators import FermionicOp
    import numpy as np
    ex_op_left = ex_operators('left',type,i,m,spin_left,lines_doc)
    ex_op_right = ex_operators('right',type,j,n,spin_right,lines_doc)
    observable = Observable(type,ex_op_left,ex_op_right,hamiltonian)
    qubit_hamiltonian = JordanWignerMapper.mode_based_mapping(observable)
    #excitation_matrix[row_num,column_num] = pEstimator().run(q_circuit,qubit_hamiltonian).result().values[0]
    #observables.append(qubit_hamiltonian)
    return qubit_hamiltonian

def matrix(type,lines_doc,N,spin_left,spin_right,hamiltonian,q_circuit,save='N'):
    
    # This is the total possible number of excitation for each site, as explained by the paper.
    lines = [values for values in lines_doc if values[1] == 0]
    N_exc = len(lines)

    # Timer start
    start = time.time() 

    # Initializing matrix
    matrix_size = N * N_exc
    excitation_matrix = np.zeros((matrix_size,matrix_size))

    # Caluculation of half the elements 
    circuits = [q_circuit] * (N * N_exc)**2
    
    print('Observables calculation...')
    with WorkerPool(n_jobs=None,shared_objects=(hamiltonian)) as pool:
        param = []
        for n in range(N_exc):
            for j in range(N):
                for m in range(N_exc):
                    for i in range(N): 
                        """The 2 following lines represents the distribution of the calculated values in the
                        matrix. For each m, there are (for 2 sites) 2 values of i. Thus, the first 2 columns
                        would be for m = 0 and for 0 <= i <= 1, then the next 2 for m = 1 and for 0 <= i <= 1.
                        The same principles for the rows, but i becomes j and m becomes n."""
                        column_num = N * m + i
                        row_num = N * n + j

                        if column_num - row_num >= 0: # This is to calculate only half of the matrix.
                            param.append((spin_left,spin_right,lines_doc,type,i,m,j,n))

        observables = pool.map(qubit_Observable,param,progress_bar=True)
                            

    # Starting quantum simulation
    job = pEstimator().run([q_circuit]*int((1/2)*N*N_exc*(N*N_exc+1)),observables)
    print('Quantum Computer simulation...')
    result = job.result()
    values = result.values # This outputs all the values of the matrix in the order they were calculated above.
          # This order is line by line, from left to right, ommiting the elements we are not calculating.

    # Fill matrix
    # This strange loop is because the list generated above makes it difficult to assign the values at the right index in the matrix.
    correction = 0
    row = 0
    for column,value in enumerate(values):
        column = column + correction
        excitation_matrix[row,column] = value
        if column == N*N_exc-1:
            row += 1
            correction += -N*N_exc + row
    
    # Symmetrizing the matrix
    excitation_matrix = np.tril(excitation_matrix.T,-1) + excitation_matrix

    end = time.time()
    print('Time:',end - start,'seconds.')

    if save.upper() == 'Y':
        if type[1] == '+':
            identifier = '_AC'
        if type[1] == '-':
            identifier = '_CA'

        np.save(os.path.join(output_directory,type[0]+identifier),excitation_matrix)

    return excitation_matrix

lines_doc = excitdef_reader(excit_document,excitation_directory)
if generate_matrix.upper() == 'ALL':
    for type in ['H+','H-','S+','S-']:
        print('##### '+type+' #####')
        print(matrix(type,lines_doc,N,spin_left,spin_right,Hamiltonian,circuit,generate_npy))
        print()
else:
    print('##### '+generate_matrix+' #####')
    print(matrix(generate_matrix,lines_doc,N,spin_left,spin_right,Hamiltonian,circuit,generate_npy))
