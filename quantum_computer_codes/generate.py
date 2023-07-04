## Packages ################################################
import numpy as np
import matplotlib.pyplot as plt
import copy,time,sys,os
from mpire import WorkerPool
from qiskit.primitives import Estimator
from qiskit_nature.second_q.mappers import JordanWignerMapper
from qiskit_nature.second_q.operators import FermionicOp
from qiskit_nature.second_q.hamiltonians import FermiHubbardModel



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
    print('Using the excitation.def files included with the code.\n')



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
checking for the presence or absence or a fermion, which translates to n and n**dag."""
def check(type,N,site,spin):
    if not type == 'presence' and not type == 'absence':
        raise Exception('Type must be either "presence" or "absence".')
    if type == 'presence':
        return create(N,site,spin) @ destroy(N,site,spin)
    if type == 'absence':
        return destroy(N,site,spin) @ create(N,site,spin)


## Calculation of excited states
"""This function calculates the product of the operators associated with one exited state.
It is important to note that this doesn't calculate the product of those operators with the GS.
This task is for the quantum computer."""
def ex_operators(type,i,m,spin,lines_doc):
    """Parameters
    type: which matrix we are calculating (H+,H-,S+ or S-).
    i,m: parameters of the calculated element of the matrix. i is the site numbers and m is the excitation label.
    spin: spin to be used on each side of the H and S matrices equation.
    lines_doc: list form of the excitation.def document. Obtained with the excitdef_reader() function.    
    """
    
    """The following block merely ensure proper values are inserted in the function."""
    valid_type = ['H+','H-','S+','S-']
    valid_spin = ['+','-']
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
        """We isolate the lines that correspond to i in lines_doc."""
        lines = [values for values in lines_doc if values[1] == i]
        
        """Here, we fetch the values of ra,rb and t from the lines_doc we inputed in the function."""
        t = lines[m][0]
        ra = lines[m][2]
        rb = lines[m][3]

        if t == 1:
            ex_ops = c_operator @ check(check_operator,N,ra,spin_op)
        elif t == 3:
            ex_ops = c_operator @ check(check_operator,N,ra,spin_op) @ check(check_operator,N,rb,spin_op)
        else:
            ex_ops = c_operator @ check(check_operator,N,ra,spin_op) @ check(check_operator,N,rb,spin)
    return ex_ops


"""This is the function that calculates the complete observable to be used by the quantum
computer. """

def Observable(type,ex_op_left,ex_op_right,hamiltonian):
    """Parameters
    type: which matrix we are calculating (H+,H-,S+ or S-).
    ex_op_left/right: excited operator calculated using the ex_operator() function.
    hamiltonian: qiskit's FermionicOp object which is, in this case, our hamiltonian.
    """

    '''Notice how ex_op_left has .transpose() and .conjugate at its end, because we want the bra, not the ket.'''
    ex_op_left = ex_op_left.transpose().conjugate()
    
    """The result is the final observable to be used by the quantum computer."""
    if type[0] == 'H':
        return ex_op_left @ hamiltonian @ ex_op_right
    if type[0] == 'S':
        return ex_op_left @ ex_op_right


## Calculation of the entire matrix using the quantum computer
"""qubit_Observable is the parallelized function."""
def qubit_Observable(hamiltonian,spin_left,spin_right,lines_doc,type,i,m,j,n):
    """ Parameters
    hamiltonian: qiskit's FermionicOp object which is, in this case, our hamiltonian.
    spin_left/right: spin to be used on each side of the H and S matrices equation.
    lines_doc: list form of the excitation.def document. Obtained with the excitdef_reader() function.    
    type: which matrix we are calculating (H+,H-,S+ or S-).
    i,m,j,n: parameters of the calculated element of the matrix. i,j are site numbers are m,n are excitation labels.
    """

    """One downside of the mpire parallelization is that it 
    requires importing the packages in the target function."""
    from qiskit_nature.second_q.mappers import JordanWignerMapper
    from qiskit_nature.second_q.operators import FermionicOp
    import numpy as np

    ex_op_left = ex_operators(type,i,m,spin_left,lines_doc)
    ex_op_right = ex_operators(type,j,n,spin_right,lines_doc)
    observable = Observable(type,ex_op_left,ex_op_right,hamiltonian)
    qubit_observable = JordanWignerMapper.mode_based_mapping(observable)

    return qubit_observable

def matrix(type,lines_doc,N,spin_left,spin_right,hamiltonian,q_circuit,save='N'):
    """ Parameters
    hamiltonian: qiskit's FermionicOp object which is, in this case, our hamiltonian.
    spin_left/right: spin to be used on each side of the H and S matrices equation.
    lines_doc: list form of the excitation.def document. Obtained with the excitdef_reader() function.    
    type: which matrix we are calculating (H+,H-,S+ or S-).
    q_circuit: qiskit's quantum circuit object. In this case, it is our ground state.
    save: wether we save or not the matrices as .npy files. (Y or N)
    """
    
    # This is the total possible number of excitation for each site, as explained by the paper.
    lines = [values for values in lines_doc if values[1] == 0]
    N_exc = len(lines)

    # Timer start
    start = time.time() 

    # Initializing matrix
    matrix_size = N * N_exc
    excitation_matrix = np.zeros((matrix_size,matrix_size))

    print('Observables calculation...')
    # This loop creates a list of the parameters to be used for each element in the matrix
    with WorkerPool(n_jobs=None,shared_objects=(hamiltonian)) as pool:
        param = []
        for n in range(N_exc):
            for j in range(N):
                for m in range(N_exc):
                    for i in range(N): 
                        """The 2 following lines represents the distribution of the calculated values in the
                        matrix. For each m, there are (for 2 sites) 2 values of i. Thus, the first 2 columns
                        would be for m = 0 and for 0 <= i <= 1, then the next 2 for m = 1 and for 0 <= i <= 1.
                        The same logic applies to the rows, but i becomes j and m becomes n."""
                        column_num = N * m + i
                        row_num = N * n + j

                        if column_num - row_num >= 0: # This is to calculate only half of the matrix, including the diagonal.
                            param.append((spin_left,spin_right,lines_doc,type,i,m,j,n))

        '''Using the list of parameters, we use the defined pool to queue all the calculations at once.
        Each element of the matrix is equivalent to one process. When a process is done, the pool
        automatically asigns a new task. The output is a list with the results being in the same order
        as the param list. The reason the mpire module is used instead of the included multiprocessing
        module is because mpire is easier to use and supports sharing objects with all the processes
        (in this case the hamiltonian). As a nice bonus, we have a progress bar.'''
        observables = pool.map(qubit_Observable,param,progress_bar=True)

    # Starting quantum simulation
    job = Estimator().run([q_circuit]*int((1/2)*N*N_exc*(N*N_exc+1)),observables)
    print('Quantum Computer simulation...')
    result = job.result()
    values = result.values # This outputs all the values of the matrix in the order they were calculated above.
          # This order is line by line, from left to right, ommiting the elements we are not calculating.
    
    # Filling matrix
    """This strange loop is because the list generated above makes it difficult to assign the values at the right index in the matrix.
    Basically, we have to figure out these indexes:
    [0 0 0 0]
    [- 0 0 0]
    [- - 0 0]
    [- - - 0]
    """
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

    if save.upper() == 'Y':
        if type[1] == '+':
            identifier = '_AC'
        if type[1] == '-':
            identifier = '_CA'

        np.save(os.path.join(output_directory,type[0]+identifier),excitation_matrix)

    end = time.time()
    print('Time:',end - start,'seconds.')
    
    return excitation_matrix



if __name__ == '__main__':
    
    lines_doc = excitdef_reader(excit_document,excitation_directory)
    if generate_matrix.upper() == 'ALL':
        for type in ['H+','H-','S+','S-']:
            print('##### '+type+' #####')
            print(matrix(type,lines_doc,N,spin_left,spin_right,Hamiltonian,circuit,generate_npy))
            print()
    else:
        print('##### '+generate_matrix+' #####')
        print(matrix(generate_matrix,lines_doc,N,spin_left,spin_right,Hamiltonian,circuit,generate_npy))
