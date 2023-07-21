## Packages ################################################
import numpy as np
import matplotlib.pyplot as plt
import copy,time,sys,os
from mpire import WorkerPool
from qiskit_ibm_runtime import QiskitRuntimeService,Estimator as QC_Estimator
from qiskit_ibm_runtime.options import Options
from qiskit_aer.primitives import Estimator as Noisy_Estimator
from qiskit.primitives import Estimator
from qiskit_nature.second_q.mappers import JordanWignerMapper
from qiskit_nature.second_q.operators import FermionicOp
from qiskit_nature.second_q.hamiltonians import FermiHubbardModel
from qiskit.quantum_info import Statevector


### Fetch parameters.py, hamiltonian_circuit.py, excitdef_reader.py, excitation_directory  ##########
module_directory = os.path.dirname(__file__)
sys.path.insert(0,module_directory)

working_directory = os.getcwd()
sys.path.insert(0,working_directory)

from parameters import (
    N,t,U,mu,spin_green,
    generate_matrix,generate_npy,
    excit_document,output_directory,
    estimator_options,
    noisy_simulation,run_on_quantum_computer,
    token,channel,backend_device,
    force_custom_circuit
)
from hamiltonian_circuit import Hamiltonian, circuit, continue_with_diag

excitation_directory = os.path.join(module_directory,'excitation_files')
sys.path.insert(0,excitation_directory)

from excitdef_reader import excitdef_reader

if os.path.exists(os.path.join(working_directory,excit_document)):
    excitation_directory = os.getcwd()
    print('Using the excitation file in this directory.\n')
else:
    print('Using the excitation.def files included with the code.\n')

from mapping import find_best_layout


### Print options ##########################################
np.set_printoptions(linewidth= 10000,precision=2,suppress=True)


### Print Circuit ##########################################
if N == 2 and force_custom_circuit is False and continue_with_diag is False:
    print('Using hardcoded circuit.')
    print(circuit)


### IBM Service ############################################
if run_on_quantum_computer is True:
    service = QiskitRuntimeService(channel=channel, token=token)
    backend = service.get_backend(backend_device)
else:
    backend = None


### FUNCTIONS ###############################################
## Creation, destruction and check(n and n_dag) operators using qiskit

def create(N,site,spin): 
    """Parameters
    N: Number of sites in total.
    site: Number of the site for which we want to create the operator. The first site is #0.
    spin: Spin to be used for the creation of the operator.

    Returns: Creation operator as a FermionicOp.
    """
    if not spin == '-' and not spin == '+':
        raise Exception('Input must be "+" or "-".')
    if spin == '+':
        spin = 1
    if spin == '-':
        spin = 0

    """The following lines create an FermionicOp, which is a class in qiskit.
    The equation in the parenthesis is simply to accomodate qiskit's notation for
    the fock space.

    | 0,up 0,down 1,up 1,down >
    NOTE: This last line is false? Currently, the notation implemented is | 0,down 0,up 1,down 1,up > and it yields the correct results.
    """
    operator = FermionicOp(
        {
            '+_' + str((2*(N-1-site))+spin): -t,
        },
        num_spin_orbitals=2*N,
        copy=False
    )
    return operator


def destroy(N,site,spin):
    """Parameters
    N: Number of sites in total.
    site: Number of the site for which we want to create the operator. The first site is #0.
    spin: Spin to be used for the creation of the operator.
    
    Returns: Annihilation operator / destroy operator as a FermionicOp.
    """
    return create(N,site,spin).transpose().conjugate()


"""The check operator is simply the name of the operator 'n', which checks if an
electron is there or not. As per qiskit's FermionicOp class, multiplication of 
operators must be written using @. The check operator has 2 variants for if we are
checking for the presence or absence or a fermion, which translates to n and n**dag."""
def check(type,N,site,spin):
    """Parameters
    type: 'presence' or 'absence'
    N: Number of sites in total.
    site: Number of the site for which we want to create the operator. The first site is #0.
    spin: Spin to be used for the creation of the operator.

    Returns: Check operator as a FermionicOp
    """

    if not type.lower() == 'presence' and not type.lower() == 'absence':
        raise Exception('Type must be either "presence" or "absence".')
    if type.lower() == 'presence':
        return create(N,site,spin) @ destroy(N,site,spin)
    if type.lower() == 'absence':
        return destroy(N,site,spin) @ create(N,site,spin)



## Calculation of excited states
def ex_operators(type,i,m,spin,lines_doc):
    """Parameters
    type: which matrix we are calculating (H+,H-,S+ or S-).
    i,m: parameters of the calculated element of the matrix. i is the site numbers and m is the excitation label.
    spin: spin to be used on each side of the H and S matrices equation.
    lines_doc: list form of the excitation.def document. Obtained with the excitdef_reader() function.    

    Returns: The product of the operators associated with an excited state.
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
            ex_ops = c_operator @ check(check_operator,N,rb,spin) @ check(check_operator,N,ra,spin_op) 
    return ex_ops


"""This is the function that calculates the complete observable to be used by the quantum
computer. """

def Observable(type,ex_op_left,ex_op_right,hamiltonian):
    """Parameters
    type: which matrix we are calculating (H+,H-,S+ or S-).
    ex_op_left/right: excited operator calculated using the ex_operator() function.
    hamiltonian: qiskit's FermionicOp object which is, in this case, our hamiltonian.

    Returns: The observable as a FermionicOp.
    """

    """Notice how ex_op_left has .transpose() and .conjugate at its end, because we want the bra, not the ket."""
    ex_op_left = ex_op_left.transpose().conjugate()
    
    """The result is the final observable to be used by the quantum computer."""
    if type[0] == 'H':
        return ex_op_left @ hamiltonian @ ex_op_right
    if type[0] == 'S':
        return ex_op_left @ ex_op_right


## Calculation of the entire matrix using the quantum computer
"""qubit_Observable is the parallelized function."""
def qubit_Observable(hamiltonian,spin,lines_doc,type,i,m,j,n):
    """ Parameters
    hamiltonian: qiskit's FermionicOp object which is, in this case, our hamiltonian.
    spin: spin to be used on each side of the H and S matrices equation.
    lines_doc: list form of the excitation.def document. Obtained with the excitdef_reader() function.    
    type: which matrix we are calculating (H+,H-,S+ or S-).
    i,m,j,n: parameters of the calculated element of the matrix. i,j are site numbers are m,n are excitation labels.

    Returns: The observable in the Jordan-Wigner basis.
    """

    """One downside of the mpire parallelization is that it 
    requires importing the packages in the target function."""
    from qiskit_nature.second_q.mappers import JordanWignerMapper
    from qiskit_nature.second_q.operators import FermionicOp
    import numpy as np

    ex_op_left = ex_operators(type,i,m,spin,lines_doc)
    ex_op_right = ex_operators(type,j,n,spin,lines_doc)
    observable = Observable(type,ex_op_left,ex_op_right,hamiltonian)
    qubit_observable = JordanWignerMapper.mode_based_mapping(observable)

    return qubit_observable


def matrix_observables(type,lines_doc,N,spin,hamiltonian):
    """ Parameters
    hamiltonian: qiskit's FermionicOp object which is, in this case, our hamiltonian.
    spin: spin to be used on each side of the H and S matrices equation.
    lines_doc: list form of the excitation.def document. Obtained with the excitdef_reader() function.    
    type: which matrix we are calculating (H+,H-,S+ or S-).
    q_circuit: qiskit's quantum circuit object. In this case, it is our ground state.
    save: wether we save or not the matrices as .npy files. (Y or N)

    Returns: List of observables as FermionicOp objects for the specified matrix type.
    """
    
    # This is the total possible number of excitation for each site, as explained by the paper.
    lines = [values for values in lines_doc if values[1] == 0]
    N_exc = len(lines)


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
                            param.append((spin,lines_doc,type,i,m,j,n))

        '''Using the list of parameters, we use the defined pool to queue all the calculations at once.
        Each element of the matrix is equivalent to one process. When a process is done, the pool
        automatically asigns a new task. The output is a list with the results being in the same order
        as the param list. The reason the mpire module is used instead of the included multiprocessing
        module is because mpire is easier to use and supports sharing objects with all the processes
        (in this case the hamiltonian). As a nice bonus, we have a progress bar.'''
        observables = pool.map(qubit_Observable,param,progress_bar=True)

    return observables


def job(observables,circuit,backend,noisy_simulation=False,run_on_quantum_computer=False):
    """Parameters
    observables: Observables obtained with matrix_observables(). List of FermionicOps.
    circuit: Single circuit to be repeated for each observables. Circuit object from qiskit.
    backend: Backend for the calculation. Corresponds to the ibm device.
    noisy_simulation: Boolean that dictates if the job should be a quantum simulation with added noise using Qiskit's Aer module.
    run_on_quantum_computer: Boolean that dictates if the job should be ran on an actual quantum computer instead of Qiskit's simulators.

    Returns: List of values corresponding to the results of the Estimator.
    """

    # Defining estimator
    if noisy_simulation is True:
        estimator = Noisy_Estimator(backend_options=estimator_options)
    elif run_on_quantum_computer is True:
        init_layout,*rest = find_best_layout(circuit=circuit,backend=backend,num_tries=10,seed=50)
        options = Options()
        options.transpilation.initial_layout=init_layout
        estimator = QC_Estimator(backend = backend,options=options)
    else:
        estimator = Estimator()

    # Running job
    job = estimator.run([circuit]*len(observables),observables)

    if run_on_quantum_computer is True:
        print(f">>> Job ID: {job.job_id()}")
        print(f">>> Job Status: {job.status()}")
    else:
        print('Quantum Computer simulation...')

    # Fetching results
    result = job.result()
    values = result.values # This outputs all the values of the matrix in the order they were calculated above.
          # This order is line by line, from left to right, ommiting the elements we are not calculating.

    return values


def matrix(type,N,values,save=True):
    """Parameters
    type: H+, H-, S+ or S-
    N: Number of sites
    values: Values obtained from the job function defined above.
    save: If the matrices are to be saved as .npy files.

    Returns: Matrix of the specififed type.
    """

    # Find N_exc using the list of observables
    number_of_values = 0
    lenght_new_row = 0    
    while not number_of_values == len(values):
        lenght_new_row += 1
        number_of_values += lenght_new_row
    N_exc = int(lenght_new_row/N)
    
    # Initializing matrix
    matrix_size = N * N_exc
    excitation_matrix = np.zeros((matrix_size,matrix_size))

    # Filling matrix
    """This strange loop is because the list generated above makes it difficult to assign the values at the right index in the matrix.
    Basically, we have to figure out these indexes:
    [i i i i]
    [- i i i]
    [- - i i]
    [- - - i]
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

    if save is True:
        if type[1] == '+':
            identifier = '_AC'
        if type[1] == '-':
            identifier = '_CA'

        np.save(os.path.join(output_directory,type[0]+identifier),excitation_matrix)

    return excitation_matrix

if __name__ == '__main__':
    lines_doc = excitdef_reader(excit_document,excitation_directory) 
    N_exc = len([values for values in lines_doc if values[1] == 0])

    start = time.time() 

    if generate_matrix.upper() == 'ALL':
        # Grouping all observables under one job
        observables_all = []
        for type in ['H+','H-','S+','S-']:
            observables = matrix_observables(type,lines_doc,N,spin_green,Hamiltonian)
            observables_all.extend(observables)
        
        # Submitting job
        x = len(observables)
        values = job(
            observables=observables_all,
            circuit=circuit,
            backend=backend,
            noisy_simulation=noisy_simulation,
            run_on_quantum_computer=run_on_quantum_computer
            )

        # Extracting the proper matrix values from the output of job() depending on the type of the matrix
        # Calculating matrix
        controller = 0
        for type in ['H+','H-','S+','S-']:
            print('##### '+type+' #####')
            matrix_values = values[x*controller:x*(controller+1)]
            controller += 1
            print(matrix(type,N,matrix_values))
            print()
        
        end = time.time()
        print('Time:',end-start,'seconds.')

    else:
        print('##### '+generate_matrix+' #####')

        observables = matrix_observables(generate_matrix,lines_doc,N,spin_green,Hamiltonian)
        values = job(
            observables=observables,
            circuit=circuit,
            backend=backend,
            noisy_simulation=noisy_simulation,
            run_on_quantum_computer=run_on_quantum_computer
            )
        print(matrix(type,N,values))

        end = time.time()
        print('Time:',end-start,'seconds.')
