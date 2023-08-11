#!/usr/bin/env python3
### Packages ################################################
import numpy as np
from mpire import WorkerPool
import copy
import time
import sys,os

### Fetch parameters.py, hamiltonian_circuit.py, excitdef_reader.py, excitation_directory,fock_class.py and hubbard_classes.py  ##########
module_directory = os.path.dirname(__file__)
sys.path.insert(0,module_directory)

working_directory = os.getcwd()
sys.path.insert(0,working_directory)

from parameters import N,t,hopping_matrix,U,mu,generate_matrix,excit_document,spin_green,spin_gs,generate_npy,output_directory,pdf_output_directory
from hamiltonian_circuit import hopping_matrix


excitation_directory = os.path.join(module_directory,'excitation_files')
sys.path.insert(0,excitation_directory)

from excitdef_reader import excitdef_reader

if os.path.exists(os.path.join(working_directory,excit_document)):
    excitation_directory = os.getcwd()
    print('Using the excitation file in this directory.\n')
else:
    print('Using included excitation.def files.\n')


sys.path.insert(0,os.path.join(module_directory,'..','second_quantization_codes'))
import fock_class as f
from hubbard_classes import hubbard


# Print options
np.set_printoptions(linewidth = 1000,precision=2,suppress=True)


### FUNCTIONS ################################################
# In fock's, returns the excited state.
def ex_state(type,i,m,spin,gs_block_hub,gs_numerical_state,lines_doc):
    """Parameters
    type: H+,H-,S+ or S-.
    i,m: i is the site, m the excitation label.
    spin: spin to be used for the calculation.
    gs_block_hub: the list of fock's basis states for the GS, obtained with the hubbard() function when manip='yes'.
    gs_numerical_state: actual ground state of the system, obtained with the hubbard() function.

    Returns: List of Fock objects representing the excited state.
    """

    # Assigns the components of the numerical ground state to its respective state in the block
    excited_block = [[state] for state in copy.deepcopy(gs_block_hub)]
    for index,state in enumerate(excited_block):
        state.append(gs_numerical_state[index])

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

    """Defining which operators to use depending on the matrix we are calculating.
    Since calculating the '-' matrices involves checking for a 'hole', the order
    in which the sub-operators of the operator 'n' are calculated is different
    (creation before destruction instead of destruction before creation).
    Hence, we have n1_operator and n2_operator to control that order.

    c**dag c**dag c |GS> is c**dag n1_operator n2_operator |GS>
    """
    if type[1] == '+':
        c_operator = 'create'
        n1_operator = 'create'
        n2_operator = 'destroy'
    if type[1] == '-':
        c_operator = 'destroy'
        n1_operator = 'destroy'
        n2_operator = 'create'
    
    # Defining the opposite of the inputted spin
    if spin == '+':
        spin_op = '-'
    if spin == '-':
        spin_op = '+'
    
    # Calculating the excited state using the excitation.def file.
    for state in excited_block:
        lines = [x for x in lines_doc if x[1] == i] 
        
        # Fetching important values from the excitation.def
        t = lines[m][0]        
        ra = lines[m][2]
        rb = lines[m][3]

        # Calculation depending on the nature of the excited state.
        if m == 0:
            state[0].op(c_operator,i,spin)
        elif t == 1:
            state[0].op(n2_operator,ra,spin_op)
            state[0].op(n1_operator,ra,spin_op)
            state[0].op(c_operator,i,spin)
        elif t == 3:
            state[0].op(n2_operator,ra,spin_op)
            state[0].op(n1_operator,ra,spin_op)
            state[0].op(n2_operator,rb,spin_op)
            state[0].op(n1_operator,rb,spin_op)
            state[0].op(c_operator,i,spin)
        else:
            state[0].op(n2_operator,ra,spin_op)
            state[0].op(n1_operator,ra,spin_op)
            state[0].op(n2_operator,rb,spin)
            state[0].op(n1_operator,rb,spin)
            state[0].op(c_operator,i,spin)
    
    # Cleaning up the excited_state to only include the actual states that survived the operations
    excited_block[:] = [state for state in excited_block if not isinstance(state[0].fock,int)]
    
    return excited_block

def choose_gs(gs_blocks,gs_numerical_states,spin_gs):
    """Parameters
    gs_blocks: List of multiple gs blocks containing each multiple states. Obtained with hubbard().
    gs_numerical_states: GS vector. Obtained with hubbard().
    spin_gs: Spin priority of the ground state.

    Returns: gs_block and gs_numerical_state, who have been chosen based on their total spin value. 
    """
    final_spin = gs_blocks[0][0].total_spin
    final_index = 0
    for index,block in enumerate(gs_blocks):
        new_spin = block[0].total_spin
        if ((spin_gs == '+' and 0 <= new_spin and new_spin < final_spin)
        or (spin_gs == '-' and 0 >= new_spin and new_spin > final_spin)):
            final_spin = new_spin
            final_index = index
    gs_block = gs_blocks[final_index]
    gs_numerical_state = gs_numerical_states[final_index]
    return gs_block,gs_numerical_state

def element(type,N,ex_state_left,ex_state_right,hubbard_output):
    """Parameters
    type: H+, H-, S+ or S-
    N: Number of sites
    ex_state_left/right: List of Fock objects that corresponds to the states that survived the n operators and the c operator.
    hubbard_ouput: The outputs of the hubbard(manip=True) function.

    Returns: Single element of the matrix we wish to calculate.
    """

    # Defining the outputs of the hubbard() function
    blocks_matrix = hubbard_output[0]
    blocks_num = hubbard_output[1]
    blocks = hubbard_output[3]
    gs_block = hubbard_output[2]
    gs_numerical_state = hubbard_output[4]
    
    # Initializing a list that will contain the results of the distribution of the multiplication of the excited states with the hamiltonian
    scalar = []
    
    # Distributing multiplication
    if ex_state_right and ex_state_left:
        for state_left in ex_state_left:
            for state_right in ex_state_right:
                
                """If we are calculating the S matrix, no need to include the hamiltonian. If the states are the same, their product is one.
                We only calculate the product of the components of the numerical gs contained in state_left and state_right, and also the product
                of the carried signs because of the permuttations of the operators in fock's base."""
                if type[0] == 'S':
                    if np.array_equal(state_left[0].fock,state_right[0].fock):
                        value = state_left[1]*state_left[0].sign*state_right[1]*state_right[0].sign
                        scalar.append(value)

                """If we are calculating the H matrices, we need to calculate <state_left|H|state_right> for each element in ex_state_right/left. The result
                of this operation is located in the hamiltonian itself. Thus we need to look into it. Finding the indexes is troublesome, hence why it is important
                to understand how the hubbard function outputs's lists are structed. For example:
                [[N=0],[N=1],[N=2],...]
                [N=2] = [[Sz=-1],[Sz=0],...]
                """
                if type[0] == 'H':
                    # This 'if' determines if the states are in the same block. If they are not, no point in finding the indexes because the result is 0.
                    if (state_left[0].num_electrons == state_right[0].num_electrons
                    and state_left[0].total_spin == state_right[0].total_spin):

                        # The states must share the same number of electrons and the same total spin. The following lines figure out these two values.
                        n_electrons = state_left[0].num_electrons
                        for block_num in blocks_num[n_electrons]:
                            if ex_state_left[0][0].num in block_num:
                                index_sp = blocks_num[n_electrons].index(block_num)
                        
                        """Using the total spin and num of electrons, the bloc containing the two states is located. We must now find what are the python indexes
                        or in other words, the row and the column of the desired value."""
                        row = blocks_num[n_electrons][index_sp].index(state_left[0].num)
                        column = blocks_num[n_electrons][index_sp].index(state_right[0].num)

                        """Upon finding the desired value in the hamitlonian, we calculate the final value/result, which is the product of the gs numerical components, 
                        the signs and the value found in the hamiltonian."""
                        value = blocks_matrix[n_electrons][index_sp][row,column]*state_left[1]*state_left[0].sign*state_right[1]*state_right[0].sign
                        scalar.append(value)

        sum_values = sum(scalar)
        return sum_values
    else:
        return 0


def parrallelized_element(shared_lists,spin,type,i,m,j,n):
    """Parameters
    shared_lists: Obtained with mpire. Corresponds to the parameters we are sharing for all processes.
    spin: Spin for which we are calculating the Green function.
    type: H+, H-, S+ or S-
    i/j: Number of the site
    m/n: Label of the excitation

    Returns: Single element of the matrix we wish to calculate.
    """
    import numpy as np

    ex_state_left = ex_state(type,i,m,spin,shared_lists[1],shared_lists[2],shared_lists[3])
    ex_state_right = ex_state(type,j,n,spin,shared_lists[1],shared_lists[2],shared_lists[3])
    ele = element(type,N,ex_state_left,ex_state_right,shared_lists[0])
    return ele
 

def matrix(type,lines_doc,N,spin,spin_gs,t,hopping_matrix,U,mu,save):
    """Parameters
    type: H+, H-, S+ or S-
    lines_doc: List form of the excitation.def file. Obtained with excitdef_reader().
    N: Number of sites
    spin: Spin for which we are calculating the Green function.
    spin_gs: Priority total spin of the ground state.
    t: Hopping energy.
    hopping_matrix: Hopping matrix that states the allowed jumps. 1 is a possible jump, 0 is not a possible jump. Diagonal should be 0.
    U: Potential energy.
    mu: Chemical potential energy.
    generate_npy: Boolean that dictates if the matrices should be saved or not as .npy files.

    Returns: Matrix of the chosen type.
    """
    
    # This is the total possible number of excitation for each site, as explained by the paper.
    lines = [values for values in lines_doc if values[1] == 0]
    N_exc = len(lines)
    start = time.time()

    # Creating empty matrix
    matrix_size = N * N_exc
    excitation_matrix = np.zeros((matrix_size,matrix_size))
    
    # Hubbard output
    hubbard_output = list(hubbard(N,t,hopping_matrix,U,mu,spin_gs=spin_gs,manip=True,qis_not=True))
    gs_blocks = hubbard_output[2]
    gs_numerical_states = hubbard_output[4]
    hubbard_output[2],hubbard_output[4] = choose_gs(gs_blocks,gs_numerical_states,spin_gs)
    gs_block = hubbard_output[2]
    gs_numerical_state = hubbard_output[4]

    # Filling half of the matrix
    shared_lists = (hubbard_output,gs_block,gs_numerical_state,lines_doc)
    with WorkerPool(n_jobs=None,shared_objects=(shared_lists)) as pool:
        param = []
        for n in range(N_exc):
            for j in range(N):
                for m in range(N_exc):
                    for i in range(N):
                        column_num = N * m + i
                        row_num = N * n + j

                        if column_num - row_num >= 0: # This is to calculate only half of the matrix, including the diagonal.
                            param.append((spin,type,i,m,j,n))

        values = pool.map(parrallelized_element,param,progress_bar=True)

    # Filling matrix
    correction = 0
    row = 0
    for column,value in enumerate(values):
        column = column + correction
        excitation_matrix[row,column] = value
        if column == N*N_exc-1:
            row += 1
            correction += -N*N_exc + row

    # Filling all of the matrix
    excitation_matrix = np.tril(excitation_matrix.T,-1) + excitation_matrix    

    end = time.time()
    print('Time:',end-start,'seconds.')

    if type[1] == '+':
        identifier = '_AC'
    if type[1] == '-':
        identifier = '_CA'
    
    if save is True:
        np.save(os.path.join(output_directory,type[0]+identifier+'_fock'),excitation_matrix)
    
    return excitation_matrix

if __name__ == "__main__":
    omega = hubbard(N,t,hopping_matrix,U,mu,spin_gs=spin_gs,qis_not=True)
    lines_doc = excitdef_reader(excit_document,excitation_directory)
    if generate_matrix.upper() == 'ALL':
        for type in ['H+','H-','S+','S-']:
            print(type+':')
            print(matrix(type,lines_doc,N,spin_green,spin_gs,t,hopping_matrix,U,mu,generate_npy))
            print()
    else:
        print(generate_matrix+':')
        print(matrix(generate_matrix,lines_doc,N,spin_green,spin_gs,t,hopping_matrix,U,mu,generate_npy))
        

    # Generating graph using graph.py
    from graph import dvmc_spectrum
    dvmc_spectrum(omega[0],verbose=1,fock_benchmarking=True)
