### Packages ################################################
import numpy as np
import copy
import time


### Fetching parameters.py, hamiltonian_circuit.py, fock_class.py, hubbard_classes.py and excitdef_reader ###
import sys,os
if len(sys.argv) == 2:
    number = sys.argv[1]
    sys.path.insert(0,os.path.join(os.path.dirname(__file__),'examples',number+'sites'))
sys.path.insert(0,os.path.join(os.path.dirname(__file__),'..','second_quantization_codes'))

from parameters import N,t_fock,U,mu,generate_matrix,excit_document,spin_left,spin_right,generate_npy,output_directory,pdf_output_directory,excitation_directory
from hamiltonian_circuit import t_fock
from excitdef_reader import excitdef_reader
import fock_class as f
import hubbard_classes as h
import hubbard_classes_destroy as hd


# Print options
np.set_printoptions(linewidth = 1000,precision=4)


### FUNCTIONS
# In fock's, returns the excited state.
def ex_state(type,i,m,spin,gs_block_hub,gs_numerical_state,lines_doc):

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
    Hence, we have n1_operator and n2_operator to control that order."""
    if type[1] == '+':
        c_operator = 'create'
        n1_operator = 'create'
        n2_operator = 'destroy'
    if type[1] == '-':
        c_operator = 'destroy'
        n1_operator = 'destroy'
        n2_operator = 'create'
        #n1_operator = 'create'
        #n2_operator = 'destroy'
    
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
        else:
            if t == 1:
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

def element(type,N,ex_state_left,ex_state_right,hubbard_output,lines_doc):

    # Defining the outputs of the hubbard() function
    blocks_matrix = hubbard_output[0]
    blocks_num = hubbard_output[1]
    gs_block = hubbard_output[2]
    blocks = hubbard_output[3]
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
                to understand how the hubbard outputs's lists are structed."""
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


def matrix(type,lines_doc,N,spin_left,spin_right,t_mat_,U,mu,generate_npy):
    
    # This is the total possible number of excitation for each site, as explained by the paper.
    lines = [values for values in lines_doc if values[1] == 0]
    N_exc = len(lines)
    start = time.time()

    # Creating empty matrix
    matrix_size = N * N_exc
    excitation_matrix = np.zeros((matrix_size,matrix_size))
    
    # Hubbard output
    if type[1] == '+':
        hubbard_output = h.hubbard(N,t_mat_,U,mu,'yes',qis_not='Y')
    else:
        hubbard_output = hd.hubbard(N,t_mat_,U,mu,'yes',qis_not='Y')
    gs_block = hubbard_output[2]
    gs_numerical_state = hubbard_output[4]
    
    # Filling half of the matrix
    for i in range(N):
        for j in range(N):
            for m in range(N_exc):
                for n in range(N_exc): 
                    column_num = N * m + i
                    row_num = N * n + j
                    if column_num - row_num >= 0:
                        ex_state_left = ex_state(type,i,m,spin_left,gs_block,gs_numerical_state,lines_doc)
                        ex_state_right = ex_state(type,j,n,spin_right,gs_block,gs_numerical_state,lines_doc)
                        excitation_matrix[row_num,column_num] = element(type,N,ex_state_left,ex_state_right,hubbard_output,excit_document)
    
    # Filling all of the matrix
    excitation_matrix = np.tril(excitation_matrix.T,-1) + excitation_matrix    

    end = time.time()
    print('Time:',end-start,'seconds.')

    if type[1] == '+':
        identifier = '_AC'
    if type[1] == '-':
        identifier = '_CA'
    
    if generate_npy.upper() == 'Y':
        np.save(os.path.join(output_directory,type[0]+identifier+'_fock'),excitation_matrix)
    
    return excitation_matrix


lines_doc = excitdef_reader(excit_document,excitation_directory)
if generate_matrix.upper() == 'ALL':
    for type in ['H+','H-','S+','S-']:
        print(type+':')
        print(matrix(type,lines_doc,N,spin_left,spin_right,t_fock,U,mu,generate_npy))
        print()
else:
    print(type+':')
    print(matrix(type,lines_doc,N,spin_left,spin_right,t_fock,U,mu,generate_npy))
    

# Generating graph using graph.py
omega = h.hubbard(N,t_fock,U,mu,qis_not='Y')
verbose_read = 1
from graph import dvmc_spectrum
dvmc_spectrum(omega[0],verbose_read,'Y')
