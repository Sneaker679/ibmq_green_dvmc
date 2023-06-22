import fock_class as f
import hubbard_classes as h
import numpy as np
import copy
import time
import sys

np.set_printoptions(linewidth = 1000,precision=4)

def excitdef_reader(document_name):

    file = open(document_name).read()
    lines_doc = file.split('\n')[5:-1]
    for line_number in range(len(lines_doc)):
        lines_doc[line_number] = lines_doc[line_number].split() 
        lines_doc[line_number] = [eval(number) for number in lines_doc[line_number]]
    
    return lines_doc


def ex_state(type,i,m,spin,gs_block_hub,gs_numerical_state,lines_doc):

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
    
    if spin == '+':
        spin_op = '-'
    if spin == '-':
        spin_op = '+'
    
    for state in excited_block:
        lines = [x for x in lines_doc if x[1] == i] 
        t = lines[m][0]
        
        ra = lines[m][2]
        rb = lines[m][3]

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

    excited_block[:] = [state for state in excited_block if not isinstance(state[0].fock,int)]
    
    return excited_block

def element(type,N,i,j,m,n,spin_left,spin_right,hubbard_output,excit_document):

    blocks_matrix = hubbard_output[0]
    blocks_num = hubbard_output[1]
    gs_block = hubbard_output[2]
    blocks = hubbard_output[3]
    gs_numerical_state = hubbard_output[4]
    
    lines_doc = excitdef_reader(excit_document)

    ex_state_left = ex_state(type,i,m,spin_left,gs_block,gs_numerical_state,lines_doc)
    ex_state_right = ex_state(type,j,n,spin_right,gs_block,gs_numerical_state,lines_doc)
    
    scalar = []
    
    if ex_state_right and ex_state_left:
        for state_left in ex_state_left:
            for state_right in ex_state_right:

                if type[0] == 'S':
                    if np.array_equal(state_left[0].fock,state_right[0].fock):
                        value = state_left[1]*state_left[0].sign*state_right[1]*state_right[0].sign
                        scalar.append(value)

                if type[0] == 'H':
                    if (state_left[0].num_electrons == state_right[0].num_electrons
                    and state_left[0].total_spin == state_right[0].total_spin):

                        n_electrons = state_left[0].num_electrons
                        for block_num in blocks_num[n_electrons]:
                            if ex_state_left[0][0].num in block_num:
                                index_sp = blocks_num[n_electrons].index(block_num)

                        row = blocks_num[n_electrons][index_sp].index(state_left[0].num)
                        column = blocks_num[n_electrons][index_sp].index(state_right[0].num)

                        value = blocks_matrix[n_electrons][index_sp][row,column]*state_left[1]*state_left[0].sign*state_right[1]*state_right[0].sign
                        scalar.append(value)
        sum_values = sum(scalar)
        return sum_values
    else:
        return 0


def matrix(type,excit_document,N,spin_left,spin_right,t,U,mu,save='Y'):
    
    # This is the total possible number of excitation for each site, as explained by the paper.
    lines_doc = excitdef_reader(excit_document)
    lines = [values for values in lines_doc if values[1] == 0]
    N_exc = len(lines)
    start = time.time()

    matrix_size = N * N_exc
    excitation_matrix = np.zeros((matrix_size,matrix_size))
    
    hubbard_output = h.hubbard(N,t,U,mu,'yes')
    
    for i in range(N):
        for j in range(N):
            for m in range(N_exc):
                for n in range(N_exc): 
                    column_num = N * m + i
                    row_num = N * n + j
                    if column_num - row_num >= 0:
                        excitation_matrix[row_num,column_num] = element(type,N,i,j,m,n,spin_left,spin_right,hubbard_output,excit_document)
    
    excitation_matrix = np.tril(excitation_matrix.T,-1) + excitation_matrix    

    end = time.time()
    print('Time:',end-start,'seconds.')

    if save.upper() == 'Y':
        if type[1] == '+':
            identifier = '_AC'
        if type[1] == '-':
            identifier = '_CA'
        np.save(type[0]+identifier,excitation_matrix)
    
    return excitation_matrix

if __name__ == "__main__":
    # Which matrix? (H+, H-, S+ or S-.)
    type = 'ALL'

    # Number of sites in total and number of neighbors the site with the least neighbor has??
    N = 2

    # Spins
    spin_left = '+'
    spin_right = '+'

    # Values
     
    t = np.matrix([
        [0,1],
        [1,0]
    ])
    
    '''
    t = np.matrix([
        [0,1,1,0],
        [1,0,0,1],
        [1,0,0,1],
        [0,1,1,0]
    ])
    '''
    '''
    t = np.matrix([
        [0,1,0,1,0,0],
        [1,0,1,0,1,0],
        [0,1,0,0,0,1],
        [1,0,0,0,1,0],
        [0,1,0,1,0,1],
        [0,0,1,0,1,0]
    ])
    '''

    U = 4
    mu = 3
    

    # Excitation document name
    excit_doc = f'excitation{N}sites.def'
    generate_npy = 'Y' # 'Y' or 'N'.

    if type.upper() == 'ALL':
        for type_ in ['H+','H-','S+','S-']:
            print(type_+':')
            print(matrix(type_,excit_doc,N,spin_left,spin_right,t,U,mu,generate_npy))
            print()
    else:
        print(type+':')
        print(matrix(type,excit_doc,N,spin_left,spin_right,t,U,mu,generate_npy))
