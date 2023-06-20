import fock_class as f
import hubbard_classes as h
import numpy as np
import copy
import time
import sys

np.set_printoptions(linewidth = 10000
        #,threshold=sys.maxsize
        )

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
    if type[1] == '-':
        c_operator = 'destroy'
    
    if spin == '+':
        spin_op = '-'
    if spin == '-':
        spin_op = '+'
    
    for state in excited_block:
        lines = [x for x in lines_doc if x[1] == i] 
        ra = lines[m][2]
        rb = lines[m][3]

        if m == 0:
            state[0].op(c_operator,i,spin)
        elif m == 1:
            state[0].op('destroy',ra,spin_op)
            state[0].op('create',ra,spin_op)
            state[0].op(c_operator,i,spin)
        else:
            state[0].op('destroy',ra,spin_op)
            state[0].op('create',ra,spin_op)
            state[0].op('destroy',rb,spin)
            state[0].op('create',rb,spin)
            state[0].op(c_operator,i,spin)

    excited_block[:] = [state for state in excited_block if not isinstance(state[0].fock,int)]
    
    return excited_block

def element(type,N,i,j,m,n,spin_left,spin_right,hubbard_output,excit_document):

    blocks_matrix = hubbard_output[0]
    blocks_num = hubbard_output[1]
    gs_block = hubbard_output[2]
    blocks = hubbard_output[3]
    gs_numerical_state = hubbard_output[4]
   
    if i==3 and j==3 and m==13 and n==13:
        print('gs',[x.num for x in gs_block])
        print()
    
    lines_doc = excitdef_reader(excit_document)

    ex_state_left = ex_state(type,i,m,spin_left,gs_block,gs_numerical_state,lines_doc)
    ex_state_right = ex_state(type,j,n,spin_right,gs_block,gs_numerical_state,lines_doc)
    
    scalar = []
    
    if i==3 and j==3 and m==13 and n==13:
        print('left',[x[0].num for x in ex_state_left])
        print()
        print('right',[x[0].num for x in ex_state_right])
        print()

        k=0
        for index1,ele1 in enumerate(blocks_matrix):
            print('N=',k)
            for index2,ele2 in enumerate(ele1):
                print([x.num for x in blocks[index1][index2]])
                print([x.fock for x in blocks[index1][index2]])
                print(ele2)
            print()
            k += 1
    

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

def matrix(type,excit_document,N,N_min,spin_left,spin_right,t,U,mu,save='Y'):

    N_exc = 2 + N_min*(N_min+1)

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

    if save.upper() == 'Y':
        if type[1] == '+':
            identifier = '_AC'
        if type[1] == '-':
            identifier = '_CA'
        np.save(type[0]+identifier,excitation_matrix)
    
    return excitation_matrix


'''gs_energy,gs_numerical_state,gs_block_matrix = h.hubbard(N,t,U,mu)
blocks_matrix,blocks_num,gs_block,blocks = h.hubbard(N,t,U,mu,'yes')'''


if __name__ == "__main__":
    # Which matrix? (H+, H-, S+ or S-.)
    type = 'H+'

    # Number of sites in total and number of neighbors the site with the least neighbor has??
    N = 4
    N_min = 3

    # Spins
    spin_left = '+'
    spin_right = '+'

    # Values
    '''
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
    U = 4
    mu = 2

    # Excitation document name
    excit_doc = f'excitation{N}sites.def'

    print(type+':')
    print(matrix(type,excit_doc,N,N_min,spin_left,spin_right,t,U,mu))
