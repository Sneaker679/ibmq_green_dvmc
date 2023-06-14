import functions as f
import hubbard_simplified as h
import numpy as np
import copy
import time

np.set_printoptions(linewidth = 1000)

def Element_H(type,N,i,j,m,n,spin_left,spin_right,t,U,mu):
    GS,E,bloc = h.hubbard(N,t,U,mu)
    blocs_matrix,blocs_num,GS_bloc_bin_init = h.hubbard(N,t,U,mu,'no','yes')

    for index in range(len(GS_bloc_bin_init)):
        GS_bloc_bin_init[index] = [GS_bloc_bin_init[index],GS[index,0],1]

    if type != 'H-' and type != 'H+':
        raise Exception('The "type" parameter has to be "H-" or "H+".')
    if m < 0 or n < 0:
        raise Exception('Parameters "m" and "n" have to be positive integers.')
    if (spin_right != '+' and spin_right != '-') or (spin_left != '-' and spin_left != '+'):
        raise Exception('Parameters "spin_left" and "spin_right" must be equal to "-" or "+" (both being strings).')

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

    GS_bloc_bin = copy.deepcopy(GS_bloc_bin_init)
    if m == 0:
        if type == 'H+':
            for index in range(len(GS_bloc_bin)):
                GS_bloc_bin[index][0],sign = f.cp(i+1,spin_left,GS_bloc_bin[index][0])
                GS_bloc_bin[index][2] = sign
        if type == 'H-':
            for index in range(len(GS_bloc_bin)):
                GS_bloc_bin[index][0],sign = f.dp(i+1,spin_left,GS_bloc_bin[index][0])
                GS_bloc_bin[index][2] = sign
    elif m == 1:
        if type == 'H+':            
            for index in range(len(GS_bloc_bin)):
                GS_bloc_bin[index][0],sign = f.dp(i+1,spin_left_op,GS_bloc_bin[index][0])
                GS_bloc_bin[index][0],sign = f.cp(i+1,spin_left_op,GS_bloc_bin[index][0],sign)
                GS_bloc_bin[index][0],sign = f.cp(i+1,spin_left,GS_bloc_bin[index][0],sign)
                GS_bloc_bin[index][2] = sign
        if type == 'H-':
            for index in range(len(GS_bloc_bin)):
                GS_bloc_bin[index][0],sign = f.dp(i+1,spin_left_op,GS_bloc_bin[index][0])
                GS_bloc_bin[index][0],sign = f.cp(i+1,spin_left_op,GS_bloc_bin[index][0],sign)
                GS_bloc_bin[index][0],sign = f.dp(i+1,spin_left,GS_bloc_bin[index][0],sign)
                GS_bloc_bin[index][2] = sign
    else:
        lines_left = [ex_state for ex_state in lines if ex_state[1] == i]
        
        ra = lines_left[m][2]
        rb = lines_left[m][3] 
         
        if type == 'H+':
            for index in range(len(GS_bloc_bin)):
                GS_bloc_bin[index][0],sign = f.dp(rb+1,spin_left,GS_bloc_bin[index][0])
                GS_bloc_bin[index][0],sign = f.cp(rb+1,spin_left,GS_bloc_bin[index][0],sign)
                GS_bloc_bin[index][0],sign = f.dp(ra+1,spin_left_op,GS_bloc_bin[index][0],sign)
                GS_bloc_bin[index][0],sign = f.cp(ra+1,spin_left_op,GS_bloc_bin[index][0],sign)
                GS_bloc_bin[index][0],sign = f.cp(i+1,spin_left,GS_bloc_bin[index][0],sign)
                GS_bloc_bin[index][2] = sign
        if type == 'H-':            
            for index in range(len(GS_bloc_bin)):
                GS_bloc_bin[index][0],sign = f.dp(ra+1,spin_left_op,GS_bloc_bin[index][0])
                GS_bloc_bin[index][0],sign = f.cp(ra+1,spin_left_op,GS_bloc_bin[index][0],sign)
                GS_bloc_bin[index][0],sign = f.dp(rb+1,spin_left,GS_bloc_bin[index][0],sign)
                GS_bloc_bin[index][0],sign = f.cp(rb+1,spin_left,GS_bloc_bin[index][0],sign)
                GS_bloc_bin[index][0],sign = f.dp(i+1,spin_left,GS_bloc_bin[index][0],sign)
                GS_bloc_bin[index][2] = sign
   
    GS_bloc_bin[:] = [state for state in GS_bloc_bin if not isinstance(state[0],int)]
    for index in range(len(GS_bloc_bin)):
        GS_bloc_bin[index][0] = f.BinToNum(GS_bloc_bin[index][0])
    ex_state_left = copy.deepcopy(GS_bloc_bin)
    

    GS_bloc_bin = copy.deepcopy(GS_bloc_bin_init)
    if n == 0:
        if type == 'H+':
            for index in range(len(GS_bloc_bin)):
                GS_bloc_bin[index][0],sign = f.cp(j+1,spin_right,GS_bloc_bin[index][0])
                GS_bloc_bin[index][2] = sign
        if type == 'H-':
            for index in range(len(GS_bloc_bin)):
                GS_bloc_bin[index][0],sign = f.dp(j+1,spin_right,GS_bloc_bin[index][0])
                GS_bloc_bin[index][2] = sign
    elif n == 1:
        if type == 'H+':            
            for index in range(len(GS_bloc_bin)):
                GS_bloc_bin[index][0],sign = f.dp(j+1,spin_right_op,GS_bloc_bin[index][0])
                GS_bloc_bin[index][0],sign = f.cp(j+1,spin_right_op,GS_bloc_bin[index][0],sign)
                GS_bloc_bin[index][0],sign = f.cp(j+1,spin_right,GS_bloc_bin[index][0],sign)
                GS_bloc_bin[index][2] = sign
        if type == 'H-':
            for index in range(len(GS_bloc_bin)):
                GS_bloc_bin[index][0],sign = f.dp(j+1,spin_right_op,GS_bloc_bin[index][0])
                GS_bloc_bin[index][0],sign = f.cp(j+1,spin_right_op,GS_bloc_bin[index][0],sign)
                GS_bloc_bin[index][0],sign = f.dp(j+1,spin_right,GS_bloc_bin[index][0],sign)
                GS_bloc_bin[index][2] = sign
    else:
        lines_right = [ex_state for ex_state in lines if ex_state[1] == j]
        
        ra = lines_right[n][2]
        rb = lines_right[n][3] 
         
        if type == 'H+':
            for index in range(len(GS_bloc_bin)):
                GS_bloc_bin[index][0],sign = f.dp(rb+1,spin_right,GS_bloc_bin[index][0])
                GS_bloc_bin[index][0],sign = f.cp(rb+1,spin_right,GS_bloc_bin[index][0],sign)
                GS_bloc_bin[index][0],sign = f.dp(ra+1,spin_right_op,GS_bloc_bin[index][0],sign)
                GS_bloc_bin[index][0],sign = f.cp(ra+1,spin_right_op,GS_bloc_bin[index][0],sign)
                GS_bloc_bin[index][0],sign = f.cp(j+1,spin_right,GS_bloc_bin[index][0],sign)
                GS_bloc_bin[index][2] = sign
        if type == 'H-':            
            for index in range(len(GS_bloc_bin)): 
                GS_bloc_bin[index][0],sign = f.dp(ra+1,spin_right_op,GS_bloc_bin[index][0])
                GS_bloc_bin[index][0],sign = f.cp(ra+1,spin_right_op,GS_bloc_bin[index][0],sign)
                GS_bloc_bin[index][0],sign = f.dp(rb+1,spin_right,GS_bloc_bin[index][0],sign)
                GS_bloc_bin[index][0],sign = f.cp(rb+1,spin_right,GS_bloc_bin[index][0],sign)
                GS_bloc_bin[index][0],sign = f.dp(j+1,spin_right,GS_bloc_bin[index][0],sign) 
                GS_bloc_bin[index][2] = sign

    GS_bloc_bin[:] = [state for state in GS_bloc_bin if not isinstance(state[0],int)]
    for index in range(len(GS_bloc_bin)):
        GS_bloc_bin[index][0] = f.BinToNum(GS_bloc_bin[index][0]) 
    ex_state_right = copy.deepcopy(GS_bloc_bin)
    del GS_bloc_bin

    scalar = []
    if ex_state_left and ex_state_right:
        for state_left in ex_state_left:
            for state_right in ex_state_right:

                num_left,sp_left = f.count(f.NumToBin(ex_state_left[0][0],N)[0])
                num_right,sp_right = f.count(f.NumToBin(ex_state_right[0][0],N)[0])
                
                if num_left == num_right and sp_left == sp_right:
                    n = num_left
                    for bloc_num in blocs_num[num_left]:
                        if ex_state_left[0][0] in bloc_num:
                            index_sp = blocs_num[num_left].index(bloc_num)
                    
                    blocs_num[n][index_sp]
                    blocs_matrix[n][index_sp]
                    
                    row = blocs_num[n][index_sp].index(state_left[0])
                    column = blocs_num[n][index_sp].index(state_right[0])
                    
                    value = blocs_matrix[n][index_sp][row,column] * state_left[1] * state_right[1] * state_left[2] * state_right[2]
                    scalar.append(value)
        sum_values = sum(scalar)
        return sum_values
    else:
        return 0

def matrix(type,N,N_min,spin_left,spin_right,t,U,mu):

    start = time.time()
    
    N_exc = 2 + N_min*(N_min+1)

    matrix_size = N * N_exc
    excitation_matrix = np.zeros((matrix_size,matrix_size))
    
    for i in range(N):
        for j in range(N):
            for m in range(N_exc):
                for n in range(N_exc): 
                    column_num = N_exc * i + m
                    row_num = N_exc * j + n
                    if column_num - row_num >= 0:
                        excitation_matrix[row_num,column_num] = Element_H(type,N,i,j,m,n,spin_left,spin_right,t,U,mu)

    excitation_matrix = np.tril(excitation_matrix.T,-1) + excitation_matrix

    end = time.time()
    print('Time:',end - start,'seconds.')

    return excitation_matrix

#matrix(S+,2,1,+,+,1,4,2)

type = input('Which matrix? (H+ or H-.)  ')
N = int(input('Number of sites in total?  '))
N_min = int(input('Number of neighbors the site with the least neighbor has?  '))
spin_left = input('Spin of left side?  ')
spin_right = input('Spin of right side  ')
t = int(input('Value of t?  '))
U = int(input('Value of U?  '))
mu = int(input('Value of mu?  '))

print(matrix(type,N,N_min,spin_left,spin_right,t,U,mu))

