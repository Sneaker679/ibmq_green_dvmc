import sys
import numpy as np

# Base de fock

def fock(N):
    return np.zeros((1,2*N)).astype(int)


# Creation

def cp(site,spin,fock,sign = 1):
    if type(fock) is not int: # Necessary for the hubbard function later.
        if spin == '-':
            sp = int(fock[0].size/2) # To access elements in the second half of matrix.
        if spin == '+':
            sp = 0
        if fock[0,site-1+sp] == 0: # -1 because sites start at one and arrays at 0.
            partial_bin = ''.join(str(x) for x in fock.tolist()[0][:(site-1+sp)])
            sign_exp = 0
            for ele in partial_bin:
                if ele == '1':
                    sign_exp += 1
            sign = sign*(-1)**sign_exp
            fock[0,site-1+sp] = 1
        else:
            fock = 0
    return fock,sign


# Destruction, same as creation with minor changes.

def dp(site,spin,fock,sign = 1):
    if type(fock) is not int:
        if spin == '-':
            sp = int(fock[0].size/2)
        if spin == '+':
            sp = 0
        if fock[0,site-1+sp] == 1:
            partial_bin = ''.join(str(x) for x in fock.tolist()[0][:(site-1+sp)])
            sign_exp = 0
            for ele in partial_bin:
                if ele == '1':
                    sign_exp += 1
            sign = sign*(-1)**sign_exp
            fock[0,site-1+sp] = 0
        else:
            fock = 0
    return fock,sign


# Transform binary matrix to number

def BinToNum(fock):
    fock = fock.tolist() # This line is necessary for the hubbard function, else error.
    return int(''.join(str(x) for x in fock[0]),2) # Converts string binary to number.


# Transform number to binary matrix

def NumToBin(num,N):
    binary = np.binary_repr(num,N*2)
    return np.matrix([int(i) for i in str(binary)]),1


# Pauli's matrixes

I = np.matrix([[1,0],[0,1]])

X = np.matrix([[0,1],[1,0]])

Y = np.matrix([[0,complex(0,-1)],[complex(0,1),0]])

Z = np.matrix([[1,0],[0,-1]])


    # Pauli's matrixes main operation

S = (X-complex(0,1)*Y)/2


        # To see better in the CMD

S = np.matrix([[0,0],[1,0]])


    # Pauli's creation operators

def c(N,site):
    result = 1
    for ele in range(N):
        if ele+1 < site:
            result = np.kron(result,Z)
        if ele+1 == site:
            result = np.kron(result,S)
        if ele+1 > site:
            result = np.kron(result,I)
    return result

def d(N,site):
    return c(site,N).transpose().conjugate()
    

    # Pauli's creation operators using qiskit's notation

def cq(N,site,spin):
    if spin == 'up':
        indice = 2*site - 1

    if spin == 'down':
        indice = 2*site

    result = 1
    for ele in range(2*N):
        if ele+1 < indice:
            result = np.kron(result,Z)
        if ele+1 == indice:
            result = np.kron(result,S)
        if ele+1 > indice:
            result = np.kron(result,I)
    return result

def dq(N,site,spin):
    return cq(N,site,spin).H

def nq(N,site,spin):
    return cq(N,site,spin)*dq(N,site,spin)


# Count electrons and spins

def count(state):
    count = np.count_nonzero(state == 1)
    spin = 0
    for index,value in np.ndenumerate(state):
        if value != 0:
            if index[1] < state.shape[1]/2:
                spin += 1
            else:
                spin += -1
    return count, spin


#test

'''
base,sign = NumToBin(5,2)
base,sign = dp(2,'up',base,sign)
print(base)
base,sign = cp(1,'up',base,sign)
print(base,sign)
print()

base,sign = NumToBin(9,2)
base,sign = dp(1,'up',base,sign)
print(base)
base,sign = cp(2,'up',base,sign)
print(base,sign)'''
