import numpy as np
import sys

# Pauli's matrixes

I = np.matrix([[1,0],[0,1]])

X = np.matrix([[0,1],[1,0]])

Y = np.matrix([[0,complex(0,-1)],[complex(0,1),0]])

Z = np.matrix([[1,0],[0,-1]])

'''
print(I)
print(X)
print(Y)
print(Z)
'''

    # Pauli's matrixes main operation

S = (X-complex(0,1)*Y)/2


        # To see better in the CMD

S = np.matrix([[0,0],[1,0]])


    # Pauli's creation operators

def c(site,N):
    result = 1
    for ele in range(N):
        if ele+1 < site:
            result = np.kron(result,Z)
        if ele+1 == site:
            result = np.kron(result,S)
        if ele+1 > site:
            result = np.kron(result,I)
    return result

def d(site,N):
    return c(site,N).transpose().conjugate()

k = np.kron

a = -0.5*k(I,k(Y,k(Z,Y))) -0.5*k(I,k(X,k(Z,X))) -2*k(I,k(I,k(I,I))) -0.5*k(Y,k(Z,k(Y,I))) -0.5*k(X,k(Z,k(X,I))) +1*k(I,k(I,k(Z,Z))) +1*k(Z,k(Z,k(I,I)))

