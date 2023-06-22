import numpy as np
import sys

np.set_printoptions(linewidth = 1000000000,precision = 1, threshold=sys.maxsize)
A = np.load('H_CA.npy').real
a,b = A.shape
tol = 1e-13
for ii in range(a):
    for jj in range(b):
        if np.abs(A[ii,jj])<tol:
            print('     .    ',end='')
        else:
            print(' % 3.2e'%A[ii,jj],end='')
    print()


