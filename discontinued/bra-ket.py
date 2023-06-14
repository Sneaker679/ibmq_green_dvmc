import numpy as np

# Complex numbers et their conjugate

com1 = complex(1,2)

com1.conjugate()

# Matrix multiplication

matrix1 = np.matrix([[1,3],[3,1]])
matrix2 = np.matrix([[2,0],[1,2]])

result = np.dot(matrix1,matrix2)


# Scalar product

matrix1 = np.array([1,2])
matrix2 = np.array([1,3])

result = np.dot(matrix1,matrix2)


# Conjugate of matrix

matrix1 = np.array([1,complex(3,4)])
matrix1.conjugate()


# ket to bra

ket = np.matrix([[1],[2],[complex(1,4)]])
bra = ket.transpose().conjugate()


# bra to ket

bra = np.matrix([1,2,3])
ket = bra.transpose()


# Tensor product (Kronecker product)

matrix1 = np.matrix([[1,0],[0,1]])
matrix2 = np.matrix([[complex(1,2),1],[1,complex(3,4)]])

result = np.kron(matrix1,matrix2)
