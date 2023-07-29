from hubbard_classes import hubbard
import numpy as np

num_columns = int(input('Number of columns: '))
num_rows = int(input('Number of rows: '))
N = num_rows * num_columns

### Automatic hopping matrix
# Generating lattice's site coordinates
lattice_fock = []
for row in range(num_rows):
    for column in range(num_columns):
        lattice_fock.append(np.matrix([row,column]))

# Generating hopping matrix using the coordinates defined before
hopping_matrix = np.zeros((N,N))
for index1,site1 in enumerate(lattice_fock):
    for index2,site2 in enumerate(lattice_fock):
        
        '''The condition for checking if a hopping is possible is here. We are substracting the coordinates
        and checking if the result is a possible hopping in a grid lattice. If it is, the appropriate term
        is added to the hopping matrix.'''
        hop = np.subtract(site1,site2)
        if (np.array_equal(hop,np.matrix([1,0]))
            or np.array_equal(hop,np.matrix([-1,0]))
            or np.array_equal(hop,np.matrix([0,1]))
            or np.array_equal(hop,np.matrix([0,-1]))):

            hopping_matrix[index1,index2] = 1

possible_mu = [-4,-3,-2,-1,0,1,2,3,4,5,6,7,8]

for mu in possible_mu:
    for spin_gs in ['+','-']:
        print()
        print(f'For mu={mu} and spin_gs={spin_gs}:')
        output = hubbard(N=N,spin_gs=spin_gs,mu=mu,hopping_matrix=hopping_matrix,manip=True)
        #print(hubbard(N=N,spin_gs=spin_gs,mu=mu,hopping_matrix=hopping_matrix,manip=False,qis_not=True))
        print("GS")
        for blocks in output[2]:
            print(' [',end='')
            for block in blocks:
                print(f'{block.num} ',end='')
            print('] ',end='')
        print()

        for blocks in output[3]:
            for block in blocks:
                size = len(block)
                Sz = block[0].total_spin
                N_ = block[0].num_electrons
                print()
                print(f'Size={size}.   (N,Sz)({N_},{Sz})   ',end='')
                if block in output[2]:
                    print("GS",end='')
        print()
        print()
