import fock_class as f
import numpy as np
import copy
import time
import sys

def hubbard(N=2,t=np.matrix([[0,-1],[-1,0]]),U=4,mu=2,manip='no',prt='no',qis_not='N'):

    bank = [state_num for state_num in range(4**N)]
    gs_energies = []

    if manip.lower() == 'yes':
        blocks_matrix = []
        blocks_num = []
        blocks = []
        for x in range(2*N + 1):
            blocks_matrix.append([])
            blocks_num.append([])
            blocks.append([])

    while bank:
        block = [f.Fock(N,bank[0],qis_not)]
        block_tuples = []
        for state in block:
            for i in range(N):
                for j in range(N):
                    if not i == j and not t[i,j] == 0 :
                        for spin in ['+','-']:
                            new_state = copy.deepcopy(state)
                            new_state.op('destroy',i,spin)
                            new_state.op('create',j,spin)

                            if not isinstance(new_state.fock,int):
                                block_tuples.append((state,new_state))
                                if not new_state.num in [x.num for x in block]:
                                    new_state_in_block = copy.deepcopy(new_state)
                                    new_state_in_block.sign = 1
                                    block.append(new_state_in_block)

        block_matrix = np.zeros((len(block),len(block)))

        def sorting(object):
            return object.num

        block.sort(key = sorting)        

        block_num = [x.num for x in block]
        for number in block_num:
            if number in bank:
                bank.remove(number)
        
        for tuple in block_tuples:
            hopping = tuple[0].fock - tuple[1].fock
            hopping = hopping.tolist()
            site1 = hopping[0].index(-1)
            site2 = hopping[0].index(1)

            if qis_not == 'N':
                if site1 > N - 1:
                    site1 = int(site1 - N)
                    site2 = int(site2 - N)
            else:
                if site1 % 2 == 0:
                    site1 = int(site1/2)
                else:
                    site1 = int((site1-1)/2)
                if site2 % 2 == 0:
                    site2 = int(site2/2)
                else:
                    site2 = int((site2-1)/2)
            block_matrix[block_num.index(tuple[0].num),block_num.index(tuple[1].num)] = t[site1,site2]*tuple[1].sign

        if not U == 0: 
            for index,num in enumerate(block_num):
                for site in range(N):
                    new_state = f.Fock(N,num,qis_not)
                    for spin in ['+','-']:
                        new_state.op('destroy',site,spin)
                        new_state.op('create',site,spin)
                                        
                    if not isinstance(new_state.fock,int):
                        block_matrix[index,index] += U                  
                    '''                    
                    new_state = f.Fock(N,num,qis_not)
                    new_state.op('destroy',site,'+')                    
                    new_state.op('create',site,'+')

                    if not isinstance(new_state.fock,int):
                        block_matrix[index,index] += -0.5*U

                    new_state = f.Fock(N,num,qis_not)
                    new_state.op('destroy',site,'-')
                    new_state.op('create',site,'-')

                    if not isinstance(new_state.fock,int):
                        block_matrix[index,index] += -0.5*U

                    block_matrix[index,index] += (0.5**2)*U
                    '''

        if not mu == 0:
            for index,num in enumerate(block_num):
                for site in range(N):
                    for spin in ['+','-']:
                        new_state = f.Fock(N,num,qis_not)
                        new_state.op('destroy',site,spin)
                        new_state.op('create',site,spin)
                        if not isinstance(new_state.fock,int):
                            block_matrix[index,index] += -mu
                     
        energies,numerical_states = np.linalg.eigh(block_matrix)
        gs_energy_block = min(energies)
        gs_energies.append(gs_energy_block)
        if min(gs_energies) >= gs_energy_block:
            gs_block_matrix = block_matrix
            gs_numerical_state = numerical_states[:,0].transpose()
            gs_energy_block = min(energies)
            gs_block_num = block_num
            if manip.lower() == 'yes':
                gs_block = block

        if prt.lower() == 'yes':
            print('States:',block_num)
            print(block_matrix)
            print('gs_energy_block:',gs_energy_block)
            print()
        
        if manip.lower() == 'yes':
            blocks_matrix[block[0].num_electrons].append(block_matrix)
            blocks_num[block[0].num_electrons].append(block_num)
            blocks[block[0].num_electrons].append(block)
    
    gs_energy = min(gs_energies)
    if prt.lower() == 'yes':
        print()
        print()
        print('gs_energy:',gs_energy)
        print('gs_state\n',gs_numerical_state)
        print('gs_bloc:\n',gs_block_matrix)
        print()

    if manip.lower() == 'no':
        return gs_energy,gs_numerical_state,gs_block_num,gs_block_matrix
    if manip.lower() == 'yes':
        return blocks_matrix,blocks_num,gs_block,blocks,gs_numerical_state

#print(hubbard(N=2,U=4,mu=0))

'''
Ref = -4.82842712474619
current_energy = 0
current_decimal = -.5
increment = 0.001

for i in range(1000000):
    current_decimal += increment
    tempE = hubbard(current_decimal,2,t,4,0)[0]
    print(tempE,"---",current_decimal)
    if (np.abs(tempE - Ref) < np.abs(current_energy - Ref)):
        current_energy = tempE
        print("ACCEPTED",np.abs(tempE-Ref))
    else:
        #current_decimal -= increment
        #increment/=2
        ...
'''            


