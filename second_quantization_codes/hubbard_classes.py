from fock_class import Fock
from math import isclose
import numpy as np
import copy
import time
import sys

"""This will probably be difficult to understand... :("""

"""The goal of this function is to calculate everything related to the hamiltonian of the 
Fermi-Hubbard Model in Fock's basis. This includes:

- The ground state energy
- The ground state vector
- The ground state block (state vectors that form the ground state in their binary number form)
- The ground state block matrix

Or when manip = 'Y':
- All the block matrices (named blocks_matrix here for simplicity).
- All the blocks (state vectors in their binary form).
- The ground state block (state vectors that form the ground state in their binary number form)
- All the blocks (in their Fock object form. Said objects are defined in fock_class.py)
- The ground state vector
"""

def hubbard(N=2,t=-1,hopping_matrix=np.matrix([[0,1],[1,0]]),U=4,mu=2,spin_gs='+',manip=False,prt=False,qis_not=False):
    """Parameters
    N: Number of sites
    t: Hopping matrix
    U: Potential energy
    mu: Chemical potential
    spin_gs: The total spin the ground state should have ('+' or '-')
    manip: Boolean parameter that dictates if the function should return way more things (list is above). This can be heavy on the memory.
    prt: Boolean parameter that dictates if the function should print some of the information in the CMD.
    qis_not: Boolean parameter that dictates if the calculation should be made using qiskit's Fock basis notation (|0up 0down 1up 1down ...>).
    """

    # The bank is a list of states that are yet to be operated on to determine which block is belongs to.
    bank = [state_num for state_num in range(4**N)]

    # Defined for later
    gs_energy = 1e5
    gs_energies = []
    gs_numerical_states = []
    gs_states = []
    gs_blocks = []

    if manip is True:
        blocks_matrix = []
        blocks_num = []
        blocks = []
        gs_blocks_matrix = []
        for x in range(2*N + 1):
            blocks_matrix.append([])
            blocks_num.append([])
            blocks.append([])

    # While the bank is not empty, we want to find all the blocks.
    while bank:
        # Initializing Block list with in it its first state we are operating on.
        block = [Fock(N,bank[0],qis_not)]
        # This will store all the hoppings made, which will determine the position of the t value in the hamitlonian,
        # as well as the phase.
        block_tuples = []

        for state in block:
            for i in range(N):
                for j in range(N):
                    if not i == j and not hopping_matrix[i,j] == 0 :
                        for spin in ['+','-']:
                            """We apply creation/annihilation operators on initial state to calculate 
                            the other states in the same block."""
                            new_state = copy.deepcopy(state)
                            new_state.op('destroy',i,spin)
                            new_state.op('create',j,spin)

                            """The Fock class is supposed to turn the numpy matrix (state) into an integer of 0
                            if the state is completely destroyed by either the creation of a particle that is 
                            already there or the annihilation of a void. Thus, if the state was NOT destroyed, it
                            counts as a hopping. Thus the following lines are executed to store the hopping in
                            block_tuples."""
                            if not isinstance(new_state.fock,int):
                                block_tuples.append((state,new_state))
                                if not new_state.num in [x.num for x in block]:
                                    # These lines make sure that the state added to the block list has a phase '+'.
                                    new_state_in_block = copy.deepcopy(new_state)
                                    new_state_in_block.sign = 1
                                    block.append(new_state_in_block)

        # Creation of the empty matrix depending on the size of the block
        block_matrix = np.zeros((len(block),len(block)))

        # Defining a rule to sort the block list from smallest to biggest binary number
        def num_sorting(object):
            return object.num
        block.sort(key = num_sorting)        

        # We duplicate the block list to have a more readable form in which the states are represented by their binary number
        block_num = [x.num for x in block]
        for number in block_num:
            if number in bank:
                # We empty the bank of the states already calculated
                bank.remove(number)
        
        # Adding t values into the hamiltonian block. Extra lines compensate for the use or not of qiskit's notation.
        for tuple in block_tuples:
            hopping = tuple[0].fock - tuple[1].fock
            hopping = hopping.tolist()
            site1 = hopping[0].index(-1)
            site2 = hopping[0].index(1)

            if qis_not is False:
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
            block_matrix[block_num.index(tuple[0].num),block_num.index(tuple[1].num)] = t*hopping_matrix[site1,site2]*tuple[1].sign

        # We add the potential term in the hamiltonian block
        if not U == 0: 
            for index,num in enumerate(block_num):
                for site in range(N):
                    new_state = Fock(N,num,qis_not)
                    for spin in ['+','-']:
                        new_state.op('destroy',site,spin)
                        new_state.op('create',site,spin)
                                        
                    if not isinstance(new_state.fock,int):
                        block_matrix[index,index] += U                  

                    '''
                    #This is to symetrize the hamiltonian. Don't bother.
                    new_state = Fock(N,num,qis_not)
                    new_state.op('destroy',site,'+')                    
                    new_state.op('create',site,'+')

                    if not isinstance(new_state.fock,int):
                        block_matrix[index,index] += -0.5*U

                    new_state = Fock(N,num,qis_not)
                    new_state.op('destroy',site,'-')
                    new_state.op('create',site,'-')

                    if not isinstance(new_state.fock,int):
                        block_matrix[index,index] += -0.5*U

                    block_matrix[index,index] += (0.5**2)*U
                    '''
        # We add the chemical potential to the hamitlonian block
        if not mu == 0:
            for index,num in enumerate(block_num):
                for site in range(N):
                    for spin in ['+','-']:
                        new_state = Fock(N,num,qis_not)
                        new_state.op('destroy',site,spin)
                        new_state.op('create',site,spin)
                        if not isinstance(new_state.fock,int):
                            block_matrix[index,index] += -mu

        # We calculate the eigen value of the now complete block matrix
        energies,numerical_states = np.linalg.eigh(block_matrix)
        gs_energy_block = min(energies)

        # This is to be able to choose which groundstate we want. The positive spin one or negative spin one.
        block_spin = block[0].total_spin

        # If multiple blocks have the same fundamental energy, but not the same ground state, this variable is declared
        # to make sure that the ground state is chosen from the biggest block.
        block_len = len(block)
        
        # We only consider the blocks that fit the ground state we are after. Positive total spin or Negative total spin.
        if (spin_gs == '+' and block_spin >= 0) or (spin_gs == '-' and block_spin <= 0):
            if gs_energy >= gs_energy_block:
                if gs_energy > gs_energy_block:
                    gs_blocks_matrix = []
                    gs_blocks_num = []
                    gs_numerical_states = []
                    gs_blocks = []
                gs_blocks_matrix.append(block_matrix)
                gs_blocks_num.append(block_num)
                gs_numerical_states.append(numerical_states[:,0].transpose())
                gs_energy = gs_energy_block
                if manip is True:
                    gs_blocks.append(block)

        if prt is True:
            print('States:',block_num)
            print(block_matrix)
            print('gs_energy_block:',gs_energy_block)
            print()
        
        if manip is True:
            blocks_matrix[block[0].num_electrons].append(block_matrix)
            blocks[block[0].num_electrons].append(block)
    
    # Rule to sort the blocks depending on their total spins
    def spin_sorting(object_list):
        return object_list[0][0].total_spin
   
    if manip is True:
        for index,blocks_ele in enumerate(blocks):
            blocks_matrix_ele = blocks_matrix[index]
            blocks[index] = [x for x, y in sorted(zip(blocks_ele, blocks_matrix_ele), key=spin_sorting)]
            blocks_matrix[index] = [y for x, y in sorted(zip(blocks_ele, blocks_matrix_ele), key=spin_sorting)]

        blocks_num = copy.deepcopy(blocks)
        for index1,block_group in enumerate(blocks_num):
            for index2,block in enumerate(block_group):
                for index3,fock_object in enumerate(block):
                    blocks_num[index1][index2][index3] = blocks_num[index1][index2][index3].num

    if prt is True:
        print()
        print()
        print('gs_energy:',gs_energy)
        print('gs_states:')
        for state in gs_numerical_states:
            print(state)
        print('gs_blocks_matrix:')
        for matrix in gs_blocks_matrix:
            print(matrix)
        print()

    if manip is False:
        return gs_energy,gs_numerical_states,gs_blocks_num,gs_blocks_matrix
    if manip is True:
        return blocks_matrix,blocks_num,gs_blocks,blocks,gs_numerical_states
