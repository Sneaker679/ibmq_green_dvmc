import numpy as np
import sys

np.set_printoptions(threshold=sys.maxsize)


# Base de fock

def fock(N):
    return np.zeros((1,2*N)).astype(int)


# Creation

def cp(site,spin,fock):
    if type(fock) is not int: # Necessary for the hubbard function later.
        if spin == 'down':
            sp = int(fock[0].size/2) # To access elements in the second half of matrix.
        if spin == 'up':
            sp = 0
        if fock[0,site-1+sp] == 0: # -1 because sites start at one and arrays at 0.
            fock[0,site-1+sp] = 1
        else:
            fock = 0
    return fock


# Destruction, same as creation with minor changes.

def dp(site,spin,fock):
    if type(fock) is not int:
        if spin == 'down':
            sp = int(fock[0].size/2)
        if spin == 'up':
            sp = 0
        if fock[0,site-1+sp] == 1:
            fock[0,site-1+sp] = 0
        else:
            fock = 0
    return fock


# Transform binary matrix to number

def BinToNum(fock):
    fock = fock.tolist() # This line is necessary for the hubbard function, else error.
    return int(''.join(str(x) for x in fock[0]),2) # Converts string binary to number.


# Transform number to binary matrix

def NumToBin(num,N):
    bin = np.binary_repr(num,N*2)
    return np.matrix([int(i) for i in str(bin)])


# Hubbard hamiltonien calculation

def hubbard(N):
    bank = [] # Creation of empty bank of states.
    for ele in range(4**N): # Fill bank.
        bank.append(ele)
    hubbard = np.zeros((4**N,4**N)).astype(int).astype(str) # Create empty matrix of strings '0'.
    
    # Calculation of potential energy of ground state for each state in the bank.
    for state in bank: 
        '''The number of array this list contains will be equal to the number of 'U' energy at the state. 
        That's because the result will always be the same for the sum after U, either 0 or 1(times U). All the results
        are stored in this list, hence the lenght of the list is equal to the number preceding U.'''
        nU = []
        for site in range(N):
            u = NumToBin(state,N) # Creation of binary number of the state used for the calculation.
            for spin in ['up','down']: # One calculation per spin.
                u = dp(site,spin,u) # We first apply the annihilation per dirac's notation.
                u = cp(site,spin,u) # Then the creation for the the same particle at the same site.
                # The possible scenarios are: 
                ''' 1. the annihilation destroys the state and sum = 0.
                    2. the annihilation destroys the particle of the state, then the creation recreates that particle,
                    leaving the final state unchanged from the initial one.'''
            if not isinstance(u,int): # If the result is not 0 (state is unchanged):
                nU.append(u) # We append the result in nU.
                '''After the for loop has ended, nU will contain all the unchanged states. 
                Hence, all the elements of the list are the same.'''
        if nU: # If nU is not empty.
            for ele in range(len(nU)): # We assign (the lenght of nU) + 'U' to the proper state in the hamiltonian.  
                if len(nU) == 1:
                    hubbard[BinToNum(nU[ele]),state] = 'U'
                else: 
                    hubbard[BinToNum(nU[ele]),state] = str(len(nU)) + 'U'      
    
    while bank: # While bank is not empty.
        bloc = [] # Creation of empty bloc list.
        bloc.append(NumToBin(bank[0],N)) # We add to bloc the first element of the bank. 
        h = 0 # This is useful later. Serves to keep track of where we are in bloc.
        
        """For each element in bloc, we need to calculate the states and ensure the bloc
        doesn't get any bigger. bloc is updated at the end of one of these loops. This prevents 
        unnecessary calculations."""
        for ele in bloc:
            # The following lines merely represent the sum for the first term of Hubbard's model.
            for spin in ['up','down']:
                for i in range(N):
                    for j in range(N):
                        if i != j:

                            '''Each loop, bloc is appended the state (or ele) for which
                            we are doing the calculation. In other word, the state at 
                            position h is duplicated at the end of bloc. We then operate
                            with cp and dp on it and check the result.'''

                            bloc.append(bloc[h]+fock(N))
                            bloc[-1] = dp(i+1,spin,bloc[-1])
                            bloc[-1] = cp(j+1,spin,bloc[-1])
                            
                            if isinstance(bloc[-1],int): # If the result is 0, delete from bloc.
                                del bloc[-1]
                                
                                '''If the result is not 0, that means the state has an energy.
                                Hence the "else":'''
                            else:
                                '''The following line adds the energy -t at the appropriate
                                place in the hubbard matrix. It is possible that in a previous
                                loop, said energy was already added. This means that this new
                                result might be a duplicate, meaning the result is already in bloc.
                                For optimisation purposes, we remove the duplicate from bloc.'''
                                hubbard[BinToNum(bloc[-1]),BinToNum(bloc[h])] = 't' 
                                for state in bloc[:(len(bloc)-1)]:
                                    if np.array_equal(bloc[-1],state):
                                        del bloc[-1]
                                        break
            """If we reach this line, it means that, if bloc contains any more unique state,
            we need to duplicate the >next< element. This next element is controlled trough the
            varaible 'h'. At the end of the loop, we simply add one. At the next loop, bloc will
            be appended at its end, the h element of bloc, on which we are doing the calculations."""
            h += 1
        '''This line removes from the bank the states that have been calculated.'''
        for ele in range(len(bloc)):
            bloc[ele] = BinToNum(bloc[ele])
            bank.remove(bloc[ele])
        
    print(hubbard) 
    
hubbard(3)
