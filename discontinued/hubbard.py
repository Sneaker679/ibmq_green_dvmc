import functions as f
import numpy as np
import copy
import time

def hubbard(N,t=1,U=4,mu=2,prt='no',manip='no',time='no'):
    
    if time == 'yes':
        start = time.time()
    
    # Creation of the bank of state
    bank = [state_num for state_num in range(4**N)]
    ground_states_energies = []       
    
    # For later
    if manip == 'yes':
        blocs_matrix = []
        blocs_num = []
        for x in range(2*N + 1):
            blocs_matrix.append([])
            blocs_num.append([])

    # Each iteration of this while loop deletes items from bank until all states
    # have been calculated.
    while bank:
        
        # States of the bloc in binary. We add the first element of bank in it.
        bloc_bin = [f.NumToBin(bank[0],N)[0]]
        
        # We initialize this list for later. States of the bloc representated
        # with numbers instead of binary.
        bloc_num = []
	
        # We store in this list all the associations between states.
        bloc_tuples = []
        
        # The sum that creates the full bloc and also all the associations.
        for state in bloc_bin:
            for i in range(N):
                for j in range(N):
                    if j != i:
                        for spin in ['+','-']:
                            current_state = copy.deepcopy(state)
                            current_state,sign = f.dp(i+1,spin,current_state)
                            current_state,sign = f.cp(j+1,spin,current_state,sign)
                            """print(state)
                            print(current_state)
                            print()"""
                            
                            if not isinstance(current_state,int):
                                bloc_bin.append(current_state)
                                for ele in bloc_bin[:(len(bloc_bin)-1)]:
                                        if np.array_equal(current_state,ele):
                                            del bloc_bin[-1]
                                bloc_tuples.append((f.BinToNum(state),
                                                    f.BinToNum(current_state),sign))        
        
        # Creating a new bloc list that lists the states by number instead.
        for bin_state in bloc_bin:
            bin_state = f.BinToNum(bin_state)
            bloc_num.append(bin_state)
            bank.remove(bin_state)
        bloc_num.sort()
        
        # Redefining the tuples to match the bloc_matrix's available indexes.
        bloc_newtuples = []
        for tuple_ in bloc_tuples:
            new_tuple1 = bloc_num.index(tuple_[0])
            new_tuple2 = bloc_num.index(tuple_[1])
            bloc_newtuples.append((new_tuple1,new_tuple2,tuple_[2]))
        del bloc_tuples
        
        # Creation of the empty bloc matrix.
        lenght_bloc = len(bloc_bin)
        bloc_matrix = np.zeros((lenght_bloc,lenght_bloc))
        
        # Adding t
        for ele in bloc_newtuples:
            bloc_matrix[ele[:2]] = -t*ele[2]
        
        # Adding U
        if U != 0:
            for state in bloc_num:
                state_index = bloc_num.index(state)
                for site in range(N):
                    current_state,sign = f.NumToBin(state,N)
                    for spin in ['+','-']:
                        current_state,sign = f.dp(site+1,spin,current_state,sign)
                        current_state,sign = f.cp(site+1,spin,current_state,sign)
                        
                    if not isinstance(current_state,int):
                        bloc_matrix[state_index,state_index] += U

        # Adding mu
        if mu != 0:
            for state in bloc_num:
                state_index = bloc_num.index(state)
                for site in range(N):
                    for spin in ['+','-']:
                        current_state,sign = f.NumToBin(state,N)
                        current_state,sign = f.dp(site+1,spin,current_state,sign)
                        current_state,sign = f.cp(site+1,spin,current_state,sign)
                        if not isinstance(current_state,int):
                            bloc_matrix[state_index,state_index] += -mu

        # Diagonalisation exacte
        E,S = np.linalg.eigh(bloc_matrix)
        ground_states_energies.append(min(E))
        
        if manip == 'yes':
            n,sp = f.count(bloc_bin[0])
            blocs_matrix[n].append(bloc_matrix)
            blocs_num[n].append(bloc_num)


        if min(ground_states_energies) >= min(E):
            GS_bloc_matrix = bloc_matrix
            #GS_bloc_num = bloc_num
            #GS_bloc_bin = bloc_bin
            if manip == 'yes':
                GS_bloc_bin = bloc_bin

        if prt == 'yes':
            print('States of the bloc:',bloc_num)
            print(bloc_matrix)
            print('Ground state energy of the matrix:',min(E))
            print()
            print()
    
    '''
    print(blocs_num[0])
    print(blocs_matrix[0][0])
    print()
    print(blocs_num[1])
    print(blocs_matrix[1][0])
    print(blocs_matrix[1][1])
    print()
    print(blocs_num[2])
    print(blocs_matrix[2][0])
    print(blocs_matrix[2][1])
    print(blocs_matrix[2][2])
    print()
    print(blocs_num[3])
    print(blocs_matrix[3][0])
    print(blocs_matrix[3][0])
    print()
    print(blocs_num[4])
    print(blocs_matrix[4][0])
    '''

    E,S = np.linalg.eigh(GS_bloc_matrix)
    if prt == 'yes':
        print()
        print()
        print('Ground state energy:',min(ground_states_energies))
        print('Ground state bloc:')
        print(GS_bloc_matrix)
        np.set_printoptions(precision=2)
        print('Energies:',E)
        print('States:')
        print(S)
        print()
    
    if time == 'yes':
        end = time.time()
        print('Time:',end - start,'seconds.')

    GS = np.asmatrix(S[:,0]).transpose()
    E = min(ground_states_energies)
    
    if not manip == 'yes':
        return GS,E,GS_bloc_matrix
    else:
        return blocs_matrix,blocs_num,GS_bloc_bin

