import sys
import numpy as np

class Fock():

    def __init__(self,N,custom_input = None):
        if N < 1:
            raise Exception('Number of sites (N) must be at least 1.')
        self.N = N
        self.sign = 1
        self.total_spin = 0
    
        if custom_input is None:
            self.fock = np.zeros((1,2*N)).astype(int)
            self.num = 0
            self.num_electrons = 0
        
        if custom_input is not None:
            if custom_input < 0:
                raise Exception('The custom input for the state number must be a positive integer.')
            self.fock = np.matrix([int(i) for i in str(np.binary_repr(custom_input,N*2))])
            self.num =  int(''.join(str(x) for x in np.nditer(self.fock)),2)
            self.num_electrons = np.count_nonzero(self.fock == 1)
            for index,value in np.ndenumerate(self.fock):
                if value != 0:
                    if index[1] < self.fock.shape[1]/2:
                        self.total_spin += 1
                    else:
                        self.total_spin += -1

    def __str__(self):
        return str(self.fock)

    def op(self,type,site,spin):
        if spin not in ['+','-']:
            raise Exception('Spin must be "+" or "-".')
        if site < 0:
            raise Exception('Site numbers start at 0.')
        if not type == 'create' and not type == 'destroy':
            raise Exception('The only two operators are "create" and "destroy".')        

        if not isinstance(self.fock,int):
            
            if type == 'create':
                checked_element = 0
                added_element = 1
                spin_controller_up = 1
                spin_controller_down = -1
            else:
                checked_element = 1
                added_element = 0
                spin_controller_up = -1
                spin_controller_down = 1

            if spin == '-':
                sp = int(self.fock[0].size/2)
                self.total_spin += spin_controller_down
            else:
                sp = 0
                self.total_spin += spin_controller_up
            
            del spin_controller_up
            del spin_controller_down

            if self.fock[0,site+sp] == checked_element:
                partial_bin = ''.join(str(x) for x in self.fock.tolist()[0][:site+sp])
                sign_power = 0
                for number in partial_bin:
                    if number == '1':
                        sign_power += 1
                self.sign = self.sign * (-1)**sign_power
                self.fock[0,site+sp] = added_element
            else:
                self.fock = 0
        
            self.num =  int(''.join(str(x) for x in np.nditer(self.fock)),2)
            self.num_electrons = np.count_nonzero(self.fock == 1)
        
        
        if isinstance(self.fock,int):
            self.num = 0
            self.sign = 0
            self.num_electrons = 0
            self.total_spin = 0


