import sys
import numpy as np

class Fock():

    def __init__(self,N,custom_input = None,qiskit_notation='N'):
        if N < 1:
            raise Exception('Number of sites (N) must be at least 1.')
        if not qiskit_notation == 'Y' and not qiskit_notation == 'N':
            raise Exception('qiskit_notation must be "Y" or "N".')
        self.N = N
        self.sign = 1
        self.total_spin = 0
        self.qiskit_notation = qiskit_notation
    
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
                    if self.qiskit_notation == 'N':
                        if index[1] < self.fock.shape[1]/2:
                            self.total_spin += 1
                        else:
                            self.total_spin += -1
                    if self.qiskit_notation == 'Y':
                        if index[1] % 2 == 0:
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
            
            checked_element = 0         if type == 'create' else 1
            added_element = 1           if type == 'create' else 0
            spin_controller_up = 1      if type == 'create' else -1
            spin_controller_down = -1   if type == 'create' else 1

            if spin == '-':
                sp = 1          if self.qiskit_notation == 'Y'  else int(self.fock[0].size/2)
            if spin == '+':
                sp = 0          if self.qiskit_notation == 'Y'  else 0
           
            index = 2*site+sp   if self.qiskit_notation == 'Y'  else site+sp

            self.total_spin += spin_controller_down if spin == '-' else spin_controller_up

            if self.fock[0,index] == checked_element:
                partial_bin = ''.join(str(x) for x in self.fock.tolist()[0][:index])
                sign_power = 0
                for number in partial_bin:
                    if number == '1':
                        sign_power += 1
                self.sign = self.sign * (-1)**sign_power
                self.fock[0,index] = added_element
            else:
                self.fock = 0
        
            self.num =  int(''.join(str(x) for x in np.nditer(self.fock)),2)
            self.num_electrons = np.count_nonzero(self.fock == 1)
        
        
        if isinstance(self.fock,int):
            self.num = 0
            self.sign = 0
            self.num_electrons = 0
            self.total_spin = 0
