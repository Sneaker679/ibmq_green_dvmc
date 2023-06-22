import pyqcm
import sys
if len(sys.argv) == 2:
    sys.path.insert(0,'./examples/'+sys.argv[1]+'sites')

from parameters import N,t,U,mu
from hamiltonian_circuit import model

import numpy as np

sec = f'R0:N{N}:S0'
model.set_target_sectors([sec])
model.set_parameters(f"""
t={-t}
U = {U}
mu = {mu}
""")
I = pyqcm.model_instance(model)
  
print(I.ground_state())


if len(sys.argv) == 2:
    result = I.cluster_spectral_function(file='./examples/'+sys.argv[1]+'sites/qcm_spectrum.pdf',eta=0.1,wmax=15)
else:
    result = I.cluster_spectral_function(file='qcm_spectrum.pdf',eta=0.1,wmax=15)

# local dos
if len(sys.argv) == 2:
    file_dos_qcm   = open('examples/'+sys.argv[1]+'sites/matrices_npy/local_dos_qcm.dat','w')
else:
    file_dos_qcm   = open('matrices_npy/local_dos_qcm.dat','w')
w_ = result[0]
local_dos = result[1]
for ii in range(len(w_)):
    file_dos_qcm.write('% 7.6f   '  %w_[ii])
    for kk in range(N):
        file_dos_qcm.write('% 7.6f '  %(local_dos[ii,kk]/np.pi))#/(total_sum)))
    file_dos_qcm.write('\n')
file_dos_qcm.close()
