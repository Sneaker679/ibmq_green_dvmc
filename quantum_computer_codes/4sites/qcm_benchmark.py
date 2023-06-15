import pyqcm
import qcm_model_2D_C4 as M
import numpy as np
from qcm_model_2D_C4 import Nc
from parameters import t,U,mu

sec = f'R0:N{Nc}:S0'
M.model.set_target_sectors([sec])
M.model.set_parameters(f"""
t={-t}
U = {U}
mu = {mu}
""")
I = pyqcm.model_instance(M.model)
  
print(I.ground_state())

result = I.cluster_spectral_function(file='qcm_spectrum.pdf',eta=0.1,wmax=15)

# local dos
file_dos_qcm   = open('local_dos_qcm.dat','w')
w_ = result[0]
local_dos = result[1]
for ii in range(len(w_)):
    file_dos_qcm.write('% 7.6f   '  %w_[ii])
    for kk in range(Nc):
        file_dos_qcm.write('% 7.6f '  %(local_dos[ii,kk]/np.pi))#/(total_sum)))
    file_dos_qcm.write('\n')
file_dos_qcm.close()
