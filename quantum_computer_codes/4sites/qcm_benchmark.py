import pyqcm
import numpy as np
import parameters as M
from parameters import t,U,mu,N

sec = f'R0:N2:S0'
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
    for kk in range(N):
        file_dos_qcm.write('% 7.6f '  %(local_dos[ii,kk]/np.pi))#/(total_sum)))
    file_dos_qcm.write('\n')
file_dos_qcm.close()
