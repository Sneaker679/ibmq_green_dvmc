### Packages ################################################
import pyqcm
import numpy as np
import sys,os

### Fetching parameters.py and hamiltonian_circuit.py #######
module_directory = os.path.dirname(__file__)
sys.path.insert(0,module_directory)

working_directory = os.getcwd()
sys.path.insert(0,working_directory)

from parameters import N,t,U,mu,output_directory,pdf_output_directory
from hamiltonian_circuit import model


# Making a list of sectors where the GS probably is.
sec = 'R0:N0:S0'
for sectors in range(2*N):
    sectors += 1
    if sectors % 2 == 0:
        sec += f'/R0:N{sectors}:S0'

# Targeting sectors and setting parameters.
model.set_target_sectors([sec])
model.set_parameters(f"""
t={-t}
U = {U}
mu = {mu}
""")

# Calculating the Ground State Energy.
I = pyqcm.model_instance(model)
print(I.ground_state())

# Generating PDF
result = I.cluster_spectral_function(file= os.path.join(pdf_output_directory,'qcm_spectrum.pdf'),eta=0.1,wmax=15)

# Generating local dos
file_dos_qcm   = open(os.path.join(output_directory,'local_dos_qcm.dat'),'w')
w_ = result[0]
local_dos = result[1]
for ii in range(len(w_)):
    file_dos_qcm.write('% 7.6f   '  %w_[ii])
    for kk in range(N):
        file_dos_qcm.write('% 7.6f '  %(local_dos[ii,kk]/np.pi))#/(total_sum)))
    file_dos_qcm.write('\n')
file_dos_qcm.close()
