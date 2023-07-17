### Packages ################################################
import numpy as np
import sys,os

### Fetching parameters.py and hamiltonian_circuit.py #######
module_directory = os.path.dirname(__file__)
sys.path.insert(0,module_directory)

working_directory = os.getcwd()
sys.path.insert(0,working_directory)

from parameters import use_qcm,N,t,U,mu,spin_green,spin_gs,output_directory,pdf_output_directory

if use_qcm is False:
    print()
    print('Skipping the qcm benchmark since use_qcm == "N" in parameters.py.')
    print()
    sys.exit()

import pyqcm
from hamiltonian_circuit import model

# Making a list of sectors where the GS probably is.
sec = ''

def new_sector(N,Ne):
    sec = f'/R0:N{Ne}:S{S}'
    if not N == Ne:
        new_Ne = 2*N-Ne
        sec += f'/R0:N{new_Ne}:S{S}'
    return sec

for Ne in range(N+1):
    S = -Ne
    for m in range(Ne+1):
        if spin_gs == '+':
            if S >= 0:
                sec += new_sector(N,Ne)
        if spin_gs == '-':
            if S <= 0:
                sec += new_sector(N,Ne)

        S = S+2
    Ne += 1
sec = sec[1:]

# Targeting sectors and setting parameters.
model.set_target_sectors([sec])
model.set_parameters(f"""
t={-t}
U = {U}
mu = {mu}
""")

# Calculating the Ground State Energy.
I = pyqcm.model_instance(model)
gs = I.ground_state(pr=False)

previous_S_value = 64
if len(gs[0][1].split('/')) > 1:
    new_sectors = gs[0][1].split('/')
    for index,sector in enumerate(new_sectors):
        sector_list = sector.split(':')
        S_value = int(sector_list[2][1:])
        if ((spin_gs == '+' and S_value < previous_S_value and S_value >= 0)
        or (spin_gs == '-' and S_value > previous_S_value and S_value <= 0)):
            desired_block_index = index
            previous_S_value = S_value

    model.set_target_sectors([f'{new_sectors[desired_block_index]}'])
    I = pyqcm.model_instance(model)
    gs = I.ground_state(pr=False)

print(gs)

# Generating PDF
if spin_green == '+':
    spin_down = False
else:
    spin_down = True

result = I.cluster_spectral_function(blocks=True,file= os.path.join(pdf_output_directory,'qcm_spectrum.pdf'),eta=0.1,wmax=15,spin_down=spin_down)

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
