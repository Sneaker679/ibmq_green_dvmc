### Pacakges ################################################
import matplotlib.pyplot as plt
import numpy as np
import sys,os

### Fetching parameters.py ##################################
working_directory = os.getcwd()
sys.path.insert(0,working_directory)

from parameters import N,t,U,mu,output_directory,pdf_output_directory,use_qcm


# Finding the paths of the .npy files to load
to_load = os.path.join(output_directory,'local_dos.dat')
if use_qcm is True:
    to_load_qcm = os.path.join(output_directory,'local_dos_qcm.dat')
to_load_fock = os.path.join(output_directory,'local_dos_fock.dat')

# Finding the number of rows (sites)
num_rows, num_cols = np.loadtxt(to_load,unpack=True).shape

for site in range(num_rows-1):
    # Step between the lines of each site
    step = 0.5
    lines = []

    # Drawing graph for ibmq
    array_y = [value+site*step for value in np.loadtxt(to_load,unpack=True,usecols=(0,site+1))[1]]
    array_x = np.loadtxt(to_load,unpack=True,usecols=(0,site+1))[0]
    line1, = plt.plot(array_x,array_y, 'r', alpha=1, linewidth=2.0, label='ibmq')
    lines.append(line1)

    # Drawing graph for qcm
    if use_qcm is True:
        array_y = [value+site*step for value in np.loadtxt(to_load_qcm,unpack=True,usecols=(0,site+1))[1]]
        array_x = np.loadtxt(to_load_qcm,unpack=True,usecols=(0,site+1))[0]
        line2, = plt.plot(array_x,array_y, 'k',linewidth=1.0,label='Exact with qcm')
        lines.append(line2)
    
    # Drawing graph for fock
    array_y = [value+site*step for value in np.loadtxt(to_load_fock,unpack=True,usecols=(0,site+1))[1]]
    array_x = np.loadtxt(to_load_fock,unpack=True,usecols=(0,site+1))[0]
    line3, = plt.plot(array_x,array_y, 'b',linewidth=0.2,label='Exact with fock')
    lines.append(line3)
    
    plt.legend(handles=lines)

plt.ylabel( r"$A_{ii} (\omega)$")
plt.xlabel('$\omega$')
plt.xlim([-10,10])
plt.title(f'N={N},t={t},U={U},mu={mu}')
plt.savefig(os.path.join(pdf_output_directory,'combined_spectrum.pdf'))
