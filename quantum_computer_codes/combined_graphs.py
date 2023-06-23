import matplotlib.pyplot as plt
import numpy as np
import sys,os
if len(sys.argv) == 2:
    number = sys.argv[1]
    sys.path.insert(0,os.path.join(os.path.dirname(__file__),'examples',number+'sites'))

from parameters import N,t,U,mu,output_directory,pdf_output_directory,use_qcm

to_load = os.path.join(output_directory,'local_dos.dat')
if use_qcm == 'Y':
    to_load_qcm = os.path.join(output_directory,'local_dos_qcm.dat')
to_load_fock = os.path.join(output_directory,'local_dos_fock.dat')

num_rows, num_cols = np.loadtxt(to_load,unpack=True).shape
if use_qcm == 'Y':
    num_rows_qcm, num_cols_qcm = np.loadtxt(to_load_qcm,unpack=True).shape
num_rows_fock, num_cols_fock = np.loadtxt(to_load_fock,unpack=True).shape

for site in range(num_rows-1):
    step = 0.5
    lines = []

    array_y = [value+site*step for value in np.loadtxt(to_load,unpack=True,usecols=(0,site+1))[1]]
    array_x = np.loadtxt(to_load,unpack=True,usecols=(0,site+1))[0]
    line1, = plt.plot(array_x,array_y, 'r', alpha=1, linewidth=2.0, label='ibmq')
    lines.append(line1)

    if use_qcm == 'Y':
        array_y = [value+site*step for value in np.loadtxt(to_load_qcm,unpack=True,usecols=(0,site+1))[1]]
        array_x = np.loadtxt(to_load_qcm,unpack=True,usecols=(0,site+1))[0]
        line2, = plt.plot(array_x,array_y, 'k',linewidth=1.0,label='Exact with qcm')
        lines.append(line2)
    
    array_y = [value+site*step for value in np.loadtxt(to_load_fock,unpack=True,usecols=(0,site+1))[1]]
    array_x = np.loadtxt(to_load_fock,unpack=True,usecols=(0,site+1))[0]
    line3, = plt.plot(array_x,array_y, 'b',linewidth=0.2,label='Exact with fock')
    lines.append(line3)
    
    plt.legend(handles=lines)

#plt.show()
plt.ylabel( r"$A_{ii} (\omega)$")
plt.xlabel('$\omega$')
plt.xlim([-10,10])
plt.title(f'N={N},t={t},U={U},mu={mu}')
plt.savefig(os.path.join(pdf_output_directory,'combined_spectrum.pdf'))

