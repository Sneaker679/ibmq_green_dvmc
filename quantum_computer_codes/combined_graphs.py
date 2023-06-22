import matplotlib.pyplot as plt
import numpy as np
from parameters import N,t,U,mu
import sys

'''
plt.plot(*np.loadtxt("local_dos.dat",unpack=True,usecols=(0,1)), linewidth=2.0)
'''

if len(sys.argv) == 2:
    to_load = './examples/'+sys.argv[1]+'sites/matrices_npy/local_dos.dat'
    to_load_qcm = './examples/'+sys.argv[1]+'sites/matrices_npy/local_dos_qcm.dat'
else:
    to_load = './matrices_npy/local_dos.dat'
    to_load_qcm = './matrices_npy/local_dos_qcm.dat'

num_rows, num_cols = np.loadtxt(to_load,unpack=True).shape
num_rows_qcm, num_cols_qcm = np.loadtxt(to_load_qcm,unpack=True).shape

if num_rows != num_rows_qcm:
    raise Exception('Number of sites is not the same for the .dat files.')

for site in range(num_rows-1):
    step = 0.5
    
    array_y = [value+site*step for value in np.loadtxt(to_load,unpack=True,usecols=(0,site+1))[1]]
    array_x = np.loadtxt(to_load,unpack=True,usecols=(0,site+1))[0]
    line1, = plt.plot(array_x,array_y, 'r', alpha=1, linewidth=2.0, label='ibmq')

    array_y = [value+site*step for value in np.loadtxt(to_load_qcm,unpack=True,usecols=(0,site+1))[1]]
    array_x = np.loadtxt(to_load_qcm,unpack=True,usecols=(0,site+1))[0]
    line2, = plt.plot(array_x,array_y, 'k',linewidth=1.0,label='Exact with qcm')
    
    plt.legend(handles=[line1, line2])

#plt.show()
plt.ylabel( r"$A_{ii} (\omega)$")
plt.xlabel('$\omega$')
plt.xlim([-10,10])

if len(sys.argv) == 2:
    plt.savefig('./examples/'+sys.argv[1]+'sites/combined_spectrum.pdf')
else:
    plt.savefig('combined_spectrum.pdf')

plt.title(f'N={N},t={t},U={U},mu={mu}')
