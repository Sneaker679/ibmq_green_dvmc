import matplotlib.pyplot as plt
import numpy as np

'''
plt.plot(*np.loadtxt("local_dos.dat",unpack=True,usecols=(0,1)), linewidth=2.0)
'''

num_rows, num_cols = np.loadtxt("local_dos_qcm.dat",unpack=True).shape
num_rows_qcm, num_cols_qcm = np.loadtxt("local_dos_qcm.dat",unpack=True).shape

if num_rows != num_rows_qcm:
    raise Exception('Number of sites is not the same for the .dat files.')

for site in range(num_rows-1):
    step = 0.5
    
    array_y = [value+site*step for value in np.loadtxt("local_dos.dat",unpack=True,usecols=(0,site+1))[1]]
    array_x = np.loadtxt("local_dos.dat",unpack=True,usecols=(0,site+1))[0]
    line1, = plt.plot(array_x,array_y, 'r', alpha=1, linewidth=2.0, label='ibmq')

    array_y = [value+site*step for value in np.loadtxt("local_dos_qcm.dat",unpack=True,usecols=(0,site+1))[1]]
    array_x = np.loadtxt("local_dos_qcm.dat",unpack=True,usecols=(0,site+1))[0]
    line2, = plt.plot(array_x,array_y, 'k',linewidth=1.0,label='Exact with qcm')
    
    plt.legend(handles=[line1, line2], loc='upper right')

#plt.show()
plt.ylabel( r"$A_{ii} (\omega)$")
plt.xlabel('$\omega$')
plt.xlim([-6,6])

plt.savefig('combined_spectrum.pdf')
