#!/usr/bin/env python3
### Pacakges ################################################
import matplotlib.pyplot as plt
import numpy as np
import sys,os

### Fetching parameters.py ##################################
working_directory = os.getcwd()
sys.path.insert(0,working_directory)

from parameters import (
    N,t,U,mu,
    output_directory,pdf_output_directory,
    use_qcm,graph_for_qcm,graph_for_ibmq,graph_for_fock,
)

# Finding the paths of the .npy files to load and fetching number of sites
if graph_for_ibmq is False and graph_for_qcm is False and graph_for_fock is False:
    print('Skipping the combining of the graphs since all graph_for... are False in parameters.py.')
    sys.exit()

if graph_for_ibmq is True:
    to_load = os.path.join(output_directory,'local_dos.dat')
    num_rows, num_cols = np.loadtxt(to_load,unpack=True).shape
if use_qcm is True and graph_for_qcm is True:
    to_load_qcm = os.path.join(output_directory,'local_dos_qcm.dat')
    num_rows, num_cols = np.loadtxt(to_load_qcm,unpack=True).shape
if graph_for_fock is True:
    to_load_fock = os.path.join(output_directory,'local_dos_fock.dat')
    num_rows, num_cols = np.loadtxt(to_load_fock,unpack=True).shape

# Find spectrum borders
aij,w=np.loadtxt(to_load_fock,unpack=True)[1,:],np.loadtxt(to_load_fock,unpack=True)[0,:]
tol = 0.001
for i in zip(aij,w):
    if i[0] > tol:
        left_border_max_frequency = i[1]
        break
for i in zip(aij[::-1],w[::-1]):
    if i[0] > tol:
        right_border_max_frequency = i[1]
        break

# Finding the number of rows (sites)
for site in range(num_rows-1):
    # Step between the lines of each site
    step = 1
    lines = []

    # Drawing graph for ibmq
    if graph_for_ibmq is True:
        array_y = [value+site*step for value in np.loadtxt(to_load,unpack=True,usecols=(0,site+1))[1]]
        array_x = np.loadtxt(to_load,unpack=True,usecols=(0,site+1))[0]
        line1, = plt.plot(array_x,array_y, 'r', alpha=1, linewidth=2.0, label='ibmq')
        lines.append(line1)

    # Drawing graph for qcm
    if use_qcm is True and graph_for_qcm is True:
        array_y = [value+site*step for value in np.loadtxt(to_load_qcm,unpack=True,usecols=(0,site+1))[1]]
        array_x = np.loadtxt(to_load_qcm,unpack=True,usecols=(0,site+1))[0]
        line2, = plt.plot(array_x,array_y, 'k',linewidth=1.0,label='Exact with qcm')
        lines.append(line2)
    
    # Drawing graph for fock
    if graph_for_fock is True:
        array_y = [value+site*step for value in np.loadtxt(to_load_fock,unpack=True,usecols=(0,site+1))[1]]
        array_x = np.loadtxt(to_load_fock,unpack=True,usecols=(0,site+1))[0]
        line3, = plt.plot(array_x,array_y, 'b',linewidth=0.2,label='Exact with fock')
        lines.append(line3)
    
    plt.legend(handles=lines)

plt.ylabel( r"$A_{ii} (\omega)$")
plt.xlabel('$\omega$')
plt.xlim([left_border_max_frequency,right_border_max_frequency])
left_border_max_frequency = -10
plt.title(f'Green Function Spectrum when N={N},t={t},U={U} and mu={mu}')
plt.savefig(os.path.join(pdf_output_directory,'combined_spectrum.pdf'))
