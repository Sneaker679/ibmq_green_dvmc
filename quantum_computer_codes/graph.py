#!/usr/bin/env python3
# zlib license:

# Copyright (c) 2019 Maxime Charlebois

# This software is provided 'as-is', without any express or implied
# warranty. In no event will the authors be held liable for any damages
# arising from the use of this software.

# Permission is granted to anyone to use this software for any purpose,
# including commercial applications, and to alter it and redistribute it
# freely, subject to the following restrictions:

# 1. The origin of this software must not be misrepresented; you must not
#    claim that you wrote the original software. If you use this software
#    in a product, an acknowledgment in the product documentation would be
#    appreciated but is not required.
# 2. Altered source versions must be plainly marked as such, and must not be
#    misrepresented as being the original software.
# 3. This notice may not be removed or altered from any source distribution.

import copy
import time
import numpy as np
from numpy import linalg as la
import matplotlib.pyplot as plt
import sys, os, re
if len(sys.argv) == 2:
    number = sys.argv[1]
    sys.path.insert(0,os.path.join(os.path.dirname(__file__),'examples',number+'sites'))

from parameters import N,U,output_directory,pdf_output_directory
from hamiltonian_circuit import omega

from qiskit.quantum_info import Pauli,Operator
from qiskit.primitives import Estimator as pEstimator
from qiskit_nature.second_q.mappers import JordanWignerMapper
from qiskit_nature.second_q.operators import FermionicOp
from qiskit import QuantumCircuit,QuantumRegister
from qiskit import transpile
from scipy.linalg import eig, eigh, ordqz
from scipy.linalg.lapack import zggev
from ctypes import cdll, c_int, c_double

print('GS energy:',omega)
full_path = os.path.realpath(__file__)
pythonPathCode, file1 = os.path.split(full_path)

#U=4
Nb=0
Nc=N

trans_invariant = False
add_noise = False
qz_decomp = False
outputDir = output_directory

spectrumparaFileName='spectrumpara.def'
verbose_read = 1

sum_rule_max_ok = 1.01
sum_rule_min_ok = 0.96

tol = 1e-10
addtl_filter = 0
pct_filter = 1.0

if (len(sys.argv)>=2):
  spectrumparaFileName=sys.argv[1]
if (len(sys.argv)>=3):
  outputDir=sys.argv[2]+'/'
if (len(sys.argv)>=4):
  tol=float(sys.argv[3])
if (len(sys.argv)>=5):
  addtl_filter=int(sys.argv[4])
if (len(sys.argv)>=6):
  pct_filter=float(sys.argv[5])  


def dvmc_spectrum(Omega,verbose=1,fock_benchmarking = 'N'):

  sum_rule_max = sum_rule_max_ok
  sum_rule_min = sum_rule_min_ok

  # defaults
  w_min_data =-15.0
  w_max_data = 15.0
  eta = 0.1
  Nw = 2000

  Nsite   = Nc + Nb
  dw = (w_max_data-w_min_data)/(Nw-1)
  w_ = np.array(range(Nw))*dw + w_min_data

  stot = time.time()
  #s = time.time()
  if fock_benchmarking == 'Y':
      S_CA = np.load(os.path.join(outputDir,'S_CA_fock.npy'))
      S_AC = np.load(os.path.join(outputDir,'S_AC_fock.npy'))
      H_CA = np.load(os.path.join(outputDir,'H_CA_fock.npy'))
      H_AC = np.load(os.path.join(outputDir,'H_AC_fock.npy'))
  else:
      S_CA = np.load(os.path.join(outputDir,'S_CA.npy'))
      S_AC = np.load(os.path.join(outputDir,'S_AC.npy'))
      H_CA = np.load(os.path.join(outputDir,'H_CA.npy'))
      H_AC = np.load(os.path.join(outputDir,'H_AC.npy'))
  
  spectrum_hole = np.zeros([Nc,Nc,Nw])
  spectrum_elec = np.zeros([Nc,Nc,Nw])

  total_sum = 0.0
  partial_sum = 0.0
  k_label = u'%3s' % 'k#'

  #symmetrize
  H_AC = 0.5*(H_AC+np.transpose(np.conjugate(H_AC)))
  S_AC = 0.5*(S_AC+np.transpose(np.conjugate(S_AC)))
  H_CA = 0.5*(H_CA+np.transpose(np.conjugate(H_CA)))
  S_CA = 0.5*(S_CA+np.transpose(np.conjugate(S_CA)))

  s = time.time()
  # eigenvalue decomposition of overlap matrix, S
  SV_AC, U_SVD_AC = eigh(S_AC)
  SV_CA, U_SVD_CA = eigh(S_CA)
  e = time.time()
  print("Diagonalization 1 : ", e-s)

  tol_AC = np.abs(SV_AC.min())
  tol_CA = np.abs(SV_CA.min())


  nc_AC = len(SV_AC)
  break_AC = False
  
  break_CA = False
  SV_AC_tmp = SV_AC[::-1]
  SV_CA_tmp = SV_CA[::-1]

  import math
  cond_num_AC = math.log10(SV_AC_tmp[0]/np.abs(SV_AC_tmp[-1]))
  cond_num_CA = math.log10(SV_CA_tmp[0]/np.abs(SV_CA_tmp[-1]))  
  print("condition number S_AC = %8.6f"%(cond_num_AC))
  print("condition number S_CA = %8.6f"%(cond_num_CA))
  
  for i in range(len(SV_AC)):
    if(SV_AC_tmp[i]<tol and not break_AC):
      nc_AC = i
      break_AC = True
    if(SV_CA_tmp[i]<tol and not break_CA):
      nc_CA = i
      break_CA = True
    if(break_AC and break_CA):
      break

  print("number of states kept AC: ", nc_AC)
  print("number of states kept CA: ", nc_CA)

  if(addtl_filter):
    nc_AC = int(pct_filter*nc_AC)
    nc_CA = int(pct_filter*nc_CA)

    print("number of states kept AC (w/addtl filtering): ", nc_AC)
    print("number of states kept CA (w/addtl filtering): ", nc_CA)
    
  s = time.time()
  D_AC_sqrt = np.sqrt(SV_AC[(S_AC.shape[0]-nc_AC):])
  b_AC = la.inv(U_SVD_AC)[(S_AC.shape[0]-nc_AC):,:]
  c = (D_AC_sqrt*b_AC.T).T
  S_AC_sqrt = U_SVD_AC[:,(S_AC.shape[0]-nc_AC):].dot(c)

  D_CA_sqrt = np.sqrt(SV_CA[(S_CA.shape[0]-nc_CA):])
  b_CA = la.inv(U_SVD_CA)[(S_CA.shape[0]-nc_CA):,:]
  c = (D_CA_sqrt*b_CA.T).T
  S_CA_sqrt = U_SVD_CA[:,(S_CA.shape[0]-nc_CA):].dot(c)

  # compute S^-1/2 = U*P^+*(P*D*P^+)^-1/2*P*U^+
  D_AC_sqrt_inv = np.reciprocal(np.sqrt(SV_AC[(S_AC.shape[0]-nc_AC):]))
  c = (D_AC_sqrt_inv*b_AC.T).T
  Sinv_AC_sqrt = U_SVD_AC[:,(S_AC.shape[0]-nc_AC):].dot(c)

  D_CA_sqrt_inv = np.reciprocal(np.sqrt(SV_CA[(S_CA.shape[0]-nc_CA):]))
  c = (D_CA_sqrt_inv*b_CA.T).T
  Sinv_CA_sqrt = U_SVD_CA[:,(S_CA.shape[0]-nc_CA):].dot(c)

  e = time.time()
  print("Forming S, Sinv : ", e-s)
  # DO BROADCASTED MULTIPLY INSTEAD OF FULL MATMUL FOR PRODUCT ABOVE

  s = time.time()
  # Compute M = Sbar^(-1/2)*H*Sbar^(-1/2)
  Sinv_sqrt_H_Sinv_sqrt_AC = Sinv_AC_sqrt.dot(H_AC).dot(Sinv_AC_sqrt)
  Sinv_sqrt_H_Sinv_sqrt_CA = Sinv_CA_sqrt.dot(H_CA).dot(Sinv_CA_sqrt)
  e = time.time()
  print("Sinv*H*Sinv : ", e-s)

  s = time.time()
  # Solve standard eigenvalue problem in reduced Hilbert space
  # M = U_M * E * U_M^-1
  e_ac,u_ac = eigh(Sinv_sqrt_H_Sinv_sqrt_AC)
  e_ca,u_ca = eigh(Sinv_sqrt_H_Sinv_sqrt_CA)
  e = time.time()
  print("Diagonalization 2 : ", e-s)
  e_ac2 = e_ac - Omega 
  e_ca2 = -e_ca + Omega

  trim_inds = np.where(np.abs(e_ac)<1e-10)
  e_ac_r = np.delete(e_ac,trim_inds)
  u_ac_r = np.delete(u_ac,trim_inds,axis=1)
  e_ac2_r = e_ac_r - Omega

  trim_inds = np.where(np.abs(e_ca)<1e-10)
  e_ca_r = np.delete(e_ca,trim_inds)
  u_ca_r = np.delete(u_ca,trim_inds,axis=1)
  e_ca2_r = -e_ca_r + Omega


  s = time.time()
  # Compute Sbar^(1/2)*U_M
  u_ac = S_AC_sqrt.dot(u_ac_r)
  u_ca = S_CA_sqrt.dot(u_ca_r)
  e = time.time()
  print("Sbar^(1/2)*U_M : ", e-s)

  u_ac2 = u_ac[:Nc,:] 
  u_ca2 = u_ca[:Nc,:]    

  # Compute G(w) = Sbar^(1/2) * U_M * ((w + i * eta +- Omega) + E)^-1 * U_M^-1 * Sbar^(1/2)
  print("computing G(w) for range of frequencies")
  s = time.time()
  for ii in range(len(w_)):
    z_ac = w_[ii] + 1.j*eta + Omega;
    z_ca = w_[ii] + 1.j*eta - Omega;

    a = np.reciprocal(z_ac-e_ac_r)
    c1 = (a*u_ac2).T
    a = np.reciprocal(z_ca+e_ca_r)
    c2 = (a*u_ca2).T 

    G_ac = u_ac2.dot(c1)  # O(x3)
    G_ca = u_ca2.dot(c2)
    
    spectrum_elec[:,:,ii] = -G_ac[:,:].imag/(np.pi)
    spectrum_hole[:,:,ii] = -G_ca[:,:].imag/(np.pi)

  totalAij = spectrum_hole + spectrum_elec
  e = time.time()
  print("Calculating spectrum : ", e-s)

  for ii in range(0,2):
    sum = (w_[1]-w_[0])*np.sum(totalAij[ii,ii,:])
    print(sum)
  e = time.time()
  print("Execution time : ", e-stot)

  '''
  # get dos
  dos = np.zeros([3,Nw], dtype='float')
  for i in range(Nw):
    dos[0,i] = np.trace(totalAij[:,:,i])
    dos[1,i] = np.trace(spectrum_elec[:,:,i])
    dos[2,i] = np.trace(spectrum_hole[:,:,i])
  file_dos   = open(outputDir+'dos.dat','w')
  for ii in range(Nw):
    file_dos.write('% 7.6f   '  %w_[ii])
    for kk in range(3):
      file_dos.write('% 7.6f '  %(dos[kk,ii]))#/(total_sum)))
    file_dos.write('\n')
  '''

  # print solution for interface with QCM
  #print_solution_for_QCM(params,u_ac2,e_ac2_r,u_ca2,e_ca2_r)

  fig, ax = plt.subplots()

  shift = 0.5
  for ii in range(Nc):
  #for ii in range(2):
    ax.plot(w_,totalAij[ii,ii,:]+ii*shift)
    ax.plot([w_.min(),w_.max()],[ii*shift,ii*shift],c='black')
  
  # local dos
  if fock_benchmarking == 'Y':
    dos_name = 'local_dos_fock.dat'
    file_name = 'fock_spectrum.pdf'
  else:
    dos_name = 'local_dos.dat'
    file_name = 'ibmq_spectrum.pdf'

  file_dos   = open(os.path.join(outputDir,dos_name),'w')
  for ii in range(Nw):
    file_dos.write('% 7.6f   '  %w_[ii])
    for kk in range(Nc):
      file_dos.write('% 7.6f '  %(totalAij[kk,kk,ii]))#/(total_sum)))
    file_dos.write('\n')
    
  #ax.set_xlim(w_.min(),w_.max())
  ax.axvline(x=0)
  ax.set_xlim(-10,10)

  #matplotlib.use('Agg')
  plt.savefig(os.path.join(pdf_output_directory,file_name))
  #plt.show()
  
  sys.exit()  

if __name__ == '__main__':
    dvmc_spectrum(omega,verbose_read)
