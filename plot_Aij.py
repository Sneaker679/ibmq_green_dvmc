import numpy as np
import pyqcm
from pyqcm import *
from pyqcm import spectral

def model1D(ns, nb, S):
    if ns%2 :
        raise ValueError('model1D should have an even number of cluster sites')
    if nb%2 :
        raise ValueError('model1D should have an even number of bath sites')
    no = ns+nb
    nb2 = nb//2

    new_cluster_model('clus', ns, nb)

    if nb != 0 :
        model_name = '{:d}L_{:d}b_gen'.format(ns,nb)
        for i in range(1, nb+1):
            new_cluster_operator('clus', 'eb{:d}'.format(i), 'one-body', [(i+ns,i+ns,1), (i+ns+no,i+ns+no,1)])
            new_cluster_operator('clus', 'tb{:d}'.format(i), 'one-body', [(1,i+ns,1), (1+no,i+ns+no,1),(ns,i+ns,S[i-1]), (ns+no,i+ns+no,S[i-1])])
    else:
        model_name = '{:d}L'.format(ns)

    add_cluster('clus', [0,0,0], [[i,0,0] for i in range(ns)])
    lattice_model(model_name, [[ns,0,0]])

    interaction_operator('U')
    hopping_operator('t', [1, 0, 0], -1)  # NN hopping

ns = 2
nb = 0

model1D(ns,nb,[])

from pyqcm.cdmft import *
import sys

Nelec = 2

set_target_sectors(['R0:N{:d}:S0'.format(Nelec)])

basic_params="""
t = 1
U = 8
mu = 4.0
"""

set_parameters(basic_params)

pyqcm.new_model_instance()

pyqcm.spectral.cluster_spectral_function(file='spectrum_rspace_ED.pdf',eta=0.1,wmax=15)
