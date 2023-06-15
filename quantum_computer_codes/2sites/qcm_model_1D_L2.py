import pyqcm
from parameters import N

Nc = N

CM = pyqcm.cluster_model(Nc)
clus = pyqcm.cluster(CM,((0,0,0),(1,0,0)))
model = pyqcm.lattice_model('1D_L2', clus, ((100,0,0),))

model.interaction_operator('U')
model.hopping_operator('t', (1,0,0), -1)  # NN hopping

