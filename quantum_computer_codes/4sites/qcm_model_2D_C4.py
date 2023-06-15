import pyqcm
from parameters import N

Nc = N

CM = pyqcm.cluster_model(Nc)
clus = pyqcm.cluster(CM,((0,0,0),(1,0,0),(0,1,0),(1,1,0)))
model = pyqcm.lattice_model('2D_C4', clus, ((100,0,0),))

model.interaction_operator('U')
model.hopping_operator('t', (1,0,0), -1)  # NN hopping
model.hopping_operator('t', (0,1,0), -1)  # NN hopping

