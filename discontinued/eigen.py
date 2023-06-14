import numpy as np
U = 4
t = -1
mu = U/2

test = np.array([
                [U-2*mu,t,t,0],
                [t,-2*mu,0,t],
                [t,0,-2*mu,t],
                [0,t,t,U-2*mu]
])
E,S = np.linalg.eigh(test)

print('E:',E)
print('S:',S)
