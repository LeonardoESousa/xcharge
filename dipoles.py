import numpy as np
import random


vecs = np.zeros((1,3))
for i in range(1000):
    theta = random.uniform(0,np.pi)
    phi   = random.uniform(0,2*np.pi)
    x = 3*np.sin(theta)*np.cos(phi)
    y = 3*np.sin(theta)*np.sin(phi)
    z = 3*np.cos(theta) 
    nvec = np.array([x,y,z])
    vecs = np.vstack((vecs,nvec))

np.savetxt('dipoles.txt', vecs[1:,:], delimiter='\t')    