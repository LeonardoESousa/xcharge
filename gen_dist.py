import numpy as np
import random
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d
import collections
from math import sqrt

nele = 100
'''
#random
for i in range(2):
	name = 'dist_'+str(i)+'.txt'
	sdist = np.random.normal(5,1,nele)
	np.savetxt(name,sdist)
'''
for i in range(2):
	name = 'dist_'+str(i)+'.txt'
	sdist = np.array([ i for j in range(nele) ])
	np.savetxt(name,sdist)
