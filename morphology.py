import numpy as np
import random
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d
import collections
from shutil import copyfile
from kmc_classes import *
################################
# SET OF FUNCS THAT GENERATE THE MORPHOLOGY OF THE SYSTEM



#### CHOOSE A FUNC TO GENERATE EXCITONS

#num: number of excitons
#X: one of the position vectors    
def gen_pair_elechole(num, X):
    Ss, species = [], []
    selection = range(X)
    chosen = []
    while len(Ss) < num:
        number = random.choice(selection)
        Ss.append(Electron(number))
        number = random.choice(selection)
        Ss.append(Hole(number))
    #while len(Ss) < num:
        #number = random.choice(selection)
        #Ss.append(Hole(number))
        #number = random.choice(selection)
        #if number not in chosen:
        #    if random.uniform(0,1) < 0.5:
        #        Ss.append(Electron(number))
        #    else:
        #        Ss.append(Hole(number))
        #    #Ss.append(Exciton('singlet',number))
        #else:
        #    pass
    return Ss
    
def gen_excitons(num, X):
    Ss, species = [], []
    selection = range(X)
    chosen = []
    while len(Ss) < num:
        number = random.choice(selection)
	Ss.append(Electron(number))
    return Ss

generate_function = gen_pair_elechole





# FUNC TO GIVE S1 AND TRIPLET ENERGIES	
def homo_lumo(s1s, t1s, mats):
    s1 = []
    t1 = []
    for i in mats:
        s1.append(np.random.normal(s1s.get(i)[0],s1s.get(i)[1]))
        t1.append(np.random.normal(t1s.get(i)[0],t1s.get(i)[1]))
    return s1, t1


energy_function = homo_lumo
















