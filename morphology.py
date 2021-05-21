import numpy as np
import random
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d
import collections
from kmc_classes import *

# SET OF FUNCS THAT GENERATE THE MORPHOLOGY OF THE SYSTEM
# note: always define the function by list (param) that contains the things needed


# reads the morphology
def read_lattice(file_name):
	X,Y,Z,Mats = [], [], [], []
	
	with open(file_name, 'r') as f:
		for line in f:
			line = line.split()
			
			x    = float(line[0])
			y    = float(line[1])
			z    = float(line[2])
			mat  = int(float(line[3]))
			
			
			X.append(x)
			Y.append(y)
			Z.append(z)
			Mats.append(mat)
	X = np.array(X)
	Y = np.array(Y)	
	Z = np.array(Z)	
	Mats = np.array(Mats)
	return X,Y,Z,Mats


#### CHOOSE A FUNC TO GENERATE EXCITONS
def gen_pair_elechole(param):
    #num: number of excitons
    #X: one of the position vectors 
    
    num    = param[0] #relabeling the input parameters 
    X      = param[1]
    
    
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
    
def gen_excitons(param):
    num    = param[0] #relabeling the input parameters 
    X      = param[1]
    
    Ss, species = [], []
    selection = range(X)
    chosen = []
    while len(Ss) < num:
        number = random.choice(selection)
        Ss.append(Exciton('singlet',number))
    return Ss

##########################################
# FUNC TO GIVE S1 AND TRIPLET ENERGIES	
def homo_lumo(param):

    s1s    = param[0] #relabeling the input parameters 
    t1s    = param[1]
    mats   = param[2]
    
    s1 = []
    t1 = []
    for i in mats:
        s1.append(np.random.normal(s1s.get(i)[0],s1s.get(i)[1]))
        t1.append(np.random.normal(t1s.get(i)[0],t1s.get(i)[1]))
    return s1, t1
############################################
# ANNI FUNCS NOTE: ALL THEM MUST HAVE THE SAME VARIABLES (system,tipos,Ss,indices,locs)

#annihilation electron-hole pair
def anni_ele_hol(system,tipos,Ss,indices,locs):
    if 'electron' in tipos and 'hole' in tipos:        
        Ss[indices[0][tipos.index('electron')]].kill('anni',system,system.s1)
        Ss[indices[0][tipos.index('hole')]].kill('anni',system,system.s1)
        if random.uniform(0,1) <= 0.75:
            system.add_particle(Exciton('triplet',locs[indices[0][0]]))
        else:
            system.add_particle(Exciton('singlet',locs[indices[0][0]]))
#annihililation exciton singlets pair
def anni_sing(system,tipos,Ss,indices,locs): # e se tiver 4 excitons no mesmo sitio?
    duplicates = set([x for x in tipos if tipos.count(x) > 1]) # checking if there is 2 occurrences of the types
    if 'singlet' in duplicates:        
        Ss[indices[0][tipos.index('singlet')]].kill('anni',system,system.s1)
############################################           
                











