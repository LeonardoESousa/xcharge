import numpy as np
import random
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d
import collections
from kmc_classes import *
from math import sqrt
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
def gen_electron(param):
    #num: number of excitons
    #X: one of the position vectors 
    
    num       = param[0] #relabeling the input parameters 
    selection = param[1]
       
    Ss, species = [], []
    chosen = []
    while len(Ss) < num:
        number = random.choice(selection)
        Ss.append(Electron(number))
    return Ss
    
def gen_hole(param):
    #num: number of excitons
    #X: one of the position vectors 
    
    num       = param[0] #relabeling the input parameters 
    selection = param[1]
       
    Ss, species = [], []

    chosen = []
    while len(Ss) < num:
        number = random.choice(selection)
        Ss.append(Hole(number))
    return Ss
    
def gen_pair_elechole(param):
    #num: number of excitons
    #X: one of the position vectors 
    
    num       = param[0] #relabeling the input parameters 
    selection = param[1]
    
    
    Ss, species = [], []
    
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
    num       = param[0] #relabeling the input parameters 
    selection = param[1]
       
    Ss, species = [], []
    
    chosen = []
    while len(Ss) < num:
        number = random.choice(selection)
        Ss.append(Exciton('singlet',number))
    return Ss

def gen_excitonsv2(param):
    num       = param[0] #relabeling the input parameters 
    selection = param[1]
       
    Ss, species = [], []
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
        id1 = Ss[indices[0][tipos.index('electron')]].identity
        id2 = Ss[indices[0][tipos.index('hole')]].identity
        Ss[indices[0][tipos.index('electron')]].kill('anni',system,system.s1)
        Ss[indices[0][tipos.index('hole')]].kill('anni',system,system.s1)
        if random.uniform(0,1) <= 0.75 and abs(id1) != abs(id2):
            system.add_particle(Exciton('triplet',locs[indices[0][0]]))
        else:
            system.add_particle(Exciton('singlet',locs[indices[0][0]]))
#annihililation exciton singlets pair
def anni_sing(system,tipos,Ss,indices,locs): # e se tiver 4 excitons no mesmo sitio?
    duplicates = set([x for x in tipos if tipos.count(x) > 1]) # checking if there is 2 occurrences of the types
    if 'singlet' in duplicates:        
        Ss[indices[0][tipos.index('singlet')]].kill('anni',system,system.s1)
############################################           
# FUNCS TO SHAPE THE GENERATION OF PARTICLES

#main function, regardless of the conditional, filters the sites according to it	
def filter_selection(X,Y,Z,Mats,shape_dic,**kwargs):

	X_cop     = X.copy()
	Y_cop     = Y.copy()
	Z_cop     = Z.copy()
	Mats_cop  = Mats.copy()
	selection = []
	
	
	mat_restriction     = kwargs['mat']   #list of materials that will generate part
	shape_type          = kwargs['shape'] #the shape that the particles will be aranged
	argument_shape_func = kwargs['argum'] #arguments of the shape func
	origin              = kwargs['origin'] #origin of the shape
	conditional = shape_dic.get(shape_type)
	
	for index_mol in range(len(X)):

		x   = X_cop[index_mol]
		y   = Y_cop[index_mol]
		z   = Z_cop[index_mol]
		mat = Mats_cop[index_mol]
		
		if conditional([x,y,z],origin,argument_shape_func):
			if (mat in mat_restriction) or (mat_restriction[0] == None):
				selection.append(index_mol)
	return selection                


#funcs that define the filtering, like inside a sphere, plane etc

def sphere_conditional(pos,r0,r_min):
	x  = pos[0] - r0[0]
	y  = pos[1] - r0[1]
	z  = pos[2] - r0[2]
	
	r_pos = sqrt(x**2 + y**2 + z**2)
	
	return r_pos <= r_min #condition to be inside the sphere
	
def plane_conditional(pos,r0,COEF):

	#Equation of plane is defined as 
	# A(X-X0) + B(Y-Y0) + C(Z-Z0) = 0
	x  = pos[0] - r0[0]
	y  = pos[1] - r0[1]
	z  = pos[2] - r0[2]
	
	A = COEF[0]
	B = COEF[1]
	C = COEF[2]
	
	cond = A*x + B*y + C*z == 0 #condition to be in the plane
	
	return cond
	
def cone_conditional(pos,r0,COEF): #hallow cone

	#Equation of plane is defined as 
	# A(X-X0) + B(Y-Y0) + C(Z-Z0) = 0
	x  = pos[0] - r0[0]
	y  = pos[1] - r0[1]
	z  = pos[2] - r0[2]
	
	A = COEF[0]
	B = COEF[1]
	C = COEF[2]
	
	cond = (x**2)/(A**2) + (y**2)/(B**2) -(z**2)/(C**2) == 0 #condition to be in the plane
	
	return cond	
	
def cilinder_conditional(pos,r0,COEF):

	#Equation of plane is defined as 
	# A(X-X0) + B(Y-Y0) + C(Z-Z0) = 0
	x  = pos[0] - r0[0]
	y  = pos[1] - r0[1]
	z  = pos[2] - r0[2]
	
	RHO_min   = COEF
	
	cond = sqrt( x**2 + y**2) <= RHO_min
	
	#print(sqrt( x**2 + y**2))
	
	return cond	
		
def no_conditional(pos,r0,COEF):
	return True			
	
######################################################





