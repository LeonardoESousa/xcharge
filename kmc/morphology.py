import numpy as np
import random
import matplotlib.pyplot as plt
from kmc.kmc_classes import *
from math import sqrt
# SET OF FUNCS THAT GENERATE THE MORPHOLOGY OF THE SYSTEM
# note: always define the function by list (param) that contains the things needed


# reads the morphology
def read_lattice(file_name):
	X,Y,Z,Mats = [], [], [], []
	
	with open(file_name, 'r') as f:
		for line in f:
			if '#' not in line:		
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

    part_sel = random.sample(selection,num)
    Ss = [Electron(number) for number in part_sel]
    
    return Ss
    
def gen_hole(param):
    #num: number of excitons
    #X: one of the position vectors 
    
    num       = param[0] #relabeling the input parameters 
    selection = param[1]

    part_sel = random.sample(selection,num)
    Ss = [Hole(number) for number in part_sel]
    
    return Ss
def gen_pair_elechole(param):
    #num: number of excitons
    #X: one of the position vectors 
    
    num       = param[0] #relabeling the input parameters 
    selection = param[1]
    
    part_sel = random.sample(selection,int(2*num))
    part_sel_ele = part_sel[0:num]
    part_sel_hol = part_sel[num+1:2*num]
    
    Ss_ele = [Electron(number) for number in part_sel_ele]
    Ss_hol = [Hole(number) for number in part_sel_hol]
    
    Ss = Ss_ele + Ss_hol
    return Ss
    
def gen_excitons(param):
    num       = param[0] #relabeling the input parameters 
    selection = param[1]
    
    part_sel = random.sample(selection,num)
    Ss = [Exciton('singlet',number) for number in part_sel]    

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
   
   
#energy functions when you have the distribuitions of t1 and s1 energies    
def s1_t1_distr(param):
    #s1_dic = {0: 'dist_s1_mat0.txt, 1:'dist_s1_mat1.txt' ...}
    s1_dic    = param[0] #relabeling the input parameters 
    t1_dic    = param[1]
    mats      = param[2]
    

    
    s1s = {}
    t1s = {}
    
    for ele in s1_dic:
    	s1s[ele] = np.loadtxt(s1_dic.get(ele))
    	t1s[ele] = np.loadtxt(t1_dic.get(ele))
    
    s1 = []
    t1 = []
    
    for i in mats:
        s1.append(random.choice(s1s.get(i)))
        t1.append(random.choice(t1s.get(i)))   
    return s1,t1
   
def homo_lumo_s1_distr(param):
    #s1_dic = {0: 'dist_s1_mat0.txt, 1:'dist_s1_mat1.txt' ...}
    homo_dic    = param[0] #relabeling the input parameters
    lumo_dic    = param[1]
    s1_dic      = param[2]
    mats        = param[3]
    
    s1s = {}
    homos = {}
    lumos = {}
    
    for ele in s1_dic:
    	s1s[ele]   = np.loadtxt(s1_dic.get(ele))
    	homos[ele] = np.loadtxt(homo_dic.get(ele))
    	lumos[ele] = np.loadtxt(lumo_dic.get(ele))

    
    homo_dist,lumo_dist,s1_dist = [], [], []
    for i in mats:
    
        s1_dist.append(random.choice(s1s.get(i)))
        homo_dist.append(random.choice(homos.get(i)))
        lumo_dist.append(random.choice(lumos.get(i)))
    return s1_dist, homo_dist, lumo_dist  
   
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
				#print(x,y,z)
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
	
def rectangle_conditional(pos,r0,COEF):
	
	x  = pos[0]
	y  = pos[1]
	z  = pos[2]	
	
	xlim = COEF[0]
	ylim = COEF[1]
	zlim = COEF[2]
	
	x_i = xlim[0]
	x_f = xlim[1]
	
	y_i = ylim[0]
	y_f = ylim[1]
	
	z_i = zlim[0]
	z_f = zlim[1]
	
	if ( ( x >= x_i) and ( x <= x_f)):
		if ( ( y >= y_i) and ( y <= y_f)):
			if ( ( z >= z_i) and ( z <= z_f)):
				return True
######################################################





