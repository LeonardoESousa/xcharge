import numpy as np
import random
from kmc.kmc_classes import *
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
    #selection: list of site's index to selected to create particles 
    
    num       = param[0] #relabeling the input parameters 
    selection = param[1]

    part_sel = random.sample(selection,num)
    Ss = [Electron(number) for number in part_sel]
    
    return Ss
    
def gen_hole(param):
    #num: number of excitons
    #selection: list of site's index to selected to create particles
    
    num       = param[0] #relabeling the input parameters 
    selection = param[1]

    part_sel = random.sample(selection,num)
    Ss = [Hole(number) for number in part_sel]
    
    return Ss

def gen_pair_elechole(param):
    #num: number of excitons
    #selection: list of site's index to selected to create particles 
    
    num       = param[0] #relabeling the input parameters 
    selection = param[1]
    
    part_sel = random.sample(selection,int(2*num))
    part_sel_ele = part_sel[:num]
    part_sel_hol = part_sel[num:]
    
    Ss_ele = [Electron(number) for number in part_sel_ele]
    Ss_hol = [Hole(number) for number in part_sel_hol]
    
    Ss = Ss_ele + Ss_hol
    return Ss
    
def gen_excitons(param):
    #num: number of excitons
    #selection: list of site's index to selected to create particles
    
    num       = param[0] #relabeling the input parameters 
    selection = param[1]
    
    part_sel = random.sample(selection,num)
    Ss       = [Exciton('singlet',number) for number in part_sel]    

    return Ss



##########################################
# FUNC TO GIVE S1 AND TRIPLET ENERGIES	
def homo_lumo(param,mats):

    s1s    = param[0] #relabeling the input parameters 
    t1s    = param[1]
    
    s1 = []
    t1 = []
    for i in mats:
        s1.append(np.random.normal(s1s.get(i)[0],s1s.get(i)[1]))
        t1.append(np.random.normal(t1s.get(i)[0],t1s.get(i)[1]))
    return s1,t1,t1,s1 #here we assume HOMO = t1, LUMO = s1

# FUNC TO GIVE S1 AND TRIPLET ENERGIES	
#def gaussian_energy(s1s,mats):
#    tipo = s1s['level']
#    s1 = []
#    for i in mats:
#        s1.append(np.random.normal(s1s[i][0],s1s[i][1]))
#    return s1, tipo


# FUNC TO GIVE S1 AND TRIPLET ENERGIES	
class Gaussian_energy():
	def __init__(self,s1s):
		self.s1s = s1s

	def assign_energy(self,mats):	
		s1 = []
		tipo = self.s1s['level']
		for i in mats:
			s1.append(np.random.normal(self.s1s[i][0],self.s1s[i][1]))
		return s1, tipo





#energy functions when you have the distribuitions of t1 and s1 energies    
def s1_t1_distr(param,mats):
    #s1_dic = {0: 'dist_s1_mat0.txt, 1:'dist_s1_mat1.txt' ...}
    s1_dic    = param[0] #relabeling the input parameters 
    t1_dic    = param[1]
    
       
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
    return s1,t1,t1,s1 #here we assume HOMO = t1, LUMO = s1
   
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
def anni_sing(system,tipos,Ss,indices,locs):
    duplicates = set([x for x in tipos if tipos.count(x) > 1]) # checking if there is 2 occurrences of the types
    if 'singlet' in duplicates:        
        Ss[indices[0][tipos.index('singlet')]].kill('anni',system,system.s1)
             
        
############################################           
# FUNCS TO SHAPE THE GENERATION OF PARTICLES

#main function, regardless of the conditional, filters the sites according to it	
def filter_selection(X,Y,Z,Mats,param_dic):

	X_cop     = X.copy()
	Y_cop     = Y.copy()
	Z_cop     = Z.copy()
	Mats_cop  = Mats.copy()
	selection = []
	
	shape_dic           = param_dic.get('shape_dic')
	mat_restriction     = param_dic.get('mat')    #list of materials that will generate part
	shape_type          = param_dic.get('shape')  #the shape that the particles will be aranged
	argument_shape_func = param_dic.get('argum')  #arguments of the shape func
	origin              = param_dic.get('origin') #origin of the shape
	conditional         = shape_dic.get(shape_type)
	
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
	
	r_pos = np.sqrt(x**2 + y**2 + z**2)
	
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
	
	cond = np.sqrt( x**2 + y**2) <= RHO_min
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
				
# dictonary of all possible geometric creation of particles				
shape_dic = {'sphere': sphere_conditional, 'plane':plane_conditional,'cone':cone_conditional,'cilinder':cilinder_conditional,'rectangle': rectangle_conditional,'free':no_conditional}				


######################################################
#Lattice funcs

def lattice(param):
	
    num_molecs =  int(param[0]) #adjusting to user's input     
    vector     =  [ float(x) for x in param[1] ]
    ps         =  [ float(x)   for x in param[2] ]
    
    
    X, Y, Z, Mats = [], [], [],[]
    ps = [i/np.sum(ps) for i in ps]
    ps = np.cumsum(ps)
    
    dx = vector[0]
    dy = vector[1]
    dz = vector[2]
    dim = []
    for elem in vector:
        if elem != 0:
            dim.append(1)
        else:
            dim.append(0)
    numx = max(dim[0]*int(num_molecs**(1/np.sum(dim))),1)
    numy = max(dim[1]*int(num_molecs**(1/np.sum(dim))),1)
    numz = max(dim[2]*int(num_molecs**(1/np.sum(dim))),1)
	  
    #Se somar +1 vai aparecer duplicados no 2D
    for nx in range(numx):
        for ny in range(numy):
            for nz in range(numz):
                X.append(nx*dx)
                Y.append(ny*dy)
                Z.append(nz*dz)
                sorte = random.uniform(0,1)
                chosen = np.where(sorte < ps)[0][0]
                Mats.append(chosen)
    X = np.array(X)
    Y = np.array(Y)
    Z = np.array(Z)
    Mats = np.array(Mats)    
    return X,Y,Z,Mats


#generates a bulk heterojunction conformation
def lattice_BHJ(param):
	
	num_molecs =  int(param[0]) #adjusting to user's input     
	vector     =  [ float(x)   for x in param[1] ]
	ps         =  [ float(x)   for x in param[2] ]
	n_loops    =  int(param[3])
	cutoff     =  float(param[4])
	
	
	X, Y, Z, Mats = [], [], [],[]
	ps = [i/np.sum(ps) for i in ps]
	ps = np.cumsum(ps)

	dx = vector[0]
	dy = vector[1]
	dz = vector[2]
	dim = []
	for elem in vector:
		if elem != 0:
			dim.append(1)
		else:
			dim.append(0)
            
	numx = max(dim[0]*int(num_molecs**(1/np.sum(dim))),1)
	numy = max(dim[1]*int(num_molecs**(1/np.sum(dim))),1)
	numz = max(dim[2]*int(num_molecs**(1/np.sum(dim))),1)
	 
	#Se somar +1 vai aparecer duplicados no 2D
	for nx in range(numx):
		for ny in range(numy):
			for nz in range(numz):
			
				X.append(nx*dx)
				Y.append(ny*dy)
				Z.append(nz*dz)
				sorte = random.uniform(0,1)
				chosen = np.where(sorte < ps)[0][0]
				Mats.append(chosen)
	X = np.array(X)
	Y = np.array(Y)
	Z = np.array(Z)
	Mats = np.array(Mats)
	
	
	#finding the neighbors around each site based on a cutoff distance
	neighbors_lattice = []
	for j in range(len(X)):
		
		r = [X[j],Y[j],Z[j]]
		neighbors = filter_mats_by_distance(r,X,Y,Z,Mats,cutoff,j)
		neighbors_lattice.append(neighbors)
	
	
	#interating the BHJ process
	for i in range(n_loops):
		Mats_new = np.copy(Mats)
		for k in range(len(X)): #finding which material is more present around the j-th site
			mats_selected  = [ Mats[ind] for ind in neighbors_lattice[k] ] 
			unique, counts = np.unique(mats_selected, return_counts=True)
			
			mat_dict = dict(zip(unique, counts))
			max_mat = max(mat_dict, key=mat_dict.get)
		
			#dealing with mats having equal max occurence
			max_occurence = mat_dict[max_mat]
			max_dict = {}
			for entry in mat_dict:
				if mat_dict[entry] == max_occurence:
					max_dict[entry] = max_occurence
			
			mat_keys = list(max_dict.keys())
			ps = np.zeros([len(mat_keys)])
			ps = ps +1
			ps = [i/np.sum(ps) for i in ps]
			ps = np.cumsum(ps)
			dice = random.uniform(0,1)
			chosen = np.where(dice < ps)[0][0]
			
			
			Mats_new[k] = mat_keys[chosen]
			
		Mats = np.copy(Mats_new)
	Mats = np.copy(Mats_new)

	unique, counts = np.unique(Mats, return_counts=True)
	mat_dict = dict(zip(unique, counts))

	return X,Y,Z,Mats
	
def multiply_lattice(lattice,n_times_ar,delta):
		
	X_lenght = delta[0] #Increments
	Y_lenght = delta[1]
	Z_lenght = delta[2]
	
		
	n_times_x = n_times_ar[0]
	n_times_y = n_times_ar[1]
	n_times_z = n_times_ar[2]	
		
	new_lattice = np.copy(lattice)
	n_sites = len(lattice)
		
	for dx in range(n_times_x):
		for dy in range(n_times_y):
			for dz in range(n_times_z):
					
				new_cell = np.copy(lattice)
				for i in range(n_sites):
					#multipling the X position in the unit cell dx times	
					#beware of the convention X,Y,Z,Mats in the lattice			
					new_cell[i][0] = new_cell[i][0] +dx*X_lenght 
					new_cell[i][1] = new_cell[i][1] +dy*Y_lenght
					new_cell[i][2] = new_cell[i][2] +dz*Z_lenght
					
				new_lattice = np.vstack((new_lattice,new_cell))
		
	return np.unique(new_lattice,axis=0)
def bilayer(par):
	n_times    = int(float(par[0]))
	axis       = par[1]
	param_latt = par[2]
		
    
	#direction in which we will join the two sides
	if axis == 'X':
		vec_joint = [1,0,0]
	if axis == 'Y':
		vec_joint = [0,1,0]
	if axis == 'Z':
		vec_joint = [0,0,1]
	
	X_left ,Y_left ,Z_left ,Mats_left  = lattice(param_latt)
	X_right,Y_right,Z_right,Mats_right = lattice(param_latt)
	
	non_duplicate_mat = set(Mats_left.tolist())
	
	#mat index of the right cell
	max_mat = max([int(mat) for mat in non_duplicate_mat ])	
	max_mat = max_mat + 1

	new_Mats_right = Mats_right.copy()
	n_mats = len(Mats_right)
	new_Mats_right = new_Mats_right + max_mat
	
	Mats_right = new_Mats_right.copy()
		
		
	#left box size	
	dX_left = np.amax(X_left) - np.amin(X_left)
	dY_left = np.amax(Y_left) - np.amin(Y_left)
	dZ_left = np.amax(Z_left) - np.amin(Z_left)
	
	#right box size
	dX_right = np.amax(X_right) - np.amin(X_right)
	dY_right = np.amax(Y_right) - np.amin(Y_right)
	dZ_right = np.amax(Z_right) - np.amin(Z_right)
	
	unit_cell_left  = [[X_left[i],Y_left[i],Z_left[i],Mats_left[i]] for i in range(len(X_left)) ]
	unit_cell_right = [[X_right[i],Y_right[i],Z_right[i],Mats_right[i]] for i in range(len(X_right)) ]
	
	unit_cell_left  = np.array(unit_cell_left)
	unit_cell_right = np.array(unit_cell_right)
	
	n_times_left = [n_times,n_times,n_times]
	multiplied_left  = multiply_lattice(unit_cell_left,n_times_left,[dX_left,dY_left,dZ_left])
	
	mult_x_left = multiplied_left[:,0]
	mult_y_left = multiplied_left[:,1]
	mult_z_left = multiplied_left[:,2]
	
	X_tot_left = np.amax(mult_x_left)-np.amin(mult_x_left)
	Y_tot_left = np.amax(mult_y_left)-np.amin(mult_y_left)
	Z_tot_left = np.amax(mult_z_left)-np.amin(mult_z_left)	
	
	#try/except in case the lattice is 1-D or 2-D
	dimension = 3
	try:
		n_times_x_right = int(X_tot_left/dX_right)
	except:
		n_times_x_right = 1
		dimension       = dimension-1
	try:
		n_times_y_right = int(Y_tot_left/dY_right)
	except:
		n_times_y_right = 1
		dimension       = dimension-1		
	try:
		n_times_z_right = int(Z_tot_left/dZ_right)
	except:
		n_times_z_right = 1
		dimension       = dimension-1
	
	n_times_right = [n_times_x_right,n_times_y_right,n_times_z_right]
	multiplied_right = multiply_lattice(unit_cell_right,n_times_right,[dX_right,dY_right,dZ_right])
	

	shift_vec = [ X_tot_left,Y_tot_left,Z_tot_left ]
	shift_vec = [ shift_vec[i]*vec_joint[i] for i in range(len(shift_vec))]
		 
	#safety distance so the sites do not overlap	 
	#lowest length in the system, will be used as reference to define a minimum distance towards the two cells
	dist = np.array([i for i in [dX_left,dY_left,dZ_left ] if i > 0])
	
	r_min = np.amin(dist)/(len(X_left)**(1/dimension))
	#shifting the entire right lattice 
	for site in multiplied_right:
		site[0] = site[0] + shift_vec[0] + vec_joint[0]*r_min
		site[1] = site[1] + shift_vec[1] + vec_joint[1]*r_min
		site[2] = site[2] + shift_vec[2] + vec_joint[2]*r_min
	
		
	hetero_lattice = np.vstack((multiplied_left,multiplied_right))
	X    = hetero_lattice[:,0]
	Y    = hetero_lattice[:,1]
	Z    = hetero_lattice[:,2]
	Mats = hetero_lattice[:,3]
	
	return 	X,Y,Z,Mats
