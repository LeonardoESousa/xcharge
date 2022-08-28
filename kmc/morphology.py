import numpy as np
import random
from kmc.particles import *
import sys
from collections import Counter

# Personalized errors
class randomized_error(Exception):
    
    def __init__(self,number):
        self.message="We were not able to include the " + str(number) + " particle(s) into the lattice. Loosen up the spacial restrictions for particle creation!"
        super().__init__(self.message)
    
    
    

# SET OF FUNCS THAT GENERATE THE MORPHOLOGY OF THE SYSTEM
# note: always define the function by list (param) that contains the things needed
#### CHOOSE A FUNC TO GENERATE PARTICLES

'''
def randomized(available, number, system, kwargs):
    mat = kwargs['mat']
    selected = random.sample(list(available),number)
    selected = [s for s in selected if system.mats[s] in mat]    
    count   = 0
    cutoff  = number*10    
    while len(selected) < number and count < cutoff:
        count = count + 1
        new = random.sample(list(available),number-len(selected))
        for n in new:
            if system.mats[n] in mat:
                selected.append(n)                      
    if count >= cutoff and len(selected) != number:
        raise randomized_error(number)                   
    return selected            
'''
    
def randomized(available, number, system, kwargs):
    mat = kwargs['mat']
    selected = random.sample(list(available),number)
    selected = [s for s in selected if system.mats[s] in mat]
    selec_size_old = len(selected)
    selec_size_new = number # could be anything that satisfies selec_size_old != selec_size_new = True, but this way is granted to work if  len(selected) < number = True
    
    count = 1
    cutoff = 10000*number*np.sqrt(len(list(available)))
    #print(count,number,selec_size_old,selec_size_new)
    while len(selected) < number and count < cutoff:
    
        count = count + 1
        selec_size_old = len(selected)
        sel_old = selected.copy()
        new      = random.sample(list(available),number-len(selected))
        new_mat  = [sel for sel in new if system.mats[sel] in mat] 
        selected = list(set(selected + new_mat)) #enforcing uniqueness                       
        selec_size_new = len(selected)
        '''
        print()
        print('--->',count,selec_size_old,selec_size_new,new, [system.mats[a] for a in new])
        print('---> sel. old',sel_old)
        print('---> sel. new',selected)
        '''
    '''    
    print('#')    
    print(count,selec_size_old,selec_size_new, len(selected) < number,selec_size_old != selec_size_new )
    print(cutoff,count)
    print(selected)
    '''
    
    if count == cutoff:#if we could not find suitable sites for particle creation
        raise randomized_error(number)                   
    return selected          
      
def interface(available, number, system, kwargs):
    mat = kwargs['mat']
    neighbors = kwargs['neigh']
    available = [i for i in available if system.mats[i] in mat]
    new_available = []
    for i in range(len(available)):
        R = np.copy(system.R) 
        dR = R - R[i,:]
        modulo = np.sqrt(np.sum(dR*dR,axis=1))
        indices = np.argsort(modulo)[1:neighbors+1]
        if np.any(system.mats[indices] != system.mats[i]):
            new_available.append(i)
        
    selected = random.sample(new_available,number)        
    return selected            

##CLASS FOR GENERATING PARTICLES IN THE SYSTEM###########################################
class Create_Particles():
    def __init__(self,kind, num, method, **kwargs):
        self.kind   = kind
        self.num    = num
        self.method = method
        self.argv   = kwargs

    def assign_to_system(self,system):
        selected = self.method(range(len(system.X)),self.num, system, self.argv)
        Particula = getattr(sys.modules[__name__], self.kind.title())
        particles = [Particula(number) for number in selected]
        system.set_particles(particles)
        
class Create_Particles_PROB():
    def __init__(self, probs, num, method, **kwargs):
        self.num    = num
        self.method = method
        self.probs  = probs
        self.argv   = kwargs
    def assign_to_system(self,system):
        coin      = [random.uniform(0,1) for num in range(self.num)]
        type_part = [part[0] for part in self.probs]
        prob      = [part[1] for part in self.probs]
        cum_sum   = np.cumsum(prob/np.sum(prob))
        x         = [ type_part[np.where(random.uniform(0,1) <= cum_sum)[0][0]] for n in range(self.num)]
        pop       = Counter(x)       
        for part in pop:		
            selected = self.method(range(len(system.X)),pop[part], system, self.argv)
            Particula = getattr(sys.modules[__name__], part)
            particles = [Particula(number) for number in selected]
            system.set_particles(particles)
#########################################################################################

##CLASSES TO ASSIGN ENERGIES TO LATTICE##################################################	
class Gaussian_energy():
    def __init__(self,s1s):
        self.s1s = s1s

    def assign_to_system(self,system):	
        type_en = self.s1s['level']
        uniq = system.uniq
        N = len(system.mats)
        
        means = np.empty(N)
        stds  = np.empty(N)
        means.fill(self.s1s[uniq[0]][0])
        stds.fill(self.s1s[uniq[0]][1])
        for m in uniq[1:]:    
            mask = system.mats == m
            means[mask] = self.s1s[m][0]
            stds[mask]  = self.s1s[m][1]
        s1 = np.random.normal(means,stds,N) 

        #if type(self.s1s[uniq[0]]) == str: #if the user provides a file that contains the energies
        #    s1 = np.random.choice(np.loadtxt(self.s1s[uniq[0]]),size = N)
        #else: #if you want to generate an on-the-fly gaussian distr.
        #    s1 = np.random.normal(self.s1s[uniq[0]][0],self.s1s[uniq[0]][1],N)
        #materials = [i for i in uniq if i != uniq[0]]
        #for m in materials:
        #    print('aqui')
        #    if type(self.s1s[m]) == str:
        #        s11 = np.random.choice(np.loadtxt(self.s1s[m]),size = N)
        #    else:
        #        s11 = np.random.normal(self.s1s[m][0],self.s1s[m][1],N)
        #    s1[system.mats == m] = s11[system.mats == m]
        #print(s1)             
        system.set_energies(s1,type_en)

#########################################################################################

# FUNCS TO SHAPE THE GENERATION OF PARTICLES


#funcs that define the filtering, like inside a sphere, plane etc
class sphere_conditional():
	def __init__(self):
		self.args = ['origin','radius_min']
	def test_cond(self,pos,args):
		r0      = args.get('origin')
		r_min   = args.get('radius_min')
		x  = pos[0] - r0[0]
		y  = pos[1] - r0[1]
		z  = pos[2] - r0[2]
		r_pos = np.sqrt(x**2 + y**2 + z**2)
		return r_pos <= r_min #condition to be inside the sphere
    
class plane_conditional():
	def __init__(self):
		self.args = ['origin','coefs']
	def test_cond(self,pos,args):
		r0      = args.get('origin')
		coefs   = args.get('coefs')		
		#Equation of plane is defined as 
		# A(X-X0) + B(Y-Y0) + C(Z-Z0) = 0
		x  = pos[0] - r0[0]
		y  = pos[1] - r0[1]
		z  = pos[2] - r0[2]

		A,B,C = coefs
		cond = A*x + B*y + C*z == 0 #condition to be in the plane
		
		return cond
    
class cone_conditional(): #hollow cone
	def __init__(self):
		self.args = ['origin','coefs']
	def test_cond(self,pos,args):
		r0      = args.get('origin')
		coefs   = args.get('coefs')
			
		#Equation of cone is defined as 
		# (X-X0)^2/A^2 + (Y-Y0)^2/B^2 + (Z-Z0)^2/C^2 = 0
		x  = pos[0] - r0[0]
		y  = pos[1] - r0[1]
		z  = pos[2] - r0[2]

		A,B,C = coefs

		    
		cond = (x**2)/(A**2) + (y**2)/(B**2) -(z**2)/(C**2) == 0 #condition to be in the cone
		    
		return cond	
    
class cilinder_conditional():
	def __init__(self):
		self.args = ['origin','radius_min']
	def test_cond(self,pos,args):
	
		r0      = args.get('origin')
		RHO_min = args.get('radius_min')
		
		x  = pos[0] - r0[0]
		y  = pos[1] - r0[1]
		z  = pos[2] - r0[2]
    
		cond = np.sqrt( x**2 + y**2) <= RHO_min
		return cond	
        
class no_conditional():
	def __init__(self):
		pass
	def test_cond(self,pos,args):
		return True			
		                
class rectangle_conditional():
	def __init__(self):
		self.args = ['limits']	
	def test_cond(self,pos,args):
		limits = args.get('limits')
		
		x,y,z          = pos	
		xlim,ylim,zlim = limits
		
		x_i,x_f = xlim
		y_i,y_f = ylim
		z_i,z_f = zlim
		
		
		if (( x >= x_i) and ( x <= x_f)):
			if (( y >= y_i) and ( y <= y_f)):
				if (( z >= z_i) and ( z <= z_f)):
					return True 
                
#main function, regardless of the conditional, filters the sites according to it	
def conditional_selection(available, number, system, kwargs):
	# dictonary of all possible geometric creation of particles				
	type_dic = {'sphere': sphere_conditional, 'plane':plane_conditional,'cone':cone_conditional,
	'cilinder':cilinder_conditional,'rectangle': rectangle_conditional,'free':no_conditional}

	selected      = []
	mat_restric   = kwargs['mat']
	shape         = kwargs['type_cond']
	conditional   = type_dic.get(shape)()
	
	
	try:#getting all kwargs given by the user
		list_kwargs   = conditional.args
		dict_kwargs   = {}
		for arg in list_kwargs:
			dict_kwargs[arg] = kwargs[arg]
		
	except:#in case that the conditional does no require additiona kwargs
		dict_kwargs = []	
	
	
	X_cop     = system.X.copy()
	Y_cop     = system.Y.copy()
	Z_cop     = system.Z.copy()
	Mats_cop  = system.mats.copy()	
	
	for index_mol in range(len(X_cop)):

		x   = X_cop[index_mol]
		y   = Y_cop[index_mol]
		z   = Z_cop[index_mol]
		mat = Mats_cop[index_mol]
		
		if conditional.test_cond([x,y,z],dict_kwargs) and (mat in mat_restric):
			selected.append(index_mol)
			
	return random.sample(selected,number)                
				


##LATTICE CLASSES####################################################################################
class ReadLattice():
    def __init__(self,file):
        self.file = file 

    def assign_to_system(self,system):       
        data = np.loadtxt(self.file,comments='#')
        system.set_morph(data[:,0],data[:,1],data[:,2],data[:,3]) 



class Lattice():
    def __init__(self,num_sites,vector,disorder,composition): #initializing the lattice class with some basic info given by the user
        self.num_sites   = num_sites
        self.vector      = vector
        self.disorder    = disorder
        self.composition = np.cumsum([i/np.sum(composition) for i in composition])
        
    def make(self): #Generating the set X,Y,Z,Mats
        dim = []
        for elem in self.vector:
            if elem != 0:
                dim.append(1)
            else:
                dim.append(0)
                
        numx = max(dim[0]*int(self.num_sites**(1/np.sum(dim))),1)
        numy = max(dim[1]*int(self.num_sites**(1/np.sum(dim))),1)
        numz = max(dim[2]*int(self.num_sites**(1/np.sum(dim))),1)

        Numx = np.random.normal(np.array(range(numx))*self.vector[0],self.disorder[0],numx) 
        Numy = np.random.normal(np.array(range(numy))*self.vector[1],self.disorder[1],numy) 
        Numz = np.random.normal(np.array(range(numz))*self.vector[2],self.disorder[2],numz) 
        
        total = numx*numy*numz
        X,Y,Z = np.meshgrid(Numx,Numy,Numz,indexing = 'ij')
        X,Y,Z = X.reshape(total,),Y.reshape(total,),Z.reshape(total,)
        
        luck = np.random.uniform(0,1,total)
        Mats = np.zeros(total)
        for i in reversed(range(len(self.composition))):
            Mats[ luck < self.composition[i] ] = i
        return X,Y,Z,Mats     
    
        
    def assign_to_system(self, system): #adding the X,Y,Z,Mats to the system
        X, Y, Z, Mats = self.make()
        system.set_morph(X,Y,Z,Mats)
        
#list sites' indexes of all neighbors of one given position
def filter_mats_by_distance(r,X,Y,Z,Mats,cutoff,r_index):
    x = r[0]
    y = r[1]
    z = r[2]
    neighbors = []
    
    for n in range(len(X)):
    
        dx   = (x - X[n])
        dy   = (y - Y[n])
        dz   = (z - Z[n])
        dist = np.sqrt(dx**2+dy**2+dz**2)
        if(dist <= cutoff and r_index != n):
            neighbors.append(n)
    #return np.array(neighbors)
    return neighbors        

class Lattice_BHJ():
    def __init__(self,n_loops,cutoff,num_sites,vector,disorder,composition):
        self.n_loops     = n_loops
        self.cutoff      = cutoff
        self.num_sites   = num_sites
        self.vector      = vector
        self.disorder    = disorder
        self.composition = composition #np.cumsum([i/np.sum(composition) for i in composition])
        self.lattice     = Lattice(self.num_sites,self.vector,self.disorder,self.composition)
    
    def make(self):
        X,Y,Z,Mats = self.lattice.make()
        
        #finding the neighbors around each site based on a cutoff distance
        neighbors_lattice = []
        for j in range(len(X)):
        
            r = [X[j],Y[j],Z[j]]
            neighbors = filter_mats_by_distance(r,X,Y,Z,Mats,self.cutoff,j)
            neighbors_lattice.append(neighbors)
    
    
        #interating the BHJ process
        for i in range(self.n_loops):
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
                    
    def assign_to_system(self, system): #adding the X,Y,Z,Mats to the system
        X, Y, Z, Mats = self.make()
        system.set_morph(X,Y,Z,Mats)
         
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
    
    
class Bilayer():    
    def __init__(self,axis,num_sites,vector,disorder,composition):
        self.axis        = axis       
        self.num_sites   = num_sites
        self.vector      = vector
        self.disorder    = disorder
        self.composition = composition #np.cumsum([i/np.sum(composition) for i in composition])
        self.lattice     = Lattice(self.num_sites,self.vector,self.disorder,self.composition)
        
    def make(self):    
        #direction in which we will join the two sides
        if self.axis == 'X':
            vec_joint = [1,0,0]
        if self.axis == 'Y':
            vec_joint = [0,1,0]
        if self.axis == 'Z':
            vec_joint = [0,0,1]
    
        
        X_left ,Y_left ,Z_left ,Mats_left  = self.lattice.make()
        X_right,Y_right,Z_right,Mats_right = self.lattice.make()   
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
        
        n_times_left = [1,1,1] #number of times that the lattice will be duplicated. Setted 1 to avoid overlap
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
         
    def assign_to_system(self, system): #adding the X,Y,Z,Mats to the system
        X, Y, Z, Mats = self.make()
        system.set_morph(X,Y,Z,Mats)

#checks if a site is in the list_set        
def is_in(list_set,n):    
    for set_indx in list_set:
        x = set_indx[0]
        y = set_indx[1]
        #print(n,x,y,n >= x,n <=y)
        if n >= x and n <=y:
            return True   		
    return False
#this function adds, to a given lattice, teeths, of size_saw size and width of width_saw   
def make_teeth(X,Y,Z,Mats,size_saw,width_saw,ps,list_set,dr,nums,disorder):
    dx,dy,dz = dr
    numx,numy,numz = nums
    
    if dx != 0:
        numx_f = numx + size_saw
        numx_i = numx 

    if dy != 0:
        numy_f = numy
        numy_i = 0 
    else:
        numy_f = 1
        numy_i = 0    		
    if dz != 0:
        numz_f = numx
        numz_i = 0 
    else:
        numz_f = 1
        numz_i = 0     

    ps = np.cumsum([i/np.sum(ps) for i in ps])
    for nnx in range(numx_i,numx_f):
        for nny in range(numy_i,numy_f):
            for nnz in range(numz_i,numz_f):    	
                if is_in(list_set,nny):
                    pass
                else:           
                    continue 
                X.append(np.random.normal(nnx*dx,disorder[0]))
                Y.append(np.random.normal(nny*dy,disorder[1]))
                Z.append(np.random.normal(nnz*dz,disorder[2]))
                sorte = random.uniform(0,1)
                chosen = np.where(sorte < ps)[0][0]
                Mats.append(chosen)             
class Saw():    
    def __init__(self,size_saw,width_saw,num_sites,vector,disorder,composition):

        self.num_sites   = num_sites
        self.vector      = vector
        self.disorder    = disorder
        self.size        = size_saw
        self.width       = width_saw
        
        
        
        self.composition_right = composition 
        self.composition_left  = [ 0 for i in range(len(composition))] + composition  
        self.lattice_right = Lattice(self.num_sites,self.vector,self.disorder,self.composition_right)        
        self.lattice_left  = Lattice(self.num_sites,self.vector,self.disorder,self.composition_left)                
        
    def make(self):    
        X_left ,Y_left ,Z_left ,Mats_left  = self.lattice_right.make()
        X_right,Y_right,Z_right,Mats_right = self.lattice_left.make()   
        
        X_left ,Y_left ,Z_left ,Mats_left    = X_left.tolist() ,Y_left.tolist() ,Z_left.tolist() ,Mats_left.tolist()
        X_right ,Y_right,Z_right ,Mats_right = X_right.tolist() ,Y_right.tolist() ,Z_right.tolist() ,Mats_right.tolist()
                
        dim = []
        for elem in self.vector:
            if elem != 0:
                dim.append(1)
            else:
                dim.append(0)
                        
        numx = max(dim[0]*int(self.num_sites**(1/np.sum(dim))),1)
        numy = max(dim[1]*int(self.num_sites**(1/np.sum(dim))),1)
        numz = max(dim[2]*int(self.num_sites**(1/np.sum(dim))),1)  
        nums = [numx,numy,numz]      
        X_max_noteeth = max(X_left)#lattice size before teeth
        
        #getting the sets of the teeths
        list_set = [[self.width*x,self.width*x+self.width-1] for x in range(int(numy/self.width)) if x%2==0] #left
        make_teeth(X_left,Y_left,Z_left,Mats_left,self.size,self.width,
        self.composition_left,list_set,self.vector,nums,self.disorder)
        X_max_teeth = max(X_left)#lattice size after teeth

        X_left = np.array(X_left)
        Y_left = np.array(Y_left)
        Z_left = np.array(Z_left)
        Mats_left = np.array(Mats_left)
        cell_left = [[X_left[i],Y_left[i],Z_left[i],Mats_left[i]] for i in range(len(X_left))]        
        
        #getting the sets of the teeths             
        list_set = [[self.width*x+self.width,self.width*x+self.width-1+self.width] for x in range(int(numy/self.width)) if x%2==0]
        make_teeth(X_right,Y_right,Z_right,Mats_right,self.size,self.width,
        self.composition_right,list_set,self.vector,nums,self.disorder)
        
        #estimating a cutoff distance for the joint
        dX_left = np.amax(X_left) - np.amin(X_left)
        dY_left = np.amax(Y_left) - np.amin(Y_left)
        dZ_left = np.amax(Z_left) - np.amin(Z_left)        
        dist = np.array([i for i in [dX_left,dY_left,dZ_left ] if i > 0])        
        r_min = np.amin(dist)/(len(X_left)**(1/sum(dim)))
        if r_min < 5:
            r_min = 5	    
        X_right = np.array(X_right)
        Y_right = np.array(Y_right)
        Z_right = np.array(Z_right)
        Mats_right = np.array(Mats_right)

	#joining the sides
        X_max      = np.amax(X_right)
        X_right    = -X_right +2*X_max -(X_max_teeth-X_max_noteeth) + r_min
        cell_right = [[X_right[i],Y_right[i],Z_right[i],Mats_right[i]] for i in range(len(X_right)) ]
        cell       = np.vstack((cell_left,cell_right))           
           
        return 	cell[:,0],cell[:,1],cell[:,2],cell[:,3]         
         
    def assign_to_system(self, system): #adding the X,Y,Z,Mats to the system
        X, Y, Z, Mats = self.make()
        system.set_morph(X,Y,Z,Mats)  

class Electric():
    def __init__(self,**kwargs):
        self.eps   = kwargs['eps']
        self.field = kwargs['field']

    def assign_to_system(self,system):
        system.set_medium(self.eps)    
        system.set_electric_field(self.field)
