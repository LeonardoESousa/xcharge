import numpy as np
import random
from kmc.particles import *
import sys
from collections import Counter

# SET OF FUNCS THAT GENERATE THE MORPHOLOGY OF THE SYSTEM
# note: always define the function by list (param) that contains the things needed
#### CHOOSE A FUNC TO GENERATE PARTICLES

def randomized(available_sites, num_sites, system, kwargs):
    """
    Selects a specified number of unique sites randomly from available sites, 
    ensuring that each selected site corresponds to a specified material type.

    Parameters:
    - available_sites (set): A set of indices representing available sites.
    - num_sites (int): The number of sites to select.
    - system (object): An object representing the system, which should have a 'mats' attribute mapping site indices to material types.
    - kwargs (dict): A dictionary of additional arguments. Should contain a key 'mat' mapping to a list of acceptable material types.

    Returns:
    - selected_sites (list): A list of the indices of the selected sites.
    """
    acceptable_materials = kwargs['mat']
    
    # Initially select random sites
    initial_selection = random.sample(available_sites, num_sites)
    selected_sites = {site for site in initial_selection if system.mats[site] in acceptable_materials}

    # Set a cutoff for number of attempts to select sites
    cutoff = 10000 * num_sites * np.sqrt(len(available_sites))

    attempts = 0
    while len(selected_sites) < num_sites and attempts < cutoff:
        attempts += 1
        
        # Select additional sites
        additional_selection = random.sample(available_sites, num_sites - len(selected_sites))
        acceptable_sites = {site for site in additional_selection if system.mats[site] in acceptable_materials}

        # Add new sites to the set of selected sites
        selected_sites.update(acceptable_sites)
    
    if attempts == cutoff:
        raise Exception(f"Could not find {num_sites} suitable sites for particle creation after {cutoff} attempts.")
    
    # Convert the set of selected sites back to a list before returning
    return list(selected_sites)

              
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
        system.set_energies(s1,type_en)

#########################################################################################

##LATTICE CLASSES####################################################################################
class ReadLattice():
    def __init__(self,file):
        self.file = file 

    def assign_to_system(self,system):       
        data = np.loadtxt(self.file,comments='#')
        system.set_morph(data[:,0],data[:,1],data[:,2],data[:,3]) 



class Lattice():
    '''Arguments:
    num_sites:   number of sites in the lattice (integer)
    vector:      3-element list with distance between sites in x,y and z directions in Å (list)
    disorder:    3-element list with standard deviation of distances between sites in x,y and z directions in Å (list)
    composition: n-element list with the proportion of each of the n different materials in the lattice  (list)'''
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

        num  = int(self.num_sites**(1/np.sum(dim)))       
        numx = max(dim[0]*num,1)
        numy = max(dim[1]*num,1)
        numz = max(dim[2]*num,1)
        Numx = np.array(range(numx))*self.vector[0]
        Numy = np.array(range(numy))*self.vector[1]
        Numz = np.array(range(numz))*self.vector[2]

        total = numx*numy*numz
        X,Y,Z = np.meshgrid(Numx,Numy,Numz,indexing = 'ij')
        X,Y,Z = X.reshape(total,),Y.reshape(total,),Z.reshape(total,)
        # add noise to the lattice
        X = X + np.random.normal(0,self.disorder[0],total)
        Y = Y + np.random.normal(0,self.disorder[1],total)
        Z = Z + np.random.normal(0,self.disorder[2],total)

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
    '''Arguments'''
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
         
class Electric():
    def __init__(self,**kwargs):
        self.eps   = kwargs['eps']
        self.field = kwargs['field']

    def assign_to_system(self,system):
        system.set_medium(self.eps)    
        system.set_electric_field(self.field)
