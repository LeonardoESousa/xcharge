import numpy as np
import random
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d
#plt.switch_backend('agg') #para cluster
import collections
from shutil import copyfile
from math import sqrt
################################
#AUXILIARY CODE TO GENERATE LATTICE FOR KMC

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


#MORPHOLOGY FUNCS (option 2)

#num_molecs: number of molecules
#vector: 3 component lattice vector. For lower dimension, make distance 0 
#ps: vector with the relative probability of each material
def lattice(param):
	
    num_molecs =  int(param[0][0]) #adjusting to user's input     
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
                #print(nz)
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
	
	num_molecs =  int(param[0][0]) #adjusting to user's input     
	vector     =  [ float(x) for x in param[1] ]
	ps         =  [ float(x)   for x in param[2] ]
	n_loops    =  int(param[3][0])
	cutoff     =  float(param[4][0])
	
	
	X, Y, Z, Mats = [], [], [],[]
	ps = [i/np.sum(ps) for i in ps]
	ps = np.cumsum(ps)
	print(ps)
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
			#print(nz)
			
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
		#print(Mats)
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
			
			#print(max_dict,mat_keys,dice)
			#print(ps)
			#print( mat_keys[chosen])
			#print()
			#input()			
			
		
		Mats = np.copy(Mats_new)
		#print("Loop %s out %s!" % (i+1,n_loops))
	Mats = np.copy(Mats_new)

	#input()
	
	unique, counts = np.unique(Mats, return_counts=True)
	mat_dict = dict(zip(unique, counts))
	print()
	print("Done! the new lattice has the following population:")
	print(mat_dict)
	return X,Y,Z,Mats
	
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
'''	
def lattice_BHJ(param):
	
	num_molecs =  int(param[0][0]) #adjusting to user's input     
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
			#print(nz)
			
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
	
	n_loops = 10
	cutoff = 1
	
	unique, counts = np.unique(Mats, return_counts=True)
	mat_dict = dict(zip(unique, counts))
	print(mat_dict,max(mat_dict), mat_dict.get(max(mat_dict)))
	
	
	for i in range(n_loops):
	
		Mats_new = Mats.copy()
		for j in range(len(X)):
		
			r_index = j
			r = [X[r_index],Y[r_index],Z[r_index]]
			mat_local = Mats[r_index]
			
			mats_selected = filter_mats_by_distance(r,X,Y,Z,Mats,cutoff,r_index)
			unique, counts = np.unique(mats_selected, return_counts=True)
			
			mat_dict = dict(zip(unique, counts))

			max_mat = max(mat_dict, key=mat_dict.get)
			
			
			Mats_new[j] = max_mat
		Mats = Mats_new.copy()
	Mats = Mats_new
	
	return X,Y,Z,Mats
	
def filter_mats_by_distance(r,X,Y,Z,Mats,cutoff,r_index):
	x = r[0]
	y = r[1]
	z = r[2]
	mats_selected = []
	R = np.zeros([len(X)])
	for n in range(len(X)):
	
		dx   = (x - X[n])
		dy   = (y - Y[n])
		dz   = (z - Z[n])
		dist = np.sqrt(dx**2+dy**2+dz**2)
		if(dist <= cutoff and r_index != n):
			mats_selected.append(Mats[n])
	return np.array(mats_selected)
'''			
def load_cif(mol_file):

	#This program reads .cif files through a .mol2 files

	print("------------------------------------------------------------------------------------------------------------")
	print("༼つಠ益ಠ༽つ                          LOAD .CIF PROGRAM (via .mol extension)                       ༼つಠ益ಠ༽つ")
	print("author: Tiago Cassiano, version: 1.2v (hoping that will sufice)")
	print("STEPS TO USE THE PROGRAM:")
	print("1) Pick a .cif geometry somewhere with your desired unit cell")
	print("2) Install Mercury software (https://www.ccdc.cam.ac.uk/support-and-resources/Downloads/)")
	print("3) Load .cif there")
	print("		3.1) On the display box (left-side bellow), check packing")
	print("		3.2) Check if the geometry seems like it should be! (Special care if there is hydrogen bonds)")
	print("		3.3) File --> Save as")
	print("		3.4) Choose Mol2 files (*.mol2 *.mol) ")
	print("4) Insert the file name for this program")
	print("(BEWARE) If your structure has disordered atoms, you must remove them. To do so, follow: ")
	print("		A) Install ChimeraX (https://www.rbvi.ucsf.edu/chimerax/) or Avogadro ")
	print("		B) Upload .mol2 file ")
	print("		C) Remove the residues (in chimerax go: Select -> Residues -> choose  then Actions -> bonds -> delete)")
	print()
	print("Not so confident? It is possible to obtain the centers of mass via avogadro. For that, do:")
	print("		A) Load .mol2 ")
	print("		B) Selection Mode: Molecule ")
	print("		C) Select a molecule, then click add center of mass ")
	print("		D) Delete de molecule ")
	print("------------------------------------------------------------------------------------------------------------")
	print()
	print()
	#input_file= input("Assuming that you already finished all these steps, give to me the .mol file name:  ")

	input_file   = mol_file[0][0]

	output_file  = input_file.split(".")[0]+".txt"
	matinfo_file = "matinfo.txt"

	periodic_table = {
	    "H"  :  1.00797,
	    "He" : 4.00260,
	    "Li" : 6.941,
	    "Be" : 9.01218,
	    "B"  : 10.81,
	    "C"  : 12.011,
	    "N"  : 14.0067,
	    "O"  : 15.9994,
	    "F"  : 18.998403,
	    "Ne" : 20.179,
	    "Na" : 22.98977,
	    "Mg" :     24.305,
	    "Al" :	26.98154,
	    "Si" :	28.0855,
	    "P"  :   30.97376,
	    "S"  :   32.06,
	    "Cl" :	35.453,
	    "K"  :   39.0983,
	    "Ar" :	39.948,
	    "Ca" :	40.08,
	    "Sc" :	44.9559,
	    "Ti" :	47.90,
	    "V"  :   50.9415,
	    "Cr" :	51.996,
	    "Mn" :	54.9380,
	    "Fe" :	55.847,
	    "Ni" :	58.70,
	    "Co" :	58.9332,
	    "Cu" :	63.546,
	    "Zn" :	65.38,
	    "Ga" :	69.72,
	    "Ge" :	72.59,
	    "As" :	74.9216,
	    "Se" :	78.96,
	    "Br" :	79.904,
	    "Kr" :	83.80,
	    "Rb" :	85.4678,
	    "Sr" :	87.62,
	    "Y"  :   88.9059,
	    "Zr" :	91.22,
	    "Nb" :	92.9064,
	    "Mo" :	95.94,
	    "Tc" :	98,
	    "Ru" :	101.07,
	    "Rh" :	102.9055,
	    "Pd" :	106.4,
	    "Ag" :	107.868,
	    "Cd" :	112.41,
	    "In" :	114.82,
	    "Sn" :	118.69,
	    "Sb" :	121.75,
	    "I"  :   126.9045,
	    "Te" :	127.60,
	    "Xe" :	131.30,
	    "Cs" :	132.9054,
	    "Ba" :	137.33,
	    "La" :	138.9055,
	    "Ce" :	140.12,
	    "Pr" :	140.9077,
	    "Nd" :	144.24,
	    "Pm" :	145,
	    "Sm" :	150.4,
	    "Eu" :	151.96,
	    "Gd" :	157.25,
	    "Tb" :	158.9254,
	    "Dy" :	162.50,
	    "Ho" :	164.9304,
	    "Er" :	167.26,
	    "Tm" :	168.9342,
	    "Yb" :	173.04,
	    "Lu" :	174.967,
	    "Hf" :	178.49,
	    "Ta" :	180.9479,
	    "W"  :   183.85,
	    "Re" :	186.207,
	    "Os" :	190.2,
	    "Ir" :	192.22,
	    "Pt" :	195.09,
	    "Au" :	196.9665,
	    "Hg" :	200.59,
	    "Tl" :	204.37,
	    "Pb" :	207.2,
	    "Bi" :	208.9804,
	    "Po" :	209,
	    "At" :	210,
	    "Rn" :	222,
	    "Fr" :	223,
	    "Ra" :	226.0254,
	    "Ac" :	227.0278,
	    "Pa" :	231.0359,
	    "Th" :	232.0381,
	    "Np" :	237.0482,
	    "U"  :   238.029,
	    "Pu" :	242,
	    "Am" :	243,
	    "Bk" :	247,
	    "Cm" :	247,
	    "No" :	250,
	    "Cf" :	251,
	    "Es" :	252,
	    "Hs" :	255,
	    "Mt" :	256,
	    "Fm" :	257,
	    "Md" :	258,
	    "Lr" :	260,
	    "Rf" :	261,
	    "Bh" :	262,
	    "Db" :	262,
	    "Sg" :	263,
	    "Uun": 	269,
	    "Uuu": 	272,
	    "Uub": 	277, 
	}



	class atom_class: # atom and its attributes
		
		def __init__(self,element,index,position,neighbors):
			self.element   = element
			self.index     = index
			self.position  = position
			self.neighbors = neighbors		

		def get_mass(self,periodic_table):
			for atom_ele in periodic_table:
				if ( atom_ele.lower() == self.element.lower() ):
					self.mass = periodic_table[atom_ele]
			
	#READING .MOL FILES	
	init_atoms_array     = [] #atoms array without object initialization and neighbors
	init_neighbors_array = []
	atom_part=False
	bond_part=False
	atom_keyword = "@<TRIPOS>ATOM"
	bond_keyword = "@<TRIPOS>BOND" 

	with open(input_file) as f:
		i = 0
		for line in f:
		
			if(atom_part and (atom_keyword not in line ) and (bond_keyword not in line)  ):
				line    = line.split()
				indx    = int(line[0])
				element = line[5].split(".")[0]
				
				x = float(line[2])
				y = float(line[3])
				z = float(line[4])
				pos = [x,y,z]
				init_atoms_array.append([indx,element,pos])
				
			
			if(bond_part and (atom_keyword not in line ) and (bond_keyword not in line)  ):
			
				if("@" in line): #avoid additional data like @<TRIPOS>SUBSTRUCTURE
					break
					
					
				line        = line.split()
				origin_atom = int(line[1])
				target_atom = int(line[2]) # the atom indx that bounds with the origin_atom
				init_neighbors_array.append([origin_atom,target_atom])
				
			if( atom_keyword in line):
				atom_part=True
			if( bond_keyword in line):
				bond_part=True		
				atom_part=False

	 # init_neighbors_array but with a list of target atoms in the same line 
	 # [ [1,2], [1,3] ... --> [ [1,[2,3]] ...  
	init_neighbors_array_joint = []
	N_at = len(init_atoms_array)

	for i in range(1,N_at+1):

		neighbor_ith = []
		for bond in init_neighbors_array:
			origin_atom = bond[0]
			target_atom = bond[1]
			if( origin_atom == i ):
				neighbor_ith.append(target_atom)
			if( target_atom == i ):
				neighbor_ith.append(origin_atom)
			
		init_neighbors_array_joint.append([i,neighbor_ith])

	atoms_array     = []	
	for i in range(N_at):	
		indx_initneigh = init_neighbors_array_joint[i][0]
		indx_initatoms = init_atoms_array[i][0]
		
		
		if(indx_initneigh == indx_initatoms):
		
			indx     = indx_initneigh # atom's index
			atom_ele = init_atoms_array[i][1]
			pos      = init_atoms_array[i][2] # [x,y,z]
			neighbor = init_neighbors_array_joint[i][1]
				

			if not neighbor: #workaround for cases where the molecue is a single atom
				neighbor = [indx]
			
	
			atoms_array.append( atom_class(atom_ele,indx,pos,neighbor) )
				 
		
	#input()	
		
		
	for atom_individual in atoms_array: #calculating the mass of each atom
		atom_individual.get_mass(periodic_table)
		




	# Finding the molecules in the crystal
	# agragates_array = [ [1,2,3,4] , [6,7,52] ...] => 1,2,3,4 are in the same molecule, 6,7 and 52 consist of another one...
	N_at = len(atoms_array)
	neighbor_arr = []
	for i in range(N_at):
		neighbor_arr.append( [atoms_array[i].index] + atoms_array[i].neighbors )
	neighbor_arr = np.array(neighbor_arr,dtype=object)
	neighbor_arr_new = neighbor_arr
	first_it = True #first iteration



	# the following lines deal with a limitation of python. If the array of neighbors is symmetric
	# (if all atoms have the same number of neighbors). On this situation, the array intialized in:
	# neighbor_arr = np.array(neighbor_arr,dtype=object) 
	# has a restrict shape (n,n_neigh) where n_neigh is the number of neighbors of each atom
	# however, we must have arrays shapping (n,) regardless of the n_neigh
	# to overcome the issue, I included a dirty list to disrupt the symmetry. Then, I removed.

	neigh_len_array = [len(atom.neighbors) for atom in atoms_array]
	n_nla           = len(neigh_len_array)
	nla0            = neigh_len_array[0]
	symetric_case   = True
	for i in range(1,n_nla):
		if (nla0 != neigh_len_array[i]):
			symetric_case =False
			break
	if(symetric_case):
		y_len = np.shape(neighbor_arr)[1]
		dirty = np.arange(y_len+1).tolist()
		neighbor_arr = []
		for i in range(N_at):
			neighbor_arr.append( [atoms_array[i].index] + atoms_array[i].neighbors )
		neighbor_arr.append(dirty)
		neighbor_arr = np.array(neighbor_arr,dtype=object)
		neighbor_arr = np.delete(neighbor_arr,-1)
		neighbor_arr_new = neighbor_arr
	########## END GAMBIARRA ###################################



	#while( ( neighbor_arr != neighbor_arr_new).any() or (first_it == True)):
	while( ( not np.array_equal(neighbor_arr,neighbor_arr_new) ) or (first_it == True)  ):
		first_it = False
		neighbor_arr = neighbor_arr_new	
		
		for j in range(N_at):
			L_j = neighbor_arr[j]
			for k in range(N_at):
				L_k = neighbor_arr[k]

				L_j_set = set(L_j)
				L_k_set = set(L_k)
					
				done=False
				
				#if the neighbor of atom J coincides with K, sum them
				for ele_j in L_j_set:
					for ele_k in L_k_set:
						if (ele_j == ele_k):
							done=True
							break
					if(done):
						break

				if((done)):



					neighbor_arr_new[j] = list(set(L_j + L_k))
					neighbor_arr_new[k] = neighbor_arr_new[j]

				# neighbor_arr_new = neighbor_arr + all first neighbors
				# if neighbor_arr_new == neighbor_arr => reached the molecules
			
	agragates_array = []
	for i in range(N_at):
		L_i = list(set(neighbor_arr[i]))
		if L_i not in agragates_array: #removing duplicates and renaiming neighbor_arr_new
			agragates_array.append(L_i)
		

	agragates_array = np.array(agragates_array, dtype=object)






		
	# changing from an array of indexes to an array of objects
	def get_atoms(agragates_array,atoms_array):
		new_agg_array = []
		N_mol = len(agragates_array)
		
		for i in range(N_mol):
			group_atoms = []
			molecule = agragates_array[i]
			#print(molecule)
			#print()
			for atom in molecule:

				for atom_test in atoms_array:
					if ( atom_test.index == atom ):
						
						group_atoms.append(atom_test)
						break
				aux = [ a.index for a in group_atoms]
				#print(aux)		
				#print(atom, atom_test.index)
				#input()
			new_agg_array.append(group_atoms)
			
		return new_agg_array
	agragates_array = get_atoms(agragates_array,atoms_array)

	def bond_signature(agragates_array):
		#returns an array with all bonds in a molecule. This will be used to characterize the materials
		# molecules with different bond signature =  differente molecules

		N = len(agragates_array)
		molecular_types = np.zeros([N],dtype=object)
		
		
		for i in range(N):
			molecule = agragates_array[i]
			agg_mol  = []
			
			for atom in molecule:
				element  = atom.element
				neighbor = atom.neighbors
				element_neighbor = []
				
			
				for atom_aux in molecule:
					if (atom_aux.index in neighbor):		
						element_neighbor.append(atom_aux.element)
				neigh_test = [atom_test.element for atom_test in atoms_array]	
				
				element_neighbor = [element] + element_neighbor #novo 
				element_neighbor.sort() #must investigate if influences the output	
				#print(element_neighbor)
				agg_mol.append(element_neighbor)
		
			agg_mol.sort()		
			
			agg_mol = np.array(agg_mol,dtype=object)
			molecular_types[i] = agg_mol
			
		return molecular_types	

		
	def signature_string_to_int(molecular_types):
		# reads [ [bonds in mol1] , [bonds in mol 2] ...] and returns [1,2,3...] = [mat1,mat2,...]
		N = len(agragates_array)
		molecular_types_int = np.zeros([N],dtype=int)
		
		
		molecular_types_int[0] = 1
		
		for i in range(1,N):
			match = False
			
			for j in range(0,i):
				
				#if ( (molecular_types[i].all() == molecular_types[j].all())  and (i != j)) : #aqui
				#if ( np.all(molecular_types[i]==molecular_types[j])  and (i != j)) : #aqui			
				if ( (np.array_equal(molecular_types[i],molecular_types[j]))  and (i != j)) :
				
					#print(molecular_types[i],molecular_types[j])
					molecular_types_int[i] = molecular_types_int[j]
					match = True
					break
				
			if (match == False):
				molecular_types_int[i] = molecular_types_int[i-1] + 1
		return molecular_types_int	

	molecular_types     = bond_signature(agragates_array)
	molecular_types_int = signature_string_to_int(molecular_types)



	reduced_moltype_list = list(set(molecular_types_int.tolist()))

	#AQUI
	samples_arr = []
	for mol_type in reduced_moltype_list: #getting samples of each molecule type

		molecule_samples_indx = [i for i in range(len(molecular_types_int)) if molecular_types_int[i] == mol_type][0]
		mol_sample            = np.copy(agragates_array[molecule_samples_indx])
		samples_arr.append([mol_type-1,mol_sample]) #adjust Leo convention


	def get_stoichiometry(array_atoms):
		histogram = collections.Counter(array_atoms)
		soich = ''
		for element in histogram:
		
			ele_count = str(histogram.get(element))
			ele_name  = element
			soich = soich +"| "+ele_count+ele_name+" "
			
		soich = soich + "|"
		return soich




		
	#getting the elements of the neighbors of each atom within a molecule
	def get_bondelements(samples_arr,atoms_array):
		bond_samples = []
		
		for sample in samples_arr: #for every molecule type sample
			mol_type = sample[0]
			mol_arr  = sample[1]
			neigh_ele_arr     = []
			origin_atom       = []
			
			
			for atom in mol_arr: #for every atom inside this sample
				neigh     = atom.neighbors #can break if the atoms has no neighbors
				len_neigh = len(neigh)
				neigh_ele = np.zeros([len_neigh],dtype=object)
				origin_atom.append(atom.element)
				for i in range(len_neigh):
					indx = neigh[i]
					for atom2 in atoms_array:
						indx2 = atom2.index
						
						if(indx == indx2):
							neigh_ele[i] = atom2.element
							break
				neigh_ele = neigh_ele.tolist()
				neigh_ele_arr.append(neigh_ele)		
			
			bond_samples.append([mol_type,origin_atom,neigh_ele_arr])
			
		return 	bond_samples
		
	bond_samples = get_bondelements(samples_arr,atoms_array)	
	def print_bondsamples(bond_samples):

		matfile          = "matinfo.txt"
		format_mat       = "Material type: %i \n"
		format_info      = "Number of atoms = %i, Stoichiometry : %s \nBond structure: \n"
		format_origin_at = "%4s -->\t"
		
		with open(matfile,'w') as file_mat:
			for moltype in bond_samples:
			
				type_m          = moltype[0]
				origin_atom_arr = moltype[1]
				neigh_ele_arr   = moltype[2]
				n_at            = len(origin_atom_arr)
				
				soich = get_stoichiometry(origin_atom_arr)
				
				file_mat.write(format_mat  % (type_m)) 
				file_mat.write(format_info % (n_at,soich)) 
				
				for i in range(n_at):
					origin_atom_i   = origin_atom_arr[i]
					neigh_ele_arr_i = neigh_ele_arr[i]
					file_mat.write(format_origin_at % (origin_atom_i)+ '\t'.join([ k for k in neigh_ele_arr_i]) + '\n')		
				file_mat.write('\n')		
	print_bondsamples(bond_samples)


	N = len(agragates_array)
	center_masses = []
	for i in range(N):
		molecule      = agragates_array[i]
		molecule_type = molecular_types_int[i]
		tot_mass = sum([ at.mass for at in molecule ])
		
		x_cm = 0
		y_cm = 0
		z_cm = 0
		
		
		for atoms in molecule:
			pos  = atoms.position
			mass = atoms.mass
			
			x   = pos[0]
			y   = pos[1]
			z   = pos[2]
			
			x_cm = x_cm + x*mass
			y_cm = y_cm + y*mass
			z_cm = z_cm + z*mass
			
		output = [molecule_type,x_cm/tot_mass,y_cm/tot_mass,z_cm/tot_mass]	
		center_masses.append(output)
		
		
	center_masses_np = np.array(center_masses)
	X_cm = center_masses_np[:,1]
	Y_cm = center_masses_np[:,2]
	Z_cm = center_masses_np[:,3]


	X_lenght = np.amax(X_cm) - np.amin(X_cm)
	Y_lenght = np.amax(Y_cm) - np.amin(Y_cm)
	Z_lenght = np.amax(Z_cm) - np.amin(Z_cm)



	print("The unit cell with %i molecules has the dimensions of %f X %f X %f angstrons" %(len(center_masses_np),X_lenght,Y_lenght,Z_lenght))
	print()
	n_times = int(input("Insert the number of times that you want to duplicate the unit cell (integer!):"))

	#n_times = 2

	#generating multiple cells
	def multiply_unitcells(center_mass,n_times):
		
		
		new_lattice = np.copy(center_mass)
		n_cm = len(center_mass)
		
		for dx in range(n_times):
			for dy in range(n_times):
				for dz in range(n_times):
					
					new_cell = np.copy(center_mass)
					for i in range(n_cm):
						#multipling the X position in the unit cell dx times	
						#print(dx,dy,dz)			
						new_cell[i][1] = new_cell[i][1] +dx*X_lenght 
						new_cell[i][2] = new_cell[i][2] +dy*Y_lenght
						new_cell[i][3] = new_cell[i][3] +dz*Z_lenght
					
					new_lattice = np.vstack((new_lattice,new_cell))
		
		return new_lattice
		
	new_lattice = multiply_unitcells(center_masses_np,n_times)
	new_lattice = np.unique(new_lattice.round(decimals=5), axis=0) #removing duplicates

	X_cm_mult = new_lattice[:,1]
	Y_cm_mult = new_lattice[:,2]
	Z_cm_mult = new_lattice[:,3]


	X_lenght_mult = np.amax(X_cm_mult) - np.amin(X_cm_mult)
	Y_lenght_mult = np.amax(Y_cm_mult) - np.amin(Y_cm_mult)
	Z_lenght_mult = np.amax(Z_cm_mult) - np.amin(Z_cm_mult)

	print("Now, the lattice, with %i molecules, has the dimensions of %f X %f X %f angstrons" %(len(new_lattice),X_lenght_mult,Y_lenght_mult,Z_lenght_mult))

	#print(new_lattice)

	#if you want to see both unit and multiplied cells
	'''
	fig = plt.figure()
	ax = plt.axes(projection='3d')
	ax.scatter3D(X_cm_mult, Y_cm_mult, Z_cm_mult,c='b',marker='^');
	ax.scatter3D(X_cm, Y_cm, Z_cm, c='r',s=50);
	plt.show()	
	'''

		
	X_cm_mult = new_lattice[:,1]
	Y_cm_mult = new_lattice[:,2]
	Z_cm_mult = new_lattice[:,3]		
	Mats      = new_lattice[:,0] #careful, others funcs have mats in the last position		
	Mats = Mats -1 #adjusting to Leo's mat. convention
		
	return X_cm_mult,Y_cm_mult,Z_cm_mult,Mats
			
			
			
####################### END FUNCTIONS ###################

#####################################################################
############# FUNCS TO MESS AROUND WITH UNIT CELLs (option 3) #########

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
	
def func_parametrize(t):
	#return [t,t,0]
	return [np.cos(t),np.sin(t),t]
	#return [np.cos(t),np.sin(t),0]
	#return [(np.cos(t))**2,np.sin(t),-t]
	#return [((t+1)/(t-1)),((t-1)/(t+1)),0]
	#return  [np.exp(t)*np.sin(t),np.exp(-t)*np.cos(t),0]
	#return  [np.sin(t),np.sin(2*t),0]
	
#incomplete	
def parametrize_mat(par):
	cell_file_1 = par[0][0] 
		
	
	N_points = 500
	t_min    = 0
	t_max    = 20
	
	t_array  = np.linspace(t_min,t_max,N_points)

	func = [ func_parametrize(t) for t in t_array ] 

	func = np.array(func)	
	xfun = func[:,0] #points of the curve
	yfun = func[:,1]
	zfun = func[:,2]

	r0 = [xfun[0],yfun[0],zfun[0]]

	'''
	#plot of the function parametrized
	fig = plt.figure()
	ax = plt.axes(projection='3d')
	ax.scatter3D(xfun, yfun, zfun,c='r',marker='^');
	'''

	#getting the lowest values in each direction to produce a normalized lattice
	try:
		x_min = min( i for i in np.abs(xfun) if i != 0)
	except:
		x_min = 1
	try:	
		y_min = abs(min( j for j in np.abs(yfun) if j != 0))
	except:
		y_min = 1
	try:
		z_min = abs(min( k for k in np.abs(zfun) if k != 0))
	except:
		z_min = 1


	#normalizing curve
	xfun_norm = np.array([ int(x/x_min) for x in xfun ])
	yfun_norm = np.array([ int(y/y_min) for y in yfun ])
	zfun_norm = np.array([ int(z/z_min) for z in zfun ])



	trajectory  = []

	for i in range(N_points-1):

		x0 = xfun_norm[i]
		y0 = yfun_norm[i]
		z0 = zfun_norm[i]
		
		dx_next = xfun_norm[i+1] -x0
		dy_next = yfun_norm[i+1] -y0
		dz_next = zfun_norm[i+1] -z0
		
		#getting the displacement vector in each step
		displ_vec = np.array([dx_next,dy_next,dz_next])

				
		trajectory.append(displ_vec)

	#this function recieves the displacement vectors and returns the normalized lattice
	def vect_to_points(r0,trajectory):
		
		traj = []
		traj.append(r0)
		
		for i in range(1,len(trajectory)):
		
			x0 = traj[i-1][0]
			y0 = traj[i-1][1]
			z0 = traj[i-1][2]	
		
			x_new = x0 + trajectory[i][0]
			y_new = y0 + trajectory[i][1]
			z_new = z0 + trajectory[i][2]
			
			r   = [x_new,y_new,z_new]
			traj.append(r)
		return traj
		
	traj = vect_to_points(r0,trajectory)
	traj = np.array(traj)
		
	print("this draw will host %s unit cells" % (len(trajectory)))	
	x = traj[:,0]
	y = traj[:,1]
	z = traj[:,2]


	X_unit,Y_unit,Z_unit,Mats_unit = read_lattice(cell_file_1)
	dX = np.amax(X_unit) - np.amin(X_unit)
	dY = np.amax(Y_unit) - np.amin(Y_unit)
	dZ = np.amax(Z_unit) - np.amin(Z_unit)

	unit_cell = [[X_unit[i],Y_unit[i],Z_unit[i],Mats_unit[i]] for i in range(len(X_unit)) ]
	unit_cell = np.array(unit_cell)

	mult_cell = unit_cell.copy()
	new_cell = unit_cell.copy()


	draw = traj
	for unit_vec in draw:

		new_cell = unit_cell.copy()
		
		for site in new_cell:
			
			site[0] = site[0] + dX*unit_vec[0]
			site[1] = site[1] + dY*unit_vec[1]
			site[2] = site[2] + dZ*unit_vec[2]
			
		mult_cell = np.vstack((mult_cell,new_cell))
		
		
	X    = mult_cell[:,0]	
	Y    = mult_cell[:,1]
	Z    = mult_cell[:,2]
	Mats = mult_cell[:,3]	
	
	#fig = plt.figure()
	#ax = plt.axes(projection='3d')
	#ax.scatter3D(X, Y, Z,c='b',marker='^');
	#plt.show()	
	
	return X,Y,Z,Mats
	
	
def heterojunction(par):
	cell_left  = par[0][0]
	cell_right = par[1][0]
	n_times    = int(float(par[2][0]))
	axis       = par[3][0]
	
	#direction in which we will join the two sides
	if axis == 'X':
		vec_joint = [1,0,0]
	if axis == 'Y':
		vec_joint = [0,1,0]
	if axis == 'Z':
		vec_joint = [0,0,1]
	
	
	X_left ,Y_left ,Z_left ,Mats_left  = read_lattice(cell_left)
	X_right,Y_right,Z_right,Mats_right = read_lattice(cell_right)
	
	
	non_duplicate_mat = set(Mats_right.tolist())
	shift_mat = input("The cell on the right has %s types of molecules. Do you wish to relabel its indexes?(y/n)" %(len(non_duplicate_mat)))
	
	
	
	#case if you want to change the labeling of right mats
	if(shift_mat == "y"):
		print("Currenly, the right input has the following mat indexes:")
		#old_index = '\t'.join([str(int(mat)) for mat in non_duplicate_mat ])
		old_index = ([str(int(mat)) for mat in non_duplicate_mat ])		
		print(old_index)
		print("Now, write the list that will substitute the indexing:")
		
		new_index = [str(int(item)) for item in input("Enter the index (separated by space): ").split()]
		print(new_index)
		#print(Mats_right)
		new_Mats_right = Mats_right.copy()
		
		n_mats = len(Mats_right)

		for i in range(n_mats):
	
			mats = Mats_right[i]
			indx = old_index.index(str(int(mats)))	
			new_Mats_right[i] = new_index[indx]

	
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
	
	n_times_x_right = int(X_tot_left/dX_right)
	n_times_y_right = int(Y_tot_left/dY_right)
	n_times_z_right = int(Z_tot_left/dZ_right)
	
	n_times_right = [n_times_x_right,n_times_y_right,n_times_z_right]
	multiplied_right = multiply_lattice(unit_cell_right,n_times_right,[dX_right,dY_right,dZ_right])
	

	shift_vec = [X_tot_left,Y_tot_left,Z_tot_left]
	shift_vec = [ shift_vec[i]*vec_joint[i] for i in range(len(shift_vec))]	 
	
	#shifting the entire right lattice 
	for site in multiplied_right:
		site[0] = site[0] + shift_vec[0]
		site[1] = site[1] + shift_vec[1] 
		site[2] = site[2] + shift_vec[2]
	
	#modeling the intersection
	#returns the number of sites along a determined axis
	def filter_lattice(lattice,axis):
		x = lattice[:,0]
		y = lattice[:,1]
		z = lattice[:,2]
		
		if axis == "X":
			y_min = np.amin(y)
			fixed_latt = lattice[y==y_min]
			z = fixed_latt[:,2]
			z_min = np.amin(z)
			fixed_latt = fixed_latt[z==z_min]
			
			return len(fixed_latt)
			
		
		if axis == "Y":#y
			x_min = np.amin(x)
			fixed_latt = lattice[x==x_min]
			z = fixed_latt[:,2]
			z_min = np.amin(z)
			fixed_latt = fixed_latt[z==z_min]
			
			return len(fixed_latt)			

		if axis == "Z":#z
			x_min = np.amin(x)
			fixed_latt = lattice[x==x_min]
			y = fixed_latt[:,1]
			y_min = np.amin(y)
			fixed_latt = fixed_latt[y==y_min]
			
			return len(fixed_latt)
	#how many sites exist along an axis?		
	n_sites_axis_left  = filter_lattice(multiplied_left,axis)
	n_sites_axis_right = filter_lattice(multiplied_right,axis)	
	
	dens_left  = dX_left/n_sites_axis_left #density of sites along an axis on the left side
	dens_right = dX_right/n_sites_axis_right #density of sites along an axis on the right side
	
	#picking the lowest dx among the two sides
	if dens_left < dens_right:
		dens = dens_left
	else:
		dens = dens_right
	
	def dist(vec1,vec2):
		return sqrt((vec1[0]-vec2[0])**2 + (vec1[1]-vec2[1])**2 +(vec1[2]-vec2[2])**2)
		
	#identifies the sites located near the heterojunction	
	def is_in_the_heterojunc(left_lattice,right_lattice,cutoff):
		right_het, left_het = [], []
		n_right = len(right_lattice)
		n_left  = len(left_lattice)
		
		for i in range(n_left):
			site_left = left_lattice[i][0:3]
			for j in range(n_right):
				site_right = right_lattice[j][0:3]

				distance = dist(site_left,site_right)
				
				#if the distance between a site on the left and right is comparable to the first neighbor
				# distance in the left side (a way to locate the sites in the heterojunc)
				if distance <= cutoff:
					#print(i,j,distance)
					right_het.append(j)
					left_het.append(i)
					
		return 	left_het, right_het		
		
	#add a disturbance in the sites located in the heterojunction	
	def shift_hetero(latt,vec,indx_ar,mean,sigma,scale):
		latt_ref = latt.copy()
		for indx in indx_ar:
			shift = np.array([ scale*pos*np.random.normal(mean,sigma) for pos in vec])
			latt[indx][0:3] = latt_ref[indx][0:3] + shift
			#print(shift)
	
	
	
	
	mean = 1
	sigma = mean
	scale = dens/1


	#getting the sites in the heterojunc
	left_het, right_het = is_in_the_heterojunc(multiplied_left,multiplied_right,dens)

	#shifting these sites
	shift_hetero(multiplied_left,vec_joint,left_het,mean,sigma,scale)
	shift_hetero(multiplied_right,vec_joint,right_het,mean,sigma,scale)

	#getting the sites that are too close
	left_het, right_het = is_in_the_heterojunc(multiplied_left,multiplied_right,dens)
	shift_hetero(multiplied_left,vec_joint,left_het,dens,0.2,1)
	shift_hetero(multiplied_right,vec_joint,right_het,dens,0.2,-1)
	# END INTERSECTION
	
		
	hetero_lattice = np.vstack((multiplied_left,multiplied_right))
	X    = hetero_lattice[:,0]
	Y    = hetero_lattice[:,1]
	Z    = hetero_lattice[:,2]
	Mats = hetero_lattice[:,3]
	
	return 	X,Y,Z,Mats
###########################################################
def write_lattice(lattice_name,func_name,output_par,X,Y,Z,latt_length):
    with open(lattice_name,'w') as f:
        f.write( ('#Func name: %s \n') %(func_name))
        f.write('\n'.join('#%s %s' % x for x in output_par))
        f.write('\n')
        for l in range(latt_length):
            line = [X[l],Y[l],Z[l],Mats[l]]
            f.write('\t'.join(["{:<10f} ".format(i) for i in line]) + '\n')
            
def draw_lattice(X,Y,Z,Mats,color_dir,fig_name):
    fig = plt.figure()
    ax = plt.axes(projection='3d')
    ax.scatter3D(X, Y, Z,c=color_dir,marker='^');
    
    try:
        plt.show()    
    except:
        plt.switch_backend('agg')
        plt.savefig(fig_name+'.png')
'''
colors_dic = {0:'black', 1:'blue', 2:'red', 3:'green', 4:'yellow'}
X,Y,Z,Mats = heterojunction([['lattice.txt'],['lattice.txt'],['2'],['X']])
colors = np.array([colors_dic.get(int(mat)) for mat in Mats])
fig = plt.figure()
ax = plt.axes(projection='3d')
ax.scatter3D(X, Y, Z,c=colors,marker='^');
plt.show()  
'''

'''
colors_dic = {0:'black', 1:'blue', 2:'red', 3:'green', 4:'yellow'}
X,Y,Z,Mats = lattice_BHJ([ [5000],[1,1,0], [0.5,0.5],[10],[1] ])
colors = np.array([colors_dic.get(int(mat)) for mat in Mats])
fig = plt.figure()
ax = plt.axes(projection='3d')
ax.scatter3D(X, Y, Z,c=colors,marker='^');
plt.show()  
exit()
'''

class morphology_function:
    def __init__(self,name,func,param_list):
    	self.name = name
    	self.func = func
    	self.param_list = param_list
    def __call__(self,param):
    	return self.func(param)	
    

#funcs option 2
lattice_dir = ["Number of Sites","displacement_vec","distribuition_vec"] #should have the same order of the function
BHJ_dir     = ["Number of Sites","displacement_vec","distribuition_vec","Number of times to loop","Cutoff distance (angstrom)"] #should have the same order of the function
loadcif_dir = ["mol_filename"]

latt_func      = morphology_function("lattice",lattice,lattice_dir)
BHJ_func       = morphology_function("bulk heterojunction (BHJ)",lattice_BHJ,BHJ_dir)
loadcif_func   = morphology_function("loadcif",load_cif,loadcif_dir)
func_list      = [latt_func,loadcif_func,BHJ_func] #list of functions to be included in option 2

#funcs option 3
parametrize_dir    = ["filename for the first unit cell"]
heterojunction_dir = ["name_left_cell","name_right_cell","number_of_reps","Joining axis (X,Y or Z)"]

paramet_func           = morphology_function("parametrized",parametrize_mat,parametrize_dir)
heterojunction_func    = morphology_function("heterojunction",heterojunction,heterojunction_dir)

func_opt3_list  = [paramet_func,heterojunction_func] #list of functions to be included in option 3

##########################################################

colors_dic = {0:'black', 1:'blue', 2:'red', 3:'green', 4:'yellow'}

print("AUX PROGRAM TO GENERATE LATTICE FOR KMC SIMULATIONS")
print("version: 1.0v, date 4/4/21 Author: Leo and Tiago")
print()
print("Select one of the options:")
print("1) I want to see some lattice.")
print("2) I want to generate a morhpology using one of our functions.")
print("3) I already have a lattice, just want to mess around with it.")

lattice_output = "lattice.txt"

choice1 = input()
if ( choice1 == "1"):
    morph_file = input("Ok! Give to me the file name of the data: ")
    X,Y,Z,Mats = read_lattice(morph_file)
    colors = np.array([colors_dic.get(int(mat)) for mat in Mats])
    draw_lattice(X,Y,Z,Mats,colors,'fig_option1.png')
    
#if you want to create a lattice     
if ( choice1 == "2"):
    list_name_func = [ func.name for func in func_list ]
    n = len(list_name_func) 
    list_name_func = '\t'.join([str(i)+"  "+str(list_name_func[i]) for i in range(n) ]) #list with the names of the functions
    print("This version has %s pre-written morphology functions. Choose one: " %(len(func_list)))
    print(list_name_func) #printing the available functions	
    choice2 = int(input())	
    func = func_list[choice2] #choosing a function	

    par_list_name = func.param_list
    
    list_name_pars = '\t'.join([str(par) for par in par_list_name ])    #list of parameters for a given func		
    print("This function has %s parameters. Namely," %(len(par_list_name)))
    print(list_name_pars)
    par = []
    
    for par_set in par_list_name: #looping through each parameter
        print()
        print(par_set) 	
        par_loc = [item for item in input("Enter the parameter's entries (separated by space): ").split()] #needed if a parameter is an array
        print(par_loc)
        par.append(par_loc)

    output_parameters = [ (par_list_name[i],par[i]) for i in range(len(par_list_name))]
       	
    
    X,Y,Z,Mats = func(par)
    colors     = np.array([colors_dic.get(int(mat)) for mat in Mats])
    draw_lattice(X,Y,Z,Mats,colors,'lattice')
    
    #writing the lattice.txt file
    write_lattice(lattice_output,func.name,output_parameters,X,Y,Z,len(X))

            


#if you want to mess with the lattice
if ( choice1 == "3"):
    list_name_func = [ func.name for func in func_opt3_list ]
    n = len(list_name_func) 
    list_name_func = '\t'.join([str(i)+"  "+str(list_name_func[i]) for i in range(n) ]) #list with the names of the functions
    print("This version has %s pre-written morphology functions. Choose one:" %(len(func_opt3_list)))
    print(list_name_func) #printing the available functions	
    choice2 = int(input())	
    func = func_opt3_list[choice2] #choosing a function	

    par_list_name = func.param_list
    
    list_name_pars = '\t'.join([str(par) for par in par_list_name ])    #list of parameters for a given func		
    print("This function has %s parameters. Namely," %(len(par_list_name)))
    print(list_name_pars)
    par = []
    
    for par_set in par_list_name: #looping through each parameter
        print()
        print(par_set) 	
        par_loc = [item for item in input("Enter the parameter's entries (separated by space): ").split()] #needed if a parameter is an array
        print(par_loc)
        par.append(par_loc)
        
    output_parameters = [ (par_list_name[i],par[i]) for i in range(len(par_list_name))]        
 
           
    X,Y,Z,Mats = func(par)
    colors = np.array([colors_dic.get(int(mat)) for mat in Mats])

    write_lattice(lattice_output,func.name,output_parameters,X,Y,Z,len(X))
    draw_lattice(X,Y,Z,Mats,colors,'lattice')
 




