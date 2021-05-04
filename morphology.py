import numpy as np
import random
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d
import collections
################################
# Code block of KMC that treats the morphology
# The file organizes as follows:
# morphology functions
# creation of its respective classes
# user-machine interaction to choose the function
# writing on lattice.txt the morphology


#MORPHOLOGY FUNCS 
#num_molecs: number of molecules
#vector: 3 component lattice vector. For lower dimension, make distance 0 
#ps: vector with the relative probability of each material
def lattice(param):
	
    num_molecs =  int(param[0][0]) #adjusting to users input     
    vector     =  [ float(x) for x in param[1] ]
    ps         =  [ int(x)   for x in param[2] ]
    
    
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


def load_cif(mol_file):

	#This program reads .cif files through a .cube.cc1 file

	print("----------------------------------------------------------------------------------------------------")
	print("༼つಠ益ಠ༽つ                          LOAD .CIF PROGRAM (via .mol extension)                       ༼つಠ益ಠ༽つ")
	print("author: Tiago Cassiano, version: 1.1v (hoping that will sufice)")
	print("STEPS TO USE THE PROGRAM:")
	print("1) Pick a .cif geometry somewhere with your desired unit cell")
	print("2) Install Mercury software (https://www.ccdc.cam.ac.uk/support-and-resources/Downloads/)")
	print("3) Load .cif there")
	print("		3.1) On the display box (left-side bellow), check packing")
	print("		3.2) Check if the geometry seems like it should be! (Special care if there is hydrogen bonds)")
	print("		3.3) File --> Save as")
	print("		3.4) Choose Mol2 files (*.mol2 *.mol) ")
	print("4) Insert the file name for this program")
	print("----------------------------------------------------------------------------------------------------")
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

	print(molecular_types_int)

	def print_onemat(indx,bond_sign,file_var):

		format_mat = "Material type: %i \n"
		file_var.write(format_mat % (indx)) 
		file_var.write("Bond structure: \n")
		for line in bond_sign:
			file_var.write('\t'.join([ i for i in line]) + '\n')		
		file_var.write('\n')
		
		
	def print_matinfo(molecular_types,molecular_types_int):
		mat_unique_indexes = []
		mat_unique         = []
		molecular_types_int_list = molecular_types_int.tolist()
		
		
		mat_unique_indexes.append(0)
		mat_unique.append(molecular_types_int[0])
		
		for ele in molecular_types_int_list:
			if ele not in mat_unique:
				indx = molecular_types_int_list.index(ele)
				mat_unique_indexes.append(indx)
				mat_unique.append(ele)

		
		mat_len = len(mat_unique_indexes)
		with open(matinfo_file, 'w') as file_var:
			for i in range(mat_len):
				mat_indx  = mat_unique[i]
				bond_sign = molecular_types[mat_unique_indexes[i]]
				print_onemat(mat_indx,bond_sign,file_var)
		
		
		
	print_matinfo(molecular_types,molecular_types_int)


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
	fig = plt.figure()
	ax = plt.axes(projection='3d')
	ax.scatter3D(X_cm_mult, Y_cm_mult, Z_cm_mult,c='b',marker='^');
	ax.scatter3D(X_cm, Y_cm, Z_cm, c='r',s=50);
	plt.show()	


		
	X_cm_mult = new_lattice[:,1]
	Y_cm_mult = new_lattice[:,2]
	Z_cm_mult = new_lattice[:,3]		
	Mats      = new_lattice[:,0]		
	Mats = Mats -1 #adjusting to Leo's mat. convention
		
	return X_cm_mult,Y_cm_mult,Z_cm_mult,Mats
			
			
			
####################### END FUNCTIONS ###################





    
class morphology_function:
    def __init__(self,name,func,param_list):
    	self.name = name
    	self.func = func
    	self.param_list = param_list
    def __call__(self,param):
    	return self.func(param)	
    


lattice_dir = ["dim","displacement_vec","distribuition_vec"] #should have the same order of the function
loadcif_dir = ["mol_filename"]

latt_func      = morphology_function("lattice",lattice,lattice_dir)
loadcif_func   = morphology_function("loadcif",load_cif,loadcif_dir)
func_list      = [latt_func,loadcif_func]

#####################################################################


print("Kinectic Monte Carlo Simulations for Excitonic Dynamics LeoKMC")
print("version: 1.0v, date 4/4/21 Author: Leo")
print()
print("Select one of the options:")
print("1) I want to run a simulation using a morphology that was generated already")
print("2) I want to generate a morhpology using one of our functions")

lattice_output = "lattice.txt"

choice1 = input()
if ( choice1 == "1"): #not implemented yet
    morph_file = input("Ok! Give to me the file name of the data:")
    exit()
     
if ( choice1 == "2"):
    list_name_func = [ func.name for func in func_list ]
    n = len(list_name_func) 
    list_name_func = '\t'.join([str(i)+"  "+str(list_name_func[i]) for i in range(n) ]) #list with the names of the functions
    print("This version has %s pre-written morphology functions. Choose one:" %(len(func_list)))
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

X,Y,Z,Mats = func(par)
N = len(X)
#writing the lattice.txt file
with open(lattice_output,'w') as f:
	for l in range(N):
		line = [X[l],Y[l],Z[l],Mats[l]]
		f.write('\t'.join(["{:<10f} ".format(i) for i in line]) + '\n')	


