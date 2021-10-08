import kmc.morphology as morphology
from kmc.kmc_classes import *

###BASIC PARAMETERS######################################################################
identifier         = '3mats_singlet' #output identifier
time_limit         = np.inf# in PS
animation_mode     = True
save_animation     = False # if you want to save the animation
animation_exten    = 'gif' # possible options ('gif' and 'mp4')
marker_type        = 1     # marker type used at the animation processs ( 0 = balls, 1 = symbols) 
pause              = False # if you want that the annimation stops in the first frame (debug purposes)
rounds             = 100   # Number of rounds
n_proc             = 6     # Number of cores to be used
#########################################################################################

###SINGLET EXCITONS######################################################################

##FORSTER RADII (Å)
r00   = 25   #Forster radius material 0 --> material 0 (Angstrom)    
r01   = 25   #material 0 --> material 1      
r10   = 25       
r11   = 25
r22   = 25
r21   = 25
r20   = 25
r12   = 25
r02   = 25
     
raios = {(0,0):r00, (0,1):r01, (1,0):r10, (1,1):r11, (0,2):r02, (2,0):r10, (2,1):r21, (1,2):r12, (2,2):r22}

##FLUORESCENCE LIFETIMES (PS)
f0 = 2900000 #lifetime of material 0
f1 = 2900    #lifetime of material 1
f2 = 2900    #lifetime of material 2
lifetimes = {0:f0,1:f1,2:f2}

##TANSITION DIPOLE MOMENTS (a.u.)
mu0 = 2.136
mu1 = 5.543
mu2 = 3
mus       = {0:mu0,1:mu1,2:mu2}

##EXCITION TRANSFER RATES
forster   = Forster(Rf=raios,life=lifetimes,mu=mus)

##FLUORESCENCE RATES
fluor     = Fluor(life=lifetimes)


##relative permitivity
relative_eps       = 3.5   




###PROCESSES#############################################################################

processes = {'singlet':[forster], 'triplet':[], 'electron':[],'hole':[]}
monomolecular = {'singlet':[fluor],'triplet':[],'electron':[],'hole':[]}
#########################################################################################

###MORPHOLOGY############################################################################

##Morphology functions

#Reading a file name that contains your lattice
#lattice_func = morphology.read_lattice
#lattice_func_par   = ["lattice_3mat.example"] # file name of the system's morphology


# Creating a new lattice at each new round
lattice_func      = morphology.lattice
displacement_vect = [ 5, 5, 5]
num_sites         = 100
distribu_vect     = [0.5,0.2,0.3]
lattice_func_par  = [num_sites,displacement_vect,distribu_vect]


##ENERGIES
#Gaussian distribuitions
t1s   = {0:(3.7,0.0), 1:(3.7,0.0), 2:(3.7,0.0), 'level':'t1'} #(Peak emission energy (eV), disperison (eV)
s1s   = {0:(6.1,0.0), 1:(6.1,0.0), 2:(6.1,0.0), 'level':'s1'} # triplet energy, disperison (eV)

ener_function      = [morphology.Gaussian_energy(s1s),morphology.Gaussian_energy(t1s)]  
#########################################################################################


##GENERATE PARTICLES#####################################################################
num_ex             = 20     #number of particles

#Type of particle
gen_function        = morphology.gen_excitons

#Choose the way that the particles will be distribuited
sel_func    = morphology.filter_selection
sel_params  = {'shape_dic': morphology.shape_dic, 'mat' : [None],
 'shape': "free", 'origin': None, 'argum' : None}
#########################################################################################

##ANNIHILATION OPTIONS###################################################################
anni               = True  # Turn on annihilation
##list of all annihi funcs that will be used
annihi_funcs_array = [morphology.anni_ele_hol,morphology.anni_sing] 
#########################################################################################

   

