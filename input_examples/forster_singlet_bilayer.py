import kmc.morphology as morphology
from kmc.rates import *
from kmc.particles import *

###BASIC PARAMETERS######################################################################
identifier         = 'forster_singlet' #output identifier
time_limit         = np.inf# in PS
animation_mode     = True
save_animation     = False # if you want to save the animation
animation_exten    = 'gif' # possible options ('gif' and 'mp4')
marker_type        = 1     # marker type used at the animation processs ( 0 = balls, 1 = symbols) 
pause              = False  # if you want that the annimation stops in the first frame (debug purposes)
rounds             = 100   # Number of rounds
n_proc             = 1     # Number of cores to be used
#########################################################################################

###SINGLET EXCITONS######################################################################

##FORSTER RADII (Å)
r00   = 25   #Forster radius material 0 --> material 0 (Angstrom)    
r01   = 25   #material 0 --> material 1      
r10   = 25       
r11   = 25     
raios = {(0,0):r00, (0,1):r01, (1,0):r10, (1,1):r11}

##FLUORESCENCE LIFETIMES (PS)
f0 = 29000 #lifetime of material 0
f1 = 2900 #lifetime of material 1
lifetimes = {0:f0,1:f1}

##TANSITION DIPOLE MOMENTS (a.u.)
mu0 = 2.136
mu1 = 5.543
mus = {0:mu0,1:mu1}

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
#lattice_func_par   = ["lattice_heterojun.example"] # file name of the system's morphology


# Creating a new lattice at each new round
lattice_func      = morphology.bilayer
n_times           = 2   #number of duplications
axis              = 'X' #axis of the junction
displacement_vect = [ 5, 5, 0]
num_sites         = 200
distribu_vect     = [1]
latt_param        = [num_sites,displacement_vect,distribu_vect]
lattice_func_par  = [n_times,axis,latt_param]


##ENERGIES
#Gaussian distribuitions
t1s   = {0:(3.7,0.0), 1:(3.7,0.0), 'level':'t1'} #(Peak emission energy (eV), disperison (eV)
s1s   = {0:(6.1,0.0), 1:(6.1,0.0), 'level':'s1'} # triplet energy, disperison (eV)

ener_function      = [morphology.Gaussian_energy(s1s),morphology.Gaussian_energy(t1s)]  
#########################################################################################


##GENERATE PARTICLES#####################################################################
num_ex             = 20     #number of particles

#Type of particle
gen_function        = morphology.gen_excitons

#Choose the way that the particles will be distribuited
sel_func    = morphology.filter_selection
sel_params  = {'shape_dic': morphology.shape_dic, 'mat' : [None],
 'shape': "rectangle", 'origin': None, 'argum' : [[130,150],[0,120],[0,1]]}
#########################################################################################

##BIMOLECULAR OPTIONS###################################################################
bimolec               = True  # Turn on annihilation
##list of all annihi funcs that will be used
bimolec_funcs_array = [morphology.ele_hol_recomb,morphology.anni_sing] 
#########################################################################################

   

