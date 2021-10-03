import kmc.morphology as morphology
from kmc.kmc_classes import *

###BASIC PARAMETERS######################################################################
identifier         = 'charge_creation' #output identifier
time_limit         = np.inf
animation_mode     = True
save_animation     = False  # if you want to save the animation
animation_exten    = 'gif' # possible options ('gif' and 'mp4')  
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
raios = {(0,0):r00, (0,1):r01, (1,0):r10, (1,1):r11}

##FLUORESCENCE LIFETIMES (PS)
f0 = 2900000 #lifetime of material 0
f1 = 2900 #lifetime of material 1
lifetimes = {0:f0,1:f1}

##TANSITION DIPOLE MOMENTS (a.u.)
mu0 = 2.136
mu1 = 5.543
mus       = {0:mu0,1:mu1}

##EXCITION TRANSFER RATES
forster   = Forster(Rf=raios,life=lifetimes,mu=mus)

##FLUORESCENCE RATES
fluor     = Fluor(life=lifetimes)


##NONRADIATIVE DECAY RATES
nonrad  = {0:0,1:0}
nonradiative = Nonrad(rate=nonrad)

##DISSOCIATION RATE
H = {(0,0):10E12,(0,1):10E12,(1,0):10E12,(1,1):10E12}
invrad = {0:0.1,1:0.1}
dissociation = Dissociation(AtH=H,invrad=invrad,T=300)


###CHARGES###############################################################################

##relative permitivity
relative_eps       = 3.5   

##ATTEMPT-TO-HOP FREQUENCY (?)
H = {(0,0):10E12,(0,1):10E12,(1,0):10E12,(1,1):10E12}

##INVERSE LOCALIZATION RADIUS ()
invrad = {0:1.5,1:1.5}

##MILLER-ABRAHAMS RATE
miller = MillerAbrahams(AtH=H,invrad=invrad,T=300)


###PROCESSES#############################################################################

processes = {'singlet':[forster], 'triplet':[], 'electron':[miller],'hole':[miller]}
monomolecular = {'singlet':[fluor],'triplet':[],'electron':[],'hole':[]}
#########################################################################################

###MORPHOLOGY############################################################################

##Morphology functions
lattice_filename   = "lattice.example" # file name of the system's morphology
X,Y,Z,Mats = morphology.read_lattice(lattice_filename)

##ENERGIES
#Gaussian distribuitions
s1s = {0:(3.7,0.0), 1:(2.85,0.0)} #(Peak emission energy (eV), disperison (eV)
t1s = {0:(6.1,0.0), 1:(5.25,0.0)} # triplet energy, disperison (eV)


#ener_function     = morphology.s1_t1_distr
ener_function      = morphology.homo_lumo
parameters_enefunc = [s1s, t1s, Mats]
s1, t1 = ener_function(parameters_enefunc)    
ene_dic = {'s1':s1, 't1':t1, 'HOMO':t1,'LUMO':s1} #careful, if you choose dissociation, you also must give HOMO and LUMO
#########################################################################################

##GENERATE PARTICLES#####################################################################
num_ex             = 10     #number of particles

#Type of particle
gen_function       = morphology.gen_pair_elechole


#Shape of the particle's generation
#Getting filter funcs from morphology		
selection = morphology.filter_selection(X,Y,Z,Mats,morphology.shape_dic,mat=[None],shape="free",origin=None,argum=None) 

parameters_genfunc = [num_ex,selection]
#########################################################################################

##ANNIHILATION OPTIONS###################################################################
anni               = True  # Turn on annihilation
##list of all annihi funcs that will be used
annihi_funcs_array = [morphology.anni_ele_hol,morphology.anni_sing] 
#########################################################################################

   

