import kmc.morphology as morphology
from kmc.rates import *
from kmc.particles import *

###BASIC PARAMETERS######################################################################
identifier         = 'charge_creation' #output identifier
time_limit         = np.inf# in PS
animation_mode     = True
save_animation     = False # if you want to save the animation
animation_exten    = 'gif' # possible options ('gif' and 'mp4')
marker_type        = 1     # marker type used at the animation processs ( 0 = balls, 1 = symbols) 
pause              = False # if you want that the annimation stops in the first frame (debug purposes)
rounds             = 100   # Number of rounds
n_proc             = 1     # Number of cores to be used
#########################################################################################

###SINGLET EXCITONS######################################################################

##FORSTER RADII (Å)
r00   = 25   #Forster radius material 0 --> material 0 (Angstrom)    
r01   = 25   #material 0 --> material 1      
r10   = 25       
r11   = 25     
radii = {(0,0):r00, (0,1):r01, (1,0):r10, (1,1):r11}

##FLUORESCENCE LIFETIMES (PS)
f0 = 2900000 #lifetime of material 0
f1 = 2900 #lifetime of material 1
lifetimes = {0:f0,1:f1}

##TANSITION DIPOLE MOMENTS (a.u.)
mu0 = 2.136
mu1 = 5.543
mus       = {0:mu0,1:mu1}

##EXCITION TRANSFER RATES
forster   = Forster(Rf=radii,life=lifetimes,mu=mus)

##FLUORESCENCE RATES
fluor     = Fluor(life=lifetimes)


##DISSOCIATION RATE
H = {(0,0):10E12,(0,1):10E12,(1,0):10E12,(1,1):10E12}
invrad = {0:0.1,1:0.1}
dissociation = Dissociation(AtH=H,invrad=invrad,T=300)

###TRIPLET EXCITONS######################################################################

##DEXTER RADII (Å)
Rds = {(0,0):10, (0,1):0, (1,0):0, (1,1):10}

##PHOSPHORESCENCE LIFETIMES (PS)
phlife = {0:5.29,1:5.29}

##SUM OF VAN DER WALLS RADII (Å)
Ls = {0:5.0,1:5.0}

##TRIPLET TRANSFER RATES
dexter = Dexter(Rd=Rds,life=phlife,L=Ls)

##PHOSPHORESCENCE RATE
phosph = Phosph(life=phlife)
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

processes = {'singlet':[forster], 'triplet':[dexter], 'electron':[miller],'hole':[miller]}
monomolecular = {'singlet':[fluor],'triplet':[phosph],'electron':[],'hole':[]}
#########################################################################################

###MORPHOLOGY############################################################################

##Morphology functions

#Reading a file name that contains your lattice
#lattice_func = morphology.read_lattice
#lattice_func_par   = ["lattice.example"] # file name of the system's morphology


# Creating a new lattice at each new round
num_sites         = 100             #number of sites of the lattice
displacement      = [5, 5, 0]       #vector of the unit cell
disorder          = [0.5,0.5,0.5]   #std deviation from avg position
composition       = [0.5,0.5]       #popuation probility Ex.: distribu_vect[0] is the prob of mat 0 appear in the lattice
lattice_func      = morphology.Lattice(num_sites,displacement,disorder,composition)



##ENERGIES
#Gaussian distribuitions
lumos = {0:(-3.7,0.0), 1:(-3.7,0.0), 'level':'lumo'} #(Peak emission energy (eV), disperison (eV)
homos = {0:(-6.1,0.0), 1:(-6.1,0.0), 'level':'homo'} # triplet energy, disperison (eV)
t1s   = {0:(3.7,0.0), 1:(3.7,0.0), 'level':'t1'} #(Peak emission energy (eV), disperison (eV)
s1s   = {0:(6.1,0.0), 1:(6.1,0.0), 'level':'s1'} # triplet energy, disperison (eV)

ener_function      = [morphology.Gaussian_energy(s1s),morphology.Gaussian_energy(t1s),morphology.Gaussian_energy(homos),morphology.Gaussian_energy(lumos)]  
#########################################################################################

##GENERATE PARTICLES#####################################################################
num_ex             = 1    #number of particles

#Type of particle
gen_function       = morphology.gen_pair_elechole
#gen_function       = morphology.gen_electron

#Choose the way that the particles will be distribuited
sel_func    = morphology.filter_selection
sel_params  = {'shape_dic': morphology.shape_dic, 'mat' : [None],
 'shape': "free", 'origin': None, 'argum' : None}
 
#########################################################################################

##BIMOLECULAR OPTIONS###################################################################
bimolec               = True  # Turn on annihilation
##list of all annihi funcs that will be used
bimolec_funcs_array = [morphology.ele_hol_recomb]#,morphology.anni_sing] 
#########################################################################################

   

