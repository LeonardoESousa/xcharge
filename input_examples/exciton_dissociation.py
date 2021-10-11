import kmc.morphology as morphology
from kmc.rates import *
from kmc.particles import *

###BASIC PARAMETERS######################################################################
identifier         = 'charge_creation' # output identifier
time_limit         = np.inf            # in ps
animation_mode     = True
save_animation     = False             # if you want to save the animation
animation_exten    = 'gif'             # possible options ('gif' and 'mp4')
marker_type        = 1                 # marker type used at the animation processs ( 0 = balls, 1 = symbols) 
pause              = False             # if you want that the annimation stops in the first frame (debug purposes)
rounds             = 100               # Number of rounds
n_proc             = 1                 # Number of cores to be used
#########################################################################################

###SINGLET EXCITONS######################################################################

##FORSTER RADII (Å)
r00   = 25   #Forster radius material 0 --> material 0 (Angstrom)    
r01   = 25   #material 0 --> material 1      
r10   = 25       
r11   = 25     
radii = {(0,0):r00, (0,1):r01, (1,0):r10, (1,1):r11}

##FLUORESCENCE LIFETIMES (ps)
f0 = 2900000 #lifetime of material 0
f1 = 2900    #lifetime of material 1
lifetimes = {0:f0,1:f1}

##TANSITION DIPOLE MOMENTS (a.u.)
mu0 = 2.136
mu1 = 5.543
mus = {0:mu0,1:mu1}

##EXCITON TRANSFER RATES
forster   = Forster(Rf=radii,life=lifetimes,mu=mus)

##FLUORESCENCE RATES
fluor     = Fluor(life=lifetimes)


##DISSOCIATION RATE
H = {(0,0):1e9,(0,1):1e9,(1,0):1e9,(1,1):1e9}
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

##ATTEMPT-TO-HOP FREQUENCY (s^-1)
H = {(0,0):1e12,(0,1):1e12,(1,0):1e12,(1,1):1e12}

##INVERSE LOCALIZATION RADIUS (Å^-1)
invrad = {0:1.5,1:1.5}

##MILLER-ABRAHAMS RATE
miller = MillerAbrahams(AtH=H,invrad=invrad,T=300)
 

###PROCESSES#############################################################################

processes = {'singlet':[forster,dissociation], 'triplet':[dexter], 'electron':[miller],'hole':[miller]}
monomolecular = {'singlet':[fluor],'triplet':[phosph],'electron':[],'hole':[]}
#########################################################################################

###MORPHOLOGY############################################################################

##Morphology functions

#Reading a file name that contains your lattice
#file = 'lattice.example'
#lattice_func = morphology.ReadLattice(file)

# Creating a new lattice at each new round
num_sites         = 400             #number of sites of the lattice
displacement      = [5, 5, 0]       #vector of the unit cell
disorder          = [0, 0, 0]       #std deviation from avg position
composition       = [0.5,0.5]       #popuation probility Ex.: distribu_vect[0] is the prob of mat 0 appear in the lattice
lattice_func      = morphology.Lattice(num_sites,displacement,disorder,composition)


#Electric Field and Dielectric constant
Electric_class    = morphology.Electric(eps=1.5,field=np.array([0,100000,0]))   

##ENERGIES
#Gaussian distribuitions
lumos = {0:(-3.7,0.0), 1:(-3.7,0.0), 'level':'lumo'} # mean energy (eV), standard deviation (eV)
homos = {0:(-6.1,0.0), 1:(-6.1,0.0), 'level':'homo'} # mean energy (eV), standard deviation (eV)
t1s   = {0:(3.7,0.0), 1:(3.7,0.0), 'level':'t1'}     # mean energy (eV), standard deviation (eV)
s1s   = {0:(6.1,0.0), 1:(6.1,0.0), 'level':'s1'}     # mean energy (eV), standard deviation (eV)

s1   = morphology.Gaussian_energy(s1s)
t1   = morphology.Gaussian_energy(t1s)
homo = morphology.Gaussian_energy(homos)
lumo = morphology.Gaussian_energy(lumos)
#########################################################################################

##GENERATE PARTICLES#####################################################################
#method    = morphology.randomized
method    = morphology.restrict_material
exciton   = morphology.Create_Particles('singlet', 2, method)

#########################################################################################

##BIMOLECULAR OPTIONS###################################################################
bimolec               = True  # Turn on annihilation
##list of all annihi funcs that will be used
bimolec_funcs_array = [morphology.ele_hol_recomb]
#########################################################################################

   

