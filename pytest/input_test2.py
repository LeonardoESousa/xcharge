import kmc.morphology as morphology
from kmc.rates import *
from kmc.particles import *
import sys
###BASIC PARAMETERS######################################################################
identifier         = 'test2_'+(sys.argv[2]) #output identifier
time_limit         = np.inf# in PS
animation_mode     = False
save_animation     = False # if you want to save the animation
animation_exten    = 'gif' # possible options ('gif' and 'mp4')
rotate             = False # True = animation rotates, False = remains fixed
marker_type        = 1     # marker type used at the animation processs ( 0 = balls, 1 = symbols) 
pause              = False # if you want that the annimation stops in the first frame (debug purposes)
rounds             = 5000  # Number of rounds
n_proc             = 5     # Number of cores to be used
frozen             = True
bimolec            = False # Turn on annihilation
periodic           = False
#########################################################################################

###MORPHOLOGY############################################################################

##Morphology functions

#Reading a file name that contains your lattice
#file = 'lattice.example'
#lattice_func = morphology.ReadLattice(file)


# Creating a new lattice at each new round
num_sites         = 100*100        #number of sites of the lattice
displacement      = [15, 15, 0]    #vector of the unit cell
disorder          = [0.0,0.0,0.0]   #std deviation from avg position
composition       = [0.5,0.5]       #popuation probility Ex.: distribu_vect[0] is the prob of mat 0 appear in the lattice
lattice_func      = morphology.Lattice(num_sites,displacement,disorder,composition)

#ENERGIES
#Gaussian distribuitions
t1s   = {0:(1.1,0.0), 1:(1.7,0.0), 'level':'t1'} #(Peak emission energy (eV), disperison (eV)
s1s   = {0:(6.1,0.0), 1:(3.7,0.0), 'level':'s1'} # triplet energy, disperison (eV)

a1 = morphology.Gaussian_energy(s1s)
a2 = morphology.Gaussian_energy(t1s) 
#########################################################################################


##GENERATE PARTICLES#####################################################################
method    = morphology.randomized
exciton   = morphology.Create_Particles('singlet', int(sys.argv[2]), method, mat=[0,1])

#########################################################################################


###SINGLET EXCITONS######################################################################

##FORSTER RADII (Å)
r00   = 20   #Forster radius material 0 --> material 0 (Angstrom)    
r01   = 20    #material 0 --> material 1      
r10   = 10     
r11   = 10    
radii = {(0,0):r00, (0,1):r01, (1,0):r10, (1,1):r11}

##FLUORESCENCE LIFETIMES (PS)
f0 = 3000 #lifetime of material 0
f1 = 6000 #lifetime of material 1
lifetimes = {0:f0,1:f1}

##TANSITION DIPOLE MOMENTS (a.u.)
mu0 = 0
mu1 = 0
mus = {0:mu0,1:mu1}

##EXCITION TRANSFER RATES
forster   = Forster(Rf=radii,life=lifetimes,mu=mus)

##FLUORESCENCE RATES
fluor     = Fluor(life=lifetimes)


#isc = ISC(rate={0:1e-4,1:1e-4})
#nonrad = Nonrad(rate={0:1e-30,1:1e-30})
###PROCESSES#############################################################################

processes = {'singlet':[forster], 'triplet':[], 'electron':[],'hole':[]}
monomolecular = {'singlet':[fluor],'triplet':[],'electron':[],'hole':[]}
#########################################################################################



   

