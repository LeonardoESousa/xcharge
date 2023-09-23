############################################################
# forster_singlet.py.
# This example illustrates the simulation of singlet excitons.
# 1 singlet exciton.
# 2 Materials.
############################################################

import kmc.morphology as morphology
from kmc.rates import *
from kmc.particles import *

###BASIC PARAMETERS######################################################################
identifier         = 'forster_singlet' #output identifier
cutoff             = 20     # cutoff distance for rates (Å)
time_limit         = np.inf # in ps
animation_mode     = True   # if you want to see the animation
save_animation     = False  # if you want to save the animation
animation_exten    = 'gif'  # possible options ('gif' and 'mp4')
rotate             = False  # True = animation rotates, False = remains fixed
marker_type        = 1      # marker type used at the animation processs ( 0 = balls, 1 = symbols) 
pause              = False  # if you want that the annimation stops in the first frame (debug purposes)
rounds             = 5000   # Number of rounds
n_proc             = 1      # Number of cores to be used
frozen_lattice     = True   # if you want for the lattice to remain the same for all rounds
periodic           = False  # if you want periodic boundary conditions
bimolec            = False  # Turn on annihilation
#########################################################################################

###SINGLET EXCITONS######################################################################

##FORSTER RADII (Å)
r00   = 25  #Forster radius material 0 --> material 0 (Angstrom)
r01   = 0   #material 0 --> material 1
r10   = 0
r11   = 25
radii = {(0,0):r00, (0,1):r01, (1,0):r10, (1,1):r11}

##FLUORESCENCE LIFETIMES (PS)
f0 = 2900 #lifetime of material 0
f1 = 290  #lifetime of material 1
lifetimes = {0:f0,1:f1}

##TANSITION DIPOLE MOMENTS (a.u.)
mu0 = 2.136
mu1 = 5.543
mus = {0:mu0,1:mu1}

##EXCITION TRANSFER RATES
forster   = Forster(Rf=radii,life=lifetimes,mu=mus)

##FLUORESCENCE RATES
fluor     = Fluor(life=lifetimes)

###PROCESSES#############################################################################

processes = {'singlet':[forster], 'triplet':[forster], 'electron':[],'hole':[]}
monomolecular = {'singlet':[fluor],'triplet':[],'electron':[],'hole':[]}
#########################################################################################

###MORPHOLOGY############################################################################

# Creating a new lattice at each new round
num_sites         = 100             #number of sites of the lattice
displacement      = [5, 5, 0]       #vector of the unit cell (x,y,z)
disorder          = [0.,0.,0.]      #std deviation from avg position
composition       = [0.5,0.5]       #population probability Ex.: composition[0] is the prob of mat 0 appear in the lattice
lattice_func      = morphology.Lattice(num_sites,displacement,disorder,composition)

#ENERGIES
#Gaussian distribuitions
t1s   = {0:(3.7,0.0), 1:(3.7,0.0), 'level':'t1'} # Peak emission energy (eV), disperison (eV)
s1s   = {0:(6.1,0.0), 1:(6.1,0.0), 'level':'s1'} # triplet energy, dispersion (eV)

a1 = morphology.GaussianEnergy(s1s)
a2 = morphology.GaussianEnergy(t1s)
#########################################################################################


##GENERATE PARTICLES#####################################################################
method    = morphology.randomized
exciton   = morphology.CreateParticles(['singlet','triplet'], [1,3], 5, method, mat=[0]) # creates 1 singlet exciton randomly at either material 0 or 1
#########################################################################################
