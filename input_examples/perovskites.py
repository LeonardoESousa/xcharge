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
identifier         = 'perovskites' #output identifier
cutoff             = 10   # in Angstroms
time_limit         = np.inf # in ps
animation_mode     = True  # if you want to see the animation
save_animation     = False  # if you want to save the animation
animation_exten    = 'gif'  # possible options ('gif' and 'mp4')
rotate             = False  # True = animation rotates, False = remains fixed
marker_type        = 1      # marker type used at the animation processs ( 0 = balls, 1 = symbols) 
pause              = False  # if you want that the annimation stops in the first frame (debug purposes)
rounds             = 1   # Number of rounds
n_proc             = 1      # Number of cores to be used
frozen_lattice     = True   # if you want for the lattice to remain the same for all rounds
periodic           = False  # if you want periodic boundary conditions
bimolec            = False  # Turn on annihilation
#########################################################################################

###DEFECTS PROCESSES######################################################################

Cs = 0
I  = 1


# FRENKEL PAIR DISSOCIAITON
k_diss = {(Cs,Cs):0, (Cs,I):1e-5, (I,Cs):0, (I,I):0} # dissociation rates for each material combination
dissociation = DissociationFP(k=k_diss)

# FRENKEL PAIR ANNIHILATION
k_annihilation = {Cs:1e-8, I:0} # annihilation rates for each material combination
annihilation = Annihilation(k=k_annihilation)

# VACANCY/INTERSTITIAL MIGRATION
k_migration = {(Cs,Cs):1e-5, (Cs,I):0, (I,Cs):0, (I,I):1e-5} # migration rates for each material combination
migration = Migration(k=k_migration)

# FRENKEL PAIR FORMATION (VACANCY + INTERSTITIAL)
k_formation = 1e-4 # formation rates for each material combination
formation = Formation(k=k_formation)

###PROCESSES#############################################################################

processes = {'vacancy':[migration, formation], 'interstitial':[migration], 'frenkelpair':[dissociation]}
monomolecular = {'vacancy':[], 'interstitial':[], 'frenkelpair':[annihilation]}
#########################################################################################

###MORPHOLOGY############################################################################

##Morphology functions

#Reading a file name that contains your lattice
#file = 'lattice.example'
#lattice_func = morphology.ReadLattice(file)


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

a1 = morphology.Gaussian_energy(s1s)
a2 = morphology.Gaussian_energy(t1s) 
#########################################################################################


##GENERATE PARTICLES#####################################################################
method    = morphology.randomized
exciton   = morphology.Create_Particles('frenkelpair', 1, method, mat=[0]) # creates 1 singlet exciton randomly at either material 0 or 1
#########################################################################################
