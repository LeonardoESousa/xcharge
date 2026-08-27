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
rounds             = 1   # Number of rounds
n_proc             = 1      # Number of cores to be used
frozen_lattice     = True   # if you want for the lattice to remain the same for all rounds
periodic           = False  # if you want periodic boundary conditions
bimolec            = False  # Turn on annihilation
time_units         = 's'
first_neighbor_rtol = 0.15
###ANIMATION SETTINGS####################################################################
square_ratio       = False
pause              = False  # if you want that the annimation stops in the first frame (debug purposes)
marker_type        = 1      # marker type used at the animation processs ( 0 = balls, 1 = symbols)
animation_mode     = True  # if you want to see the animation
save_animation     = False  # if you want to save the animation
animation_exten    = 'gif'  # possible options ('gif' and 'mp4')
rotate             = False  # True = animation rotates, False = remains fixed
plot_type          = "scatter" # options: scatter or sphere
clean_vis          = False  # vizualization with no background and axis
scatter_alpha      = 0.5 # alpha of the scatter plot
sizes_dic          = {0:150, 1:50} # size of the materials in the animation mode
colors_dic         = {0:'#1FF56D', 1:'#723E9A'} # colors of the materials
#########################################################################################

###DEFECTS PROCESSES######################################################################
material_label  = {'Cs':0,'I':1} #site indexing
Cs         = material_label['Cs']
I          = material_label['I']


# FRENKEL PAIR DISSOCIAITON
k_diss = {Cs:0, I:1e-5}
dissociation = DissociationFP(k=k_diss)

# FRENKEL PAIR ANNIHILATION
k_annihilation = {Cs:1e-8, I:0} # annihilation rates for each material combination
annihilation = Annihilation(k=k_annihilation)

# VACANCY/INTERSTITIAL MIGRATION
k_migration = {(Cs,Cs):1e-5, (Cs,I):0, (I,Cs):0, (I,I):1e-5} # migration rates for each material combination
migration = Migration(k=k_migration)

# FRENKEL PAIR FORMATION (VACANCY + INTERSTITIAL)
k_formation = {(Cs,Cs):0, (Cs,I):1e-4, (I,Cs):1e-4, (I,I):0}
formation = Formation(k=k_formation)

###PROCESSES#############################################################################

processes = {'vacancy':[migration, formation], 'interstitial':[migration, formation], 'frenkelpair2':[]}
monomolecular = {'vacancy':[], 'interstitial':[], 'frenkelpair2':[annihilation, dissociation]}
#########################################################################################

###MORPHOLOGY############################################################################

##Morphology functions

# building the lattice through cif file
file = 'cifexample.cif' #filename
remove_species = ['Br']
multiply_cell  = [2,3,4] # multiplication of the unit cel by a vector v = a + b + c
lattice_func   = morphology.ReadCIF(file,multiply_cell,material_label,remove_species)





#ENERGIES
#Gaussian distribuitions
t1s   = {0:(3.7,0.0), 1:(3.7,0.0), 'level':'t1'} # Peak emission energy (eV), disperison (eV)
s1s   = {0:(6.1,0.0), 1:(6.1,0.0), 'level':'s1'} # triplet energy, dispersion (eV)

a1 = morphology.Gaussian_energy(s1s)
a2 = morphology.Gaussian_energy(t1s) 
#########################################################################################


##GENERATE PARTICLES#####################################################################
generation_pairs = [[I, Cs]]
exciton = morphology.Create_ParticlesFP(num=2, pairs=generation_pairs)
#########################################################################################
