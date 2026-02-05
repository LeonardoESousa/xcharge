############################################################
# forster_singlet.py.                                        
# This example illustrates the simulation of singlet excitons.
# 1 singlet exciton.
# 2 Materials.
############################################################

import kmc.morphology as morphology
from kmc.rates import *
from kmc.particles import *
import sys

###BASIC PARAMETERS######################################################################
K_GEN = float(sys.argv[2])
identifier         = 'generation_'+sys.argv[2] #output identifier
cutoff             = 10   # in Angstroms
#time_limit         = np.inf # in ps
rounds             = 10000   # Number of rounds
n_proc             = 10      # Number of cores to be used
frozen_lattice     = True   # if you want for the lattice to remain the same for all rounds
periodic           = False  # if you want periodic boundary conditions
bimolec            = False  # Turn on annihilation
time_limit         = 1E6
#time_limit         = 2E6
###ANIMATION SETTINGS####################################################################
square_ratio       = False
pause              = False  # if you want that the annimation stops in the first frame (debug purposes)
marker_type        = 1      # marker type used at the animation processs ( 0 = balls, 1 = symbols)
animation_mode     = False  # if you want to see the animation
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
kd = 0.0#1E-5
k_diss = {(Cs,Cs):0, (Cs,I):kd, (I,Cs):0, (I,I):0} # dissociation rates for each material combination
dissociation = DissociationFP(k=k_diss)

# FRENKEL PAIR ANNIHILATION
kaniq = 0#1E-3
k_annihilation = {Cs:kaniq, I:kaniq} # annihilation rates for each material combination
annihilation = Annihilation(k=k_annihilation)

# VACANCY/INTERSTITIAL MIGRATION
km=0
k_migration = {(Cs,Cs):km, (Cs,I):0, (I,Cs):0, (I,I):km} # migration rates for each material combination
migration = Migration(k=k_migration)

# FRENKEL PAIR FORMATION (VACANCY + INTERSTITIAL)
k_formation = 0 # formation rates for each material combination
formation = Formation(k=k_formation)

k_generation =K_GEN#{Cs:1e-8, I:0} # annihilation rates for each material combination
gen = FP_generation(k=k_generation)
###PROCESSES#############################################################################
processes = {'vacancy':[migration,formation], 'interstitial':[migration], 'frenkelpair':[dissociation]}
monomolecular = {'vacancy':[], 'interstitial':[], 'frenkelpair':[annihilation]}
generation = [gen]
#########################################################################################

###MORPHOLOGY############################################################################

##Morphology functions

# building the lattice through cif file
file = 'cifexample.cif' #filename
remove_species = ['Br']
multiply_cell  = [10,10,10] # multiplication of the unit cel by a vector v = a + b + c
lattice_func   = morphology.ReadCIF(file,multiply_cell,material_label,remove_species)





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
