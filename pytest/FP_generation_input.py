############################################################
# FP_generation.py.                                        
# This example tests the frenkel pair generation algorithm
############################################################

import kmc.morphology as morphology
from kmc.rates import *
from kmc.particles import *
import sys

###BASIC PARAMETERS######################################################################
K_GEN = float(sys.argv[2])
identifier         = 'generation_'+sys.argv[2] #output identifier
cutoff             = 10      # in Angstroms
FPcutoff           = 10    #Frenkel pair maximum size allowed. Should not be higher than cutoff
rounds             = 10000      # Number of rounds
n_proc             = 10      # Number of cores to be used
frozen_lattice     = True   # if you want for the lattice to remain the same for all rounds
periodic           = False  # if you want periodic boundary conditions
bimolec            = False  # Turn on annihilation
time_limit = 1E6
particle_condition = False # program continues even if there are 0 particles, relevant only with generation active
print_site_position = True
###ANIMATION SETTINGS####################################################################
square_ratio       = False
pause              = False  # if you want that the annimation stops in the first frame (debug purposes)
marker_type        = 1      # marker type used at the animation processs ( 0 = balls, 1 = symbols)
animation_mode     = False  # if you want to see the animation
save_animation     = False  # if you want to save the animation
animation_exten    = 'gif'  # possible options ('gif' and 'mp4')
rotate             = False  # True = animation rotates, False = remains fixed
plot_type          = "scatter" # options: scatter or sphere
clean_vis          = True  # vizualization with no background and axis
scatter_alpha      = 0.2 # alpha of the scatter plot
sizes_dic          = {0:150, 1:50} # size of the materials in the animation mode
colors_dic         = {0:'#1FF56D', 1:'#723E9A'} # colors of the materials
#########################################################################################

k_ann  = 0
k_gen  = K_GEN
k_dis  = 0
k_form = 0
k_mig  = 0
#print(k_ann,k_gen,k_dis,k_form)

###DEFECTS PROCESSES######################################################################
material_label  = {'Cs':0,'I':1} #site indexing
Cs         = material_label['Cs']
I          = material_label['I']


# FRENKEL PAIR DISSOCIAITON
k_diss = {Cs:k_dis, I:k_dis}#{(Cs,Cs):0, (Cs,I):k_dis, (I,Cs):0, (I,I):0} # dissociation rates for each material combination
dissociation = DissociationFP(k=k_diss)

# FRENKEL PAIR ANNIHILATION
k_annihilation = {Cs:k_ann, I:k_ann} # annihilation rates for each material combination
annihilation = Annihilation(k=k_annihilation)

# VACANCY/INTERSTITIAL MIGRATION
k_migration = {(Cs,Cs):k_mig, (Cs,I):0, (I,Cs):0, (I,I):k_mig} # migration rates for each material combination
migration = Migration(k=k_migration)

# FRENKEL PAIR FORMATION (VACANCY + INTERSTITIAL)
k_formation = {(Cs,Cs):0, (Cs,I):k_form, (I,Cs):k_form, (I,I):0}
formation = Formation(k=k_formation)


#convention:
#pair = [x,y] -> the frenkelpair is the combination of a vacacy at site of type x with a intersitial of y
generation_pairs = [[I,Cs]]
gen = FP_generation(k=k_gen,pairs=generation_pairs)
###PROCESSES#############################################################################

processes = {'vacancy':[migration, formation], 'interstitial':[migration], 'frenkelpair0':[]}
monomolecular = {'vacancy':[], 'interstitial':[], 'frenkelpair0':[annihilation,dissociation]}
generation = [gen]
#########################################################################################

###MORPHOLOGY############################################################################

##Morphology functions

# building the lattice through cif file
file = "pytest/cifexample.cif" #filename
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
FP   = morphology.Create_ParticlesFP(num=1,pairs=generation_pairs)
#########################################################################################
