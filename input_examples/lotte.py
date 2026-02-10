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
identifier         = 'lotte' #output identifier
cutoff             = 10   # in Angstroms
time_limit         = np.inf # in ps
rounds             = 1   # Number of rounds
n_proc             = 1      # Number of cores to be used
frozen_lattice     = True   # if you want for the lattice to remain the same for all rounds
periodic           = False  # if you want periodic boundary conditions
bimolec            = False  # Turn on annihilation
time_limit = 1E8
particle_condition = False # program continues even if there are 0 particles, relevant only with generation active
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

T        = 298
h        = 4.135667e-15
kB       = 8.617333262e-5
Ea_ann   = 0.14
Ea_gen   = 0.8
Ea_diss  = 0.13
Ea_form  = 0.78  

def arrhenius(Ea):
  return (kB*T/h)*np.exp(-Ea/(kB*T))

k_ann = arrhenius(Ea_ann)
k_gen = arrhenius(Ea_gen)
k_dis = arrhenius(Ea_diss)
k_form = arrhenius(Ea_form)


k_ann  = 1E-6
k_gen  = 1E-6#1E-5
k_dis  = 1E-2
k_form = 1E-6
k_mig  = 1E-5
print(k_ann,k_gen,k_dis,k_form)

###DEFECTS PROCESSES######################################################################
material_label  = {'Cs':0,'I':1} #site indexing
Cs         = material_label['Cs']
I          = material_label['I']


# FRENKEL PAIR DISSOCIAITON
k_diss = {(Cs,Cs):0, (Cs,I):k_dis, (I,Cs):0, (I,I):0} # dissociation rates for each material combination
dissociation = DissociationFP(k=k_diss)

# FRENKEL PAIR ANNIHILATION
k_annihilation = {Cs:k_ann, I:k_ann} # annihilation rates for each material combination
annihilation = Annihilation(k=k_annihilation)

# VACANCY/INTERSTITIAL MIGRATION
k_migration = {(Cs,Cs):k_mig, (Cs,I):0, (I,Cs):0, (I,I):k_mig} # migration rates for each material combination
migration = Migration(k=k_migration)

# FRENKEL PAIR FORMATION (VACANCY + INTERSTITIAL)
formation = Formation(k=k_form)

gen = FP_generation(k=k_gen)
###PROCESSES#############################################################################

processes = {'vacancy':[migration, formation], 'interstitial':[migration], 'frenkelpair':[dissociation]}
monomolecular = {'vacancy':[], 'interstitial':[], 'frenkelpair':[annihilation]}
generation = [gen]
#########################################################################################

###MORPHOLOGY############################################################################

##Morphology functions

# building the lattice through cif file
file = 'cifexample.cif' #filename
remove_species = ['Br']
multiply_cell  = [3,3,4] # multiplication of the unit cel by a vector v = a + b + c
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
