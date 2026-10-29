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
time_units         = 's'
identifier         = 'lotte' #output identifier
cutoff             = 10      # in Angstroms
FPcutoff           = 5    #Frenkel pair maximum size allowed. Should not be higher than cutoff
time_limit         = np.inf # in ps
rounds             = 10      # Number of rounds
n_proc             = 10      # Number of cores to be used
frozen_lattice     = True   # if you want for the lattice to remain the same for all rounds
periodic           = False  # if you want periodic boundary conditions
bimolec            = False  # Turn on annihilation
#time_limit = 5
time_limit         = 1E3
particle_condition = False # program continues even if there are 0 particles, relevant only with generation active
###ANIMATION SETTINGS####################################################################
square_ratio       = True
pause              = False  # if you want that the annimation stops in the first frame (debug purposes)
marker_type        = 1      # marker type used at the animation processs ( 0 = balls, 1 = symbols)
animation_mode     = True  # if you want to see the animation
save_animation     = False  # if you want to save the animation
animation_exten    = 'gif'  # possible options ('gif' and 'mp4')
rotate             = True  # True = animation rotates, False = remains fixed
plot_type          = "scatter" # options: scatter or sphere
clean_vis          = True  # vizualization with no background and axis
scatter_alpha      = 0.2 # alpha of the scatter plot
sizes_dic          = {0:150, 1:50} # size of the materials in the animation mode
colors_dic         = {0:'#1FF56D', 1:'#723E9A'} # colors of the materials
print_site_position = False
#########################################################################################

T         = 298
h         = 4.135667e-15
kB        = 8.617333262e-5
Ea_ann    = 0.14
Ea_gen    = 0.8
Ea_diss   = 0.13
Ea_form   = 0.78  
Ea_mig    = 0.1
Ea_mig_I  = 0.8
Ea_mig_V  = 0.8

def arrhenius(Ea):
  return np.exp(-Ea/(kB*T))#(kB*T/h)*np.exp(-Ea/(kB*T))

k_ann   = arrhenius(Ea_ann)
k_gen   = arrhenius(Ea_gen)
k_dis   = arrhenius(Ea_diss)
k_form  = arrhenius(Ea_form)
k_mig   =  arrhenius(Ea_mig)
k_mig_I =  arrhenius(Ea_mig_I)
k_mig_V =  arrhenius(Ea_mig_V)

#k_ann = 0

print(f"k_ann={k_ann:.2e}\nk_gen={k_gen:.2e}\nk_dis={k_dis:.2e}\nk_form={k_form:.2e}\nk_mig_I={k_mig_I:.2e}\nk_mig_V={k_mig_V:.2e}")

###DEFECTS PROCESSES######################################################################
material_label  = {'Cs':0,'I':1} #site indexing
Cs         = material_label['Cs']
I          = material_label['I']


# FRENKEL PAIR DISSOCIAITON
k_diss = {Cs:k_dis, I:k_dis}
dissociation = DissociationFP(k=k_diss)

# FRENKEL PAIR ANNIHILATION
k_annihilation = {Cs:k_ann, I:k_ann} # annihilation rates for each material combination
annihilation = Annihilation(k=k_annihilation)

# VACANCY/INTERSTITIAL MIGRATION
k_migration_V = {(Cs,Cs):k_mig_V, (Cs,I):0, (I,Cs):0, (I,I):k_mig_V} # migration rates for each material combination
migration_V = Migration(k=k_migration_V)
k_migration_I = {(Cs,Cs):k_mig_I, (Cs,I):0, (I,Cs):0, (I,I):k_mig_I} # migration rates for each material combination
migration_I = Migration(k=k_migration_I)


# FRENKEL PAIR FORMATION (VACANCY + INTERSTITIAL)
k_formation = {(Cs,Cs):0, (Cs,I):k_form, (I,Cs):k_form, (I,I):0}
formation = Formation(k=k_formation)

generation_pairs = [[I,Cs]]
gen = FP_generation(k=k_gen,pairs=generation_pairs)
###PROCESSES#############################################################################

processes = {'vacancy':[migration_V, formation], 'interstitial':[migration_I], 'frenkelpair':[]}
monomolecular = {'vacancy':[], 'interstitial':[], 'frenkelpair':[annihilation,dissociation]}
generation = [gen]
#########################################################################################

###MORPHOLOGY############################################################################

##Morphology functions

# building the lattice through cif file
file = 'CsPbI3_bulk.cif' #filename
remove_species = ['Br','Pb']
multiply_cell  = [2,2,2] # multiplication of the unit cel by a vector v = a + b + c
lattice_func   = morphology.ReadCIF(file,multiply_cell,material_label,remove_species)





#ENERGIES
#Gaussian distribuitions
t1s   = {0:(3.7,0.0), 1:(3.7,0.0), 'level':'t1'} # Peak emission energy (eV), disperison (eV)
s1s   = {0:(6.1,0.0), 1:(6.1,0.0), 'level':'s1'} # triplet energy, dispersion (eV)

a1 = morphology.Gaussian_energy(s1s)
a2 = morphology.Gaussian_energy(t1s) 
#########################################################################################


##GENERATE PARTICLES#####################################################################
method  = morphology.randomized
#vac   = morphology.Create_Particles('vacancy', 3, method, mat=[I])
#inter   = morphology.Create_Particles('interstitial', 3, method, mat=[Cs])
FP      = morphology.Create_ParticlesFP(num=10,pairs=generation_pairs)
#########################################################################################
