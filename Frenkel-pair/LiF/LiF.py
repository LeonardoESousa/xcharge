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
identifier          = 'LiF' #output identifier
cutoff              = 2.9       # in Angstroms
FPcutoff            = 2.9     #Frenkel pair maximum size allowed. Should not be higher than cutoff
rounds              = 500      # Number of rounds
n_proc              = 10      # Number of cores to be used
frozen_lattice      = True   # if you want for the lattice to remain the same for all rounds
periodic            = True  # if you want periodic boundary conditions
bimolec             = False  # Turn on annihilation
time_limit          = 1E-4#1E-2
particle_condition  = False # program continues even if there are 0 particles, relevant only with generation active
print_site_position = True
###ANIMATION SETTINGS####################################################################
square_ratio       = True
pause              = True  # if you want that the annimation stops in the first frame (debug purposes)
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

T        = 298
h        = 4.135667e-15
kB       = 8.617333262e-5
Ea_ann   = 0.395
Ea_gen   = 0.0
Ea_diss  = 0.0
Ea_form  = 0.0  
Ea_mig_INTER = 0.1
Ea_mig_VACAN = 1
def arrhenius(Ea):
  return (kB*T/h)*np.exp(-Ea/(kB*T))

k_ann       = 0.0 #arrhenius(Ea_ann)
k_gen       = 0.0#arrhenius(Ea_gen)
k_dis       = 0.0#arrhenius(Ea_diss)
k_form      = 0.0#arrhenius(Ea_form)
k_mig       = 1E-3
k_mig_INTER = arrhenius(Ea_mig_INTER)
k_mig_VACAN = arrhenius(Ea_mig_VACAN)

print(k_ann,k_gen,k_dis,k_form)
print(k_mig_INTER,k_mig_VACAN)

k_mig_INTER = 2*1E5
k_mig_VACAN = 2*1E5

###DEFECTS PROCESSES######################################################################
material_label  = {'Li':0,'F':1} #site indexing
Li         = material_label['Li']
F          = material_label['F']


# FRENKEL PAIR DISSOCIAITON
k_diss = {Li:k_dis, F:k_dis}
dissociation = DissociationFP(k=k_diss)

# FRENKEL PAIR ANNIHILATION
k_annihilation = {Li:k_ann, F:k_ann}
annihilation = Annihilation(k=k_annihilation)

# VACANCY/INTERSTITIAL MIGRATION
k_migration_inter = {(Li,Li):k_mig_INTER, (Li,F):0, (F,Li):0, (F,F):k_mig_INTER} 
migration_inter  = Migration(k=k_migration_inter)

k_migration_vac = {(Li,Li):k_mig_VACAN, (Li,F):0, (F,Li):0, (F,F):k_mig_VACAN} 
migration_vac   = Migration(k=k_migration_vac )
# FRENKEL PAIR FORMATION (VACANCY + INTERSTITIAL)
k_formation = {(Li,Li):0, (Li,F):0, (F,Li):0, (F,F):0}
formation = Formation(k=k_formation)


#convention:
#pair = [x,y] -> the frenkelpair is the combination of a vacacy at site of type x with a intersitial of y
generation_pairs = [[F,Li]]
gen = FP_generation(k=k_gen,pairs=generation_pairs)
###PROCESSES#############################################################################

#processes = {'vacancy':[migration_vac, formation], 'interstitial':[migration_inter], 'frenkelpair':[]}
#monomolecular = {'vacancy':[], 'interstitial':[], 'frenkelpair':[annihilation,dissociation]}
#generation = [gen]

processes = {'vacancy':[migration_vac], 'interstitial':[migration_inter], 'frenkelpair':[]}
monomolecular = {'vacancy':[], 'interstitial':[], 'frenkelpair':[]}
generation = [gen]
#########################################################################################

###MORPHOLOGY############################################################################

##Morphology functions

# building the lattice through cif file
file = 'LiF.cif' #filename
remove_species = []
multiply_cell  = [3,3,3] # multiplication of the unit cel by a vector v = a + b + c
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
#Vac       = morphology.Create_Particles('vacancy', 1, method, mat=[Li])
Inter     = morphology.Create_Particles('interstitial', 1, method, mat=[F])

#########################################################################################
