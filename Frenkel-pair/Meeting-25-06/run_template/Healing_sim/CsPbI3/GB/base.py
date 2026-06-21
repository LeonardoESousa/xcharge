############################################################
# FP_generation.py.                                        
# This example tests the frenkel pair generation algorithm
############################################################
from itertools import product
import kmc.morphology as morphology
from kmc.rates import *
from kmc.particles import *
import sys

nFP = int(sys.argv[-1])
###BASIC PARAMETERS######################################################################
identifier         = 'CsPbI3_GB' #output identifier
cutoff             = 6     # in Angstroms
FPcutoff           = 6    #Frenkel pair maximum size allowed. Should not be higher than cutoff
rounds             = 100      # Number of rounds
n_proc             = 10# Number of cores to be used
frozen_lattice     = True   # if you want for the lattice to remain the same for all rounds
periodic           = True  # if you want periodic boundary conditions
bimolec            = False  # Turn on annihilation
time_limit         = 1E-6
particle_condition  = True # program continues even if there are 0 particles, relevant only with generation active
creation_threshold  = 1E+1
###ANIMATION SETTINGS####################################################################
print_site_position= False
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
sizes_dic          = {0:150, 1:50,2:150, 3:50} # size of the materials in the animation mode
colors_dic         = {0:'#1FF56D', 1:'#723E9A',2:"#FF0000",3:"#00FFF7"} # colors of the materials
#########################################################################################
T         = 298
h         = 4.135667e-15
kB        = 8.617333262e-5
##############################
#Ea_gen        = 0.11

Ea_ann_FP0    = 0.11
Ea_ann_FP2    = 0.19
Ea_ann_FP0_GB = 0.09
Ea_ann_FP2_GB = 0.55


Ea_diss_V     = 0.2
Ea_form_V     = 0.3 

Ea_diss_V_GB  = 0.1
Ea_form_V_GB  = 0.26

Ea_diss_I     = 0.48
Ea_form_I     = 0.25 

Ea_diss_I_GB  = 0.48 #mirrowing bulk because we dont have the data
Ea_form_I_GB  = 0.25 #mirrowing bulk because we dont have the data 

Ea_mig_V      = [0.31,0.57,0.34,0.76]
Ea_mig_I      = 0.26
Ea_mig_V_GB1  = [0.27,0.13,0.21,0.26]
Ea_mig_V_GB2  = [0.5,0.71,0.1,0.23]
Ea_mig_I_GB1  = [0.27,0.13,0.21,0.26]


def arrhenius(Ea):
  return (kB*T/h)*np.exp(-Ea/(kB*T))#(kB*T/h)*np.exp(-Ea/(kB*T))
def seq_rates(Ea_list):
  return 1/sum([1/arrhenius(E) for E in Ea_list])

k_gen   = 0 #arrhenius(Ea_gen)

k_ann_FP0      = arrhenius(Ea_ann_FP0)
k_ann_FP2      = arrhenius(Ea_ann_FP2)
k_ann_FP0_GB   = arrhenius(Ea_ann_FP0_GB)
k_ann_FP2_GB   = arrhenius(Ea_ann_FP2_GB)

k_dis    = arrhenius(Ea_diss_V)+arrhenius(Ea_diss_I)      #independent events -> sum
k_dis_GB = arrhenius(Ea_diss_V_GB)+arrhenius(Ea_diss_I_GB)#independent events -> sum

k_mig_I    = arrhenius(Ea_mig_I)
k_mig_I_GB = seq_rates(Ea_mig_I_GB1) #+ seq_rates(Ea_mig_I_GB2)
k_mig_V    = seq_rates(Ea_mig_V) #effective rate for sequential event
k_mig_V_GB = seq_rates(Ea_mig_V_GB1) + seq_rates(Ea_mig_V_GB2)

k_form_I = arrhenius(Ea_form_I)
k_form_V = arrhenius(Ea_form_V)
k_form_I_GB = arrhenius(Ea_form_I_GB)
k_form_V_GB = arrhenius(Ea_form_V_GB)


###DEFECTS PROCESSES######################################################################
material_label  = {'Cs':0,'I':1,'Cs_GB':2,'I_GB':3} #site indexing
Cs         = material_label['Cs']
I          = material_label['I']
Cs_GB      = material_label['Cs_GB']
I_GB       = material_label['I_GB']

# FRENKEL PAIR DISSOCIAITON
k_diss = {Cs:k_dis, Cs_GB:k_dis_GB, I:k_dis, I_GB:k_dis_GB} # dissociation rates for each material combination
dissociation = DissociationFP(k=k_diss)

# FRENKEL PAIR ANNIHILATION
k_annihilation = {Cs:k_ann_FP0, Cs_GB:k_ann_FP0_GB, I:k_ann_FP0, I_GB:k_ann_FP0_GB} # annihilation rates for each material combination
annihilation_FP0 = Annihilation(k=k_annihilation)

k_annihilation = {Cs:k_ann_FP2, Cs_GB:k_ann_FP2_GB, I:k_ann_FP2, I_GB:k_ann_FP2_GB} # annihilation rates for each material combination
annihilation_FP2 = Annihilation(k=k_annihilation)


def dictionary_pairs(x,y):
  kdict = {}
  for i, j in product(x,y):
    kdict[(i,j)] = 0
  return kdict
# VACANCY/INTERSTITIAL MIGRATION


kmigration_I = dictionary_pairs(material_label.values(), material_label.values())
kmigration_I.update({
    (Cs, Cs): k_mig_I,
    (I, I): k_mig_I,

    (Cs, Cs_GB): k_mig_I_GB,
    (Cs_GB, Cs_GB): k_mig_I_GB,
    (Cs_GB, Cs): k_mig_I_GB,

    (I, I_GB): k_mig_I_GB,
    (I_GB, I): k_mig_I_GB,
    (I_GB, I_GB): k_mig_I_GB}
    
    )
migration_I = Migration(k=kmigration_I)


kmigration_V = dictionary_pairs(material_label.values(), material_label.values())
kmigration_V.update({
    (Cs, Cs): k_mig_V,
    (I, I): k_mig_V,

    (Cs, Cs_GB): k_mig_V_GB,
    (Cs_GB, Cs_GB): k_mig_V_GB,
    (Cs_GB, Cs): k_mig_V_GB,
    
    (I, I_GB): k_mig_V_GB,
    (I_GB, I): k_mig_V_GB,
    (I_GB, I_GB): k_mig_V_GB})
migration_V = Migration(k=kmigration_V)
###################
# FRENKEL PAIR FORMATION (VACANCY + INTERSTITIAL)
kformation_I = dictionary_pairs(material_label.values(), material_label.values())
kformation_I.update({
   (Cs,I):k_form_I,
   (I,Cs):k_form_I,
   (I_GB,Cs):k_form_I_GB,
   (I,Cs_GB):k_form_I_GB,
   (Cs_GB,I_GB):k_form_I_GB})

kformation_V = dictionary_pairs(material_label.values(), material_label.values())
kformation_V.update({
   (Cs,I):k_form_V,
   (I,Cs):k_form_V,
   (I_GB,Cs):k_form_V,
   (I,Cs_GB):k_form_V,
   (Cs_GB,I_GB):k_form_V})

formation_I = Formation(k=kformation_I)
formation_V = Formation(k=kformation_V)

#convention:
#pair = [x,y] -> the frenkelpair is the combination of a vacacy at site of type x with a intersitial of y
generation_pairs = [[I,Cs],[I,Cs_GB],[I_GB,Cs],[I_GB,Cs_GB]]
gen = FP_generation(k=k_gen,pairs=generation_pairs)
###PROCESSES#############################################################################

processes = {'vacancy':[migration_V, formation_V], 'interstitial':[migration_I,formation_I], 'frenkelpair0':[], 'frenkelpair2':[]}
monomolecular = {'vacancy':[], 'interstitial':[], 'frenkelpair0':[annihilation_FP0], 'frenkelpair2':[annihilation_FP2,dissociation]}
generation = [gen]
#########################################################################################

###MORPHOLOGY############################################################################

##Morphology functions

# building the lattice through cif file
#file = "../CsPbI3_bulk.cif" #filename
#file = "../../CsPbI3_bulk.cif" #filename
file = "../../CsPbI3_KMC_unitcell.cif" #filename
remove_species = ['Pb']
multiply_cell  = [1,1,1]#[4,4,4] # multiplication of the unit cel by a vector v = a + b + c
lattice_func   = morphology.ReadCIF(file,multiply_cell,material_label,remove_species)





#ENERGIES
#Gaussian distribuitions
t1s   = {0:(3.7,0.0), 1:(3.7,0.0),2:(3.7,0.0),3:(3.7,0.0), 'level':'t1'} # Peak emission energy (eV), disperison (eV)
s1s   = {0:(6.1,0.0), 1:(6.1,0.0),2:(3.7,0.0),3:(3.7,0.0), 'level':'s1'} # triplet energy, dispersion (eV)

a1 = morphology.Gaussian_energy(s1s)
a2 = morphology.Gaussian_energy(t1s) 
#########################################################################################


##GENERATE PARTICLES#####################################################################
method    = morphology.randomized
#Vac       = morphology.Create_Particles('vacancy', 1, method, mat=[I])
#Int       = morphology.Create_Particles('interstitial', 1, method, mat=[Cs_GB])
FP   = morphology.Create_ParticlesFP(num=nFP,pairs=generation_pairs) #ONLY CREATES FP2+
#########################################################################################
