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
identifier         = 'marcus_electron' #output identifier
cutoff             = np.inf     # cutoff distance for rates (Å)
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


coupling = {(0,0):0.1, (0,1):0.1, (1,0):0.1, (1,1):0.1} # in eV
reorg = {0:0.1,1:0.1} # reorganization energy in eV
decay = 0.5 # coupling exponential decay constant in 1/Å
marcus = Marcus(coupling=coupling,reorg=reorg,level='lumo', temperature=300,decay=decay)

###PROCESSES#############################################################################

processes = {'singlet':[], 'triplet':[], 'electron':[marcus],'hole':[]}
monomolecular = {'singlet':[],'triplet':[],'electron':[],'hole':[]}
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
lumo   = {0:(3.7,0.1), 1:(3.7,0.1), 'level':'lumo'} # Peak emission energy (eV), disperison (eV)

a1 = morphology.GaussianEnergy(lumo)
#########################################################################################


##GENERATE PARTICLES#####################################################################
method    = morphology.randomized
exciton   = morphology.CreateParticles(['electron'], [1], 1, method, mat=[0]) # creates 1 singlet exciton randomly at either material 0 or 1
#########################################################################################
