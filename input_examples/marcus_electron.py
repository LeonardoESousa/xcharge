############################################################
# marcus_electron.py.
# This example illustrates the simulation of charge transfers with Marcus rate.
# 1 electron.
# 2 Materials.
############################################################
import numpy as np
import kmc.morphology as morphology
from kmc.rates import Marcus

###BASIC PARAMETERS######################################################################
identifier         = 'marcus_electron' #output identifier
cutoff             = 6 # cutoff distance for rates (Å)
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

###MORPHOLOGY############################################################################

# Lattice parameters
num_sites         = 100             # number of sites in the lattice
displacement      = [5, 5, 0]       # vector of the unit cell (x,y,z) (Å)
disorder          = [0.,0.,0.]      # std deviation from avg position (Å)
composition       = [0.5,0.5]       # population probability Ex.: composition[0] is the prob of mat 0 appear in the lattice
lattice_func      = morphology.Lattice(num_sites,displacement,disorder,composition)

# Energy levels
#Make LUMO values for electrons
level  = 'lumo'                                    # energy level used in the transfer rate
lumo   = {0:(3.7,0.1), 1:(3.7,0.1), 'level':level} # {material:(mean (eV), std deviation (eV)),'level':identify the level}
a1 = morphology.GaussianEnergy(lumo)
#########################################################################################


##GENERATE PARTICLES#####################################################################
particles = ['electron'] # list of particles to be added
probs     = [1]          # list of probabilities for each particle
mat       = [0]          # list of materials where particle may appear
number    = 1            # number of particles to be added
method    = morphology.randomized # method to generate particles
particle  = morphology.CreateParticles(particles, probs, number, method, mat=mat)
#########################################################################################


##TRANSFER RATES#########################################################################
coupling = {(0,0):0.1, (0,1):0.1, (1,0):0.1, (1,1):0.1} # {(material_donor,material_acceptor):coupling (eV)}
reorg = {0:0.1,1:0.1}                                   # {material:reorganization energy (eV)}
level  = 'lumo'                                         # energy level used in the transfer rate
temperature = 300                                       # in K
decay = 0.5                                             # coupling exponential decay constant in 1/Å
marcus = Marcus(coupling=coupling,reorg=reorg,level=level, temperature=temperature,decay=decay)

###PROCESSES#############################################################################
# Assign processes to each particle
# transfer processes
processes = {'singlet':[], 'triplet':[], 'electron':[marcus],'hole':[]}
# monomolecular processes
monomolecular = {'singlet':[],'triplet':[],'electron':[],'hole':[]}
#########################################################################################
