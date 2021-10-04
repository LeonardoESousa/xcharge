import kmc.morphology as morphology
from kmc.kmc_classes import *

###BASIC PARAMETERS######################################################################
identifier         = 'New' #output identifier
time_limit         = np.inf
animation_mode     = False  
pause              = False # if you want that the annimation stops in the first frame (debug purposes)
rounds             = 100   # Number of rounds
n_proc             = 6     # Number of cores to be used
#########################################################################################

###SINGLET EXCITONS######################################################################

##FORSTER RADII (Å)
r00   = 25   #Forster radius material 0 --> material 0 (Angstrom)    
r01   = 25   #material 0 --> material 1      
r10   = 25       
r11   = 25     
raios = {(0,0):r00, (0,1):r01, (1,0):r10, (1,1):r11}

##FLUORESCENCE LIFETIMES (PS)
f0 = 2900000 #lifetime of material 0
f1 = 2900 #lifetime of material 1
lifetimes = {0:f0,1:f1}

##TANSITION DIPOLE MOMENTS (a.u.)
mu0 = 2.136
mu1 = 5.543
mus       = {0:mu0,1:mu1}

##EXCITION TRANSFER RATES
forster   = Forster(Rf=raios,life=lifetimes,mu=mus)

##FLUORESCENCE RATES
fluor     = Fluor(life=lifetimes)

##ISC
isc_rates     = {0:2.4E10*1E-12,1:0}
isc = ISC(rate=isc_rates)

##NONRADIATIVE DECAY RATES
nonrad  = {0:0,1:0}
nonradiative = Nonrad(rate=nonrad)

##DISSOCIATION RATE
H = {(0,0):10E12,(0,1):10E12,(1,0):10E12,(1,1):10E12}
invrad = {0:0.1,1:0.1}
dissociation = Dissociation(AtH=H,invrad=invrad,T=300)

##FORSTER WITHOUT AVERAGED ORIENTATION FACTOR
###ForsterKappa
#forster   = ForsterKappa(Rf=raios,life=lifetimes,mu=mus)
#dipoles = np.loadtxt('dipoles.txt')
#########################################################################################

###TRIPLET EXCITONS######################################################################

##DEXTER RADII (Å)
Rds = {(0,0):10, (0,1):0, (1,0):0, (1,1):10}

##PHOSPHORESCENCE LIFETIMES (PS)
phlife = {0:5.29,1:5.29}

##SUM OF VAN DER WALLS RADII (Å)
Ls = {0:5.0,1:5.0}

##TRIPLET TRANSFER RATES
dexter = Dexter(Rd=Rds,life=phlife,L=Ls)

##PHOSPHORESCENCE RATE
phosph = Phosph(life=phlife)
#########################################################################################

###CHARGES###############################################################################

##relative permitivity
relative_eps       = 3.5   

##ATTEMPT-TO-HOP FREQUENCY (?)
H = {(0,0):10E12,(0,1):10E12,(1,0):10E12,(1,1):10E12}

##INVERSE LOCALIZATION RADIUS ()
invrad = {0:1.5,1:1.5}

##MILLER-ABRAHAMS RATE
miller = MillerAbrahams(AtH=H,invrad=invrad,T=300)


###PROCESSES#############################################################################

processes = {'singlet':[forster], 'triplet':[dexter], 'electron':[miller],'hole':[miller]}
monomolecular = {'singlet':[fluor],'triplet':[phosph],'electron':[],'hole':[]}
#########################################################################################

###MORPHOLOGY############################################################################

##Morphology functions
lattice_filename   = "lattice.txt" # file name of the system's morphology
X,Y,Z,Mats = morphology.read_lattice(lattice_filename)

##ENERGIES
#Gaussian distribuitions
s1s = {0:(3.7,0.0), 1:(2.85,0.0)} #(Peak emission energy (eV), Desvio padrao emissao (eV)
t1s = {0:(6.1,0.0), 1:(5.25,0.0)} # triplet energy, disperison (eV)

#If you have your own distribuition already
#s1s = {0:'s1_mat0.txt', 1:'s1_mat1.txt'} 
#t1s = {0:'t1_mat0.txt', 1:'t1_mat1.txt'}

#ener_function     = morphology.s1_t1_distr
ener_function      = morphology.homo_lumo
parameters_enefunc = [s1s, t1s, Mats]
s1, t1 = ener_function(parameters_enefunc)    
ene_dic = {'s1':s1, 't1':t1, 'HOMO':t1,'LUMO':s1} #careful, if you choose dissociation, you also must give HOMO and LUMO
#########################################################################################

##GENERATE PARTICLES#####################################################################

##INITIAL NUMBER OF PARTICLES
num_ex             = 3     

#Type of particle
#gen_function       = morphology.gen_pair_elechole
gen_function        = morphology.gen_excitons  #TRIPLETS OR SINGLETS??
#gen_function       = morphology.gen_electron
#gen_function       = morphology.gen_hole

#Shape of the particle's generation
#Getting filter funcs from morphology		
shape_dic = {'sphere': morphology.sphere_conditional, 'plane':morphology.plane_conditional,'cone':morphology.cone_conditional,'cilinder':morphology.cilinder_conditional,'rectangle': morphology.rectangle_conditional,'free':morphology.no_conditional}
selection = morphology.filter_selection(X,Y,Z,Mats,shape_dic,mat=[None],shape="free",origin=None,argum=None) #if you dont want to mess with the formation
#selection = morphology.filter_selection(X,Y,Z,Mats,shape_dic,mat=[None],shape="rectangle",origin=None,argum=[[40,50],[0,60],[0,40]])
#selection = morphology.filter_selection(X,Y,Z,Mats,shape_dic,mat=[None],shape="sphere",origin=[30,30,30],argum=15)
#selection = morphology.filter_selection(X,Y,Z,Mats,shape_dic,mat=[None],shape="plane",origin=[10,10,10],argum=[60,60,0])

parameters_genfunc = [num_ex,selection]
#########################################################################################

##ANNIHILATION OPTIONS###################################################################

##TURN ON ANNIHILATION
anni               = True

##list of all annihi funcs that will be used
annihi_funcs_array = [morphology.anni_ele_hol,morphology.anni_sing] 
#########################################################################################


# PROCEED ONLY IF YOU KNOW WHAT YOU ARE DOING! 
def make_system():
    system = System(X,Y,Z,Mats)   
    system.set_basic_info(monomolecular,processes,identifier,animation_mode,time_limit,pause,anni,annihi_funcs_array) 
    system.set_energies(ene_dic)
    try:
        system.set_dipoles(dipoles)
    except:
        pass
    system.set_medium(relative_eps)
    excitons = gen_function(parameters_genfunc)
    system.set_particles(excitons)

    return system
   

