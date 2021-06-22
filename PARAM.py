import sys

bib_path = "XXX"
sys.path.insert(1, bib_path)
import morphology
from kmc_classes import *



# KMC PARAMETERS 
#BASIC PARAMETERS
identifier         = 'New'
time_limit         = np.inf
animation_mode     = True
anni               = True
pause              = False # if you want to annimation stops in the first frame (debug purposes)
rounds             = 10 #number of rounds
n_proc             = 1
num_ex             = 300 #number of excitons
relative_eps       = 3.5 #relative permitivity
lattice_filename   = "lattice.txt"

 
###SINGLET RATES
r00 = 25   #Forster radius material 0 --> material 0 (Angstrom)    
r01 = 25   #material 0 --> material 1      
r10 = 25       
r11 = 25     

f0 = 2900000 #lifetime of material 0
f1 = 2900 #lifetime of material 1

#dipoles (a.u.)
mu0 = 2.136
mu1 = 5.543

#EXCITION SINGLET RATES
raios     = {(0,0):r00, (0,1):r01, (1,0):r10, (1,1):r11}
lifetimes = {0:f0,1:f1}
mus       = {0:mu0,1:mu1}

forster   = Forster(Rf=raios,life=lifetimes,mu=mus)
fluor     = Fluor(life=lifetimes)

#ISC
isc_rates     = {0:2.4E10*1E-12,1:0}
isc = ISC(rate=isc_rates)

#NONRADIATIVE RATES
nonrad  = {0:0,1:0}
nonradiative = Nonrad(rate=nonrad)

#ENERGIES
s1s = {0:(3.7,0.0), 1:(2.85,0.0)} #(Peak emission energy (eV), Desvio padrao emissao (eV)
t1s = {0:(6.1,0.0), 1:(5.25,0.0)} # triplet energy, disperison (eV)

#TRIPLET RATES
Rds = {(0,0):10, (0,1):0, (1,0):0, (1,1):10}
phlife = {0:5.29,1:5.29}
Ls = {0:5.0,1:5.0}

dexter = Dexter(Rd=Rds,life=phlife,L=Ls)
phosph = Phosph(life=phlife)

###MillerAbrahams RATES
H = {(0,0):10E12,(0,1):10E12,(1,0):10E12,(1,1):10E12}
invrad = {0:1.5,1:1.5}
miller = MillerAbrahams(AtH=H,invrad=invrad,T=300)

###Dissociation
H = {(0,0):10E12,(0,1):10E12,(1,0):10E12,(1,1):10E12}
invrad = {0:0.1,1:0.1}
dissociation = Dissociation(AtH=H,invrad=invrad,T=300)


###ForsterKappa
forster   = Forster(Rf=raios,life=lifetimes,mu=mus)
#forster   = ForsterKappa(Rf=raios,life=lifetimes,mu=mus)

#PROCESSES
processes = {'singlet':[forster], 'triplet':[dexter], 'electron':[miller],'hole':[miller]}
monomolecular = {'singlet':[fluor],'triplet':[phosph],'electron':[],'hole':[]}

#Morphology functions
X,Y,Z,Mats = morphology.read_lattice(lattice_filename)

#Type of particle
#gen_function       = morphology.gen_pair_elechole
gen_function       = morphology.gen_excitons
#gen_function       = morphology.gen_electron
#gen_function       = morphology.gen_hole

#Shape of the particle's generation
#Getting filter funcs from morphology		
shape_dic = {'sphere': morphology.sphere_conditional, 'plane':morphology.plane_conditional,'cone':morphology.cone_conditional,'cilinder':morphology.cilinder_conditional,'rectangle': morphology.rectangle_conditional,'free':morphology.no_conditional}



selection = range(len(X)) #if you dont want to mess with the formation
#selection = morphology.filter_selection(X,Y,Z,Mats,shape_dic,mat=[None],shape="free",origin=None,argum=None)
#selection = morphology.filter_selection(X,Y,Z,Mats,shape_dic,mat=[None],shape="rectangle",origin=None,argum=[[40,50],[0,60],[0,40]])
#selection = morphology.filter_selection(X,Y,Z,Mats,shape_dic,mat=[None],shape="sphere",origin=[30,30,30],argum=15)
#selection = morphology.filter_selection(X,Y,Z,Mats,shape_dic,mat=[None],shape="plane",origin=[10,10,10],argum=[60,60,0])

parameters_genfunc = [num_ex,selection]


ener_function      = morphology.homo_lumo
parameters_enefunc = [s1s, t1s, Mats]

annihi_funcs_array = [morphology.anni_ele_hol,morphology.anni_sing]#list of all annihi funcs that will be used


#### GENERATE THE SYSTEM

s1, t1 = ener_function(parameters_enefunc)
dipoles = np.loadtxt('dipoles.txt')

ene_dic = {'s1':s1, 't1':t1, 'HOMO':t1,'LUMO':s1} #careful, if you choosed dissociation, you also must give HOMO and LUMO
#ene_dic = {'s1':s1, 't1':t1}

def make_system():

    system = System(X,Y,Z,Mats)    
    system.set_energies(ene_dic)
    system.set_dipoles(dipoles)
    system.set_medium(relative_eps)
    excitons = gen_function(parameters_genfunc)
    system.set_particles(excitons)
    return system
    


