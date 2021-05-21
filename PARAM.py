from kmc_classes import *
import morphology

# KMC PARAMETERS 

#BASIC PARAMETERS
animation_mode = True
anni = True
lattice_filename = "lattice.txt"
rounds = 1 #number of rounds
num_ex = 20 #number of excitons

###SINGLET RATES
r00 = 20.8   #Forster radius material 0 --> material 0 (Angstrom)    
r01 = 22.2   #material 0 --> material 1      
r10 = 21.2        
r11 = 26.5     

f0 = 2940.0 #lifetime of material 0
f1 = 747.27 #lifetime of material 1

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
s1s = {0:(1.9,0.0), 1:(1.9,0.00)} #(Peak emission energy (eV), Desvio padrao emissao (eV)
t1s = {0:(1.2,0.0), 1:(1.2,0.00)} # triplet energy, disperison (eV)

#TRIPLET RATES
Rds = {(0,0):10, (0,1):0, (1,0):0, (1,1):10}
phlife = {0:5.291E11,1:np.inf}
Ls = {0:5.0,1:5.0}

dexter = Dexter(Rd=Rds,life=phlife,L=Ls)
phosph = Phosph(life=phlife)

###MillerAbrahams RATES
H = {(0,0):0.1,(0,1):0.1,(1,0):0.1,(1,1):0.1}
invrad = {0:10.5,1:10.5}
miller = MillerAbrahams(H=H,invrad=invrad,T=300)

###Dissociation
H = {(0,0):10000.0,(0,1):10000.0,(1,0):100.0,(1,1):100.0}
invrad = {0:0,1:0}
dissociation = Dissociation(H=H,invrad=invrad,T=300)


#PROCESSES

processes = {'singlet':[forster,dissociation], 'triplet':[dexter], 'electron':[miller],'hole':[miller]}
monomolecular = {'singlet':[],'triplet':[phosph],'electron':[],'hole':[]}
#fluor,isc,nonradiative

#Morphology functions
X,Y,Z,Mats = morphology.read_lattice(lattice_filename)

#gen_function       = morphology.gen_pair_elechole
gen_function       = morphology.gen_excitons
parameters_genfunc = [num_ex,len(X)]

ener_function      = morphology.homo_lumo
parameters_enefunc = [s1s, t1s, Mats]

#annihi_funcs_array = [morphology.anni_ele_hol] 
annihi_funcs_array = [morphology.anni_sing]
#annihi_funcs_array = [morphology.anni_ele_hol,morphology.anni_sing]#list of all annihi funcs that will be used






