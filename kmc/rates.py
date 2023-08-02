
import numpy as np
import random
from kmc.particles import *
import kmc.utils

epsilon_vaccum = 8.854187e-12    #Permitivity in C/Vm
e              = -1.60217662e-19 #Electron charge    
kb             = 8.617e-5        #Boltzmann constant
hbar           = 6.582e-16       #Reduced Planck's constant



###RATES#################################################################################    

##FUNCTION FOR SETTING RADII#############################################################

def raios(num,Rf, mat, lifetime, mats):
  # Initialize the Raios array with the value of Rf[(mat,mat)]
  Raios = np.empty(num)
  Raios.fill(Rf[(mat,mat)])

  # Use NumPy's where function to set the values of Raios for the other materials
  for m in lifetime.keys():
    if m != mat:
      Raios = np.where(mats == m, Rf[(mat,m)], Raios)

  return Raios


def raios_dist(num,Rf,mat,lifetime,mats):
    Raios = np.array(random.choices(Rf[(mat,mat)][:,0],Rf[(mat,mat)][:,1],k=num))
    materiais = [i for i in lifetime.keys() if i != mat]
    for m in materiais:
        R2 = np.array(random.choices(Rf[(mat,m)][:,0],Rf[(mat,m)][:,1],k=num))
        Raios[mats == m] = R2[mats == m]
    return Raios


#function to convert dictionary with (i,j) keys to ixj array
def dict_to_array(d):
    keys = d.keys()
    num_keys = len(set(key[0] for key in keys))
    radius = np.empty((num_keys,num_keys))
    for key in keys:
        radius[key[0],key[1]] = d[key]
    return radius

#########################################################################################

##STANDARD FORSTER TRANSFER RATE#########################################################
class Forster:
    def __init__(self,**kwargs):
        self.kind = 'jump'
        self.Rf = dict_to_array(kwargs['Rf'])
        self.lifetime = kwargs['life']
        self.mu = kwargs['mu']
        self.alpha = 1.15*0.53


    def rate(self,**kwargs):
        r      = kwargs['r']
        system = kwargs['system']
        ex     = kwargs['particle']
        mats   = kwargs['mats']
        local  = ex.position  
        mat    = kwargs['matlocal']
        num = len(mats)
        taxa = kmc.utils.forster(self.Rf[mat,:],mats,num,self.alpha*self.mu[mat], r,1/self.lifetime[mat])
        return taxa


    def action(self,particle,system,local):
        particle.move(local,system)
#########################################################################################

##TRIPLET TO SINGLET FORSTER TRANSFER####################################################
class ForsterT:
    def __init__(self,**kwargs):
        self.kind = 'jump'
        self.Rf = dict_to_array(kwargs['Rf'])
        self.lifetime = kwargs['life']
        self.mu = kwargs['mu']
        self.alpha = 1.15*0.53

    def rate(self,**kwargs):
        r      = kwargs['r']
        system = kwargs['system']
        ex     = kwargs['particle']
        mats   = kwargs['mats']
        local  = ex.position  
        mat    = kwargs['matlocal']
        num = len(mats)
        taxa = kmc.utils.forster(self.Rf[mat,:],mats,num,self.alpha*self.mu[mat], r,1/self.lifetime[mat])
        return taxa

    def action(self,particle,system,local):
        particle.move(local,system)
        particle.kill('tts',system,system.t1,'converted')
        system.set_particles([Singlet(local)])
#########################################################################################

##FORSTER ANNIHILATION RADIUS#########################################################
class Forster_Annirad:
    def __init__(self,**kwargs):
        self.kind = 'jump'
        self.Rf = dict_to_array(kwargs['Rf'])
        self.lifetime = kwargs['life']
        self.mu = kwargs['mu']
        self.alpha = 1.15*0.53
        self.anni_rad = kwargs['anni_rad']

    def rate(self,**kwargs):
        r      = kwargs['r']
        system = kwargs['system']
        ex     = kwargs['particle']
        mats   = kwargs['mats']
        cut    = kwargs['cut']
        local  = ex.position  
        mat    = kwargs['matlocal']
        num = len(mats)
        relevant_particles = [p for p in system.particles if p.identity != ex.identity and p.position in cut]
        ss = [(np.where(cut == p.position)[0][0],self.anni_rad[(mat,system.mats[p.position])][p.species]) for p in relevant_particles]
        replace_pos   = np.array([ele[0] for ele in ss],dtype=np.int32)
        replace_raios = np.array([ele[1] for ele in ss],dtype=np.double)
        mum = len(replace_pos)
        taxa = kmc.utils.forster_anni(self.Rf[mat,:],mats,num,self.alpha*self.mu[mat], r,1/self.lifetime[mat], replace_pos, replace_raios, mum)
        return taxa


    def action(self,particle,system,local):
        particle.move(local,system)       
#########################################################################################

##STANDARD DEXTER TRANSFER RATE##########################################################
class Dexter:
    def __init__(self,**kwargs):
        self.kind = 'jump'
        self.Rd = dict_to_array(kwargs['Rd'])
        self.lifetime = kwargs['life']
        self.L = kwargs['L']
       
    def rate(self,**kwargs):
        r      = kwargs['r']
        system = kwargs['system']
        ex     = kwargs['particle']
        mats   = kwargs['mats']
        local  = ex.position  
        mat    = kwargs['matlocal']
        num = len(mats)
        taxa = kmc.utils.dexter(self.Rd[mat,:],1/self.L[mat],1/self.lifetime[mat], mats,num, r)
        return taxa

    def action(self,particle,system,local):
        particle.move(local,system)
#########################################################################################       

##EXCITON DISSOCIATION BY ELECTRON HOP RATE##############################################
class Dissociation_electron:
    def __init__(self,**kwargs):
        self.kind = 'dissociation_e'
        self.AtH = kwargs['AtH']
        self.inv      = kwargs['invrad']
        self.T        = kwargs['T']

    def rate(self,**kwargs):
        system   = kwargs['system']
        r        = kwargs['r']
        particle = kwargs['particle']
        mats   = kwargs['mats']
        mat    = kwargs['matlocal']
        num      = len(mats)

        lumos = np.copy(system.lumo)
        homos = np.copy(system.homo)
        if particle.species   == 'singlet':
            s1s   = np.copy(system.s1) 
        elif particle.species == 'triplet':
            s1s   = np.copy(system.t1)

        AtH        = raios(num,self.AtH,mat,self.inv,mats)
        in_loc_rad = self.inv[mat]
        
        DEe = lumos - (homos[local] + s1s[local])
        taxae = (1e-12)*(AtH)*np.exp(-2*in_loc_rad*r -(DEe+abs(DEe))/(2*kb*self.T))
        taxae[r == 0] = 0
        return taxae
            
     
    def action(self,particle,system,local):
        e = Electron(local)
        h = Hole(particle.position)
        h.identity = -1*e.identity   
        system.set_particles([e,h])
        particle.kill(self.kind,system,system.s1,'converted')
#########################################################################################

##EXCITON DISSOCIATION BY HOLE HOP RATE##################################################
class Dissociation_hole:
    def __init__(self,**kwargs):
        self.kind = 'dissociation_h'
        self.AtH = kwargs['AtH']
        self.inv      = kwargs['invrad']
        self.T        = kwargs['T']

    def rate(self,**kwargs):
        system   = kwargs['system']
        r        = kwargs['r']
        particle = kwargs['particle']
        local    = particle.position 
        mats   = kwargs['mats']
        mat    = kwargs['matlocal']
        num      = len(mats)
        
        lumos = np.copy(system.lumo)
        homos = np.copy(system.homo)
        if particle.species   == 'singlet':
            s1s   = np.copy(system.s1) 
        elif particle.species == 'triplet':
            s1s   = np.copy(system.t1)

        AtH        = raios(num,self.AtH,mat,self.inv,mats)
        in_loc_rad = self.inv[mat]
        
        DEh = (lumos[local] - s1s[local]) - homos  
        taxah = (1e-12)*(AtH)*np.exp(-2*in_loc_rad*r -(DEh+abs(DEh))/(2*kb*self.T))
        taxah[r == 0] = 0
        return taxah
                 
    def action(self,particle,system,local):
        e = Electron(particle.position)
        h = Hole(local) 
        h.identity = -1*e.identity   
        system.set_particles([e,h])
        particle.kill(self.kind,system,system.s1,'converted')
#########################################################################################


def corrected_energies(system,s,r,dx,dy,dz):
    potential = np.copy(system.electrostatic())
    potential += -1*(system.field[0]*dx + system.field[1]*dy + system.field[2]*dz)
    r[r == 0]  = np.inf 
    potential -= s.charge*abs(e)/(4*np.pi*system.epsilon*r)
    indices_e  = np.array([x.position for x in system.particles if x.charge == -1 and x.position != s.position]).astype(int)
    indices_h  = np.array([x.position for x in system.particles if x.charge ==  1 and x.position != s.position]).astype(int)

    if s.species == 'electron':
        engs  = np.copy(system.lumo)
        engs[indices_h] = system.homo[indices_h]
        engs[indices_e] = np.inf                
        engs += -1*potential
        DE = (engs - engs[s.position]) + abs(engs - engs[s.position])
    elif s.species == 'hole':
        engs  = np.copy(system.homo)
        engs[indices_e] = system.lumo[indices_e] 
        engs[indices_h] = -np.inf     
        engs += -1*potential
        DE = (engs[s.position] - engs) + abs(engs[s.position] - engs)
    return DE



##MILLER-ABRHAMS RATE####################################################################
class MillerAbrahams:
    def __init__(self,**kwargs):
        self.kind = 'miller-abrahams'
        self.AtH  = kwargs['AtH']
        self.inv  = kwargs['invrad']
        self.T    = kwargs['T']

    def rate(self,**kwargs):
        system    = kwargs['system']
        r         = kwargs['r']
        dx        = kwargs['dx']
        dy        = kwargs['dy']
        dz        = kwargs['dz']
        particle  = kwargs['particle']
        mats   = kwargs['mats']
        mat    = kwargs['matlocal']
        num      = len(mats)        
        
        AtH        = raios(len(r),self.AtH,mat,self.inv,mats)
        in_loc_rad = self.inv[mat]

        DE = corrected_energies(system,particle,r,dx,dy,dz)                            	               
        taxa = (1e-12)*(AtH)*np.exp(-2*in_loc_rad*r -DE/(2*kb*self.T)) 
        
        taxa[r == 0] = 0
        return taxa
 
    def action(self,particle,system,local):
        particle.move(local,system)

######################################################################################### 



#MONOMOLECULAR RATES#####################################################################

##FLUORESCENCE RATE######################################################################
class Fluor:
    def __init__(self,**kwargs):
        self.kind = 'fluor'
        self.lifetime = kwargs['life']

    def rate(self,**kwargs):
        return 1/self.lifetime[kwargs['material']]
     
    def action(self,particle,system,local):
        particle.kill(self.kind,system,system.s1,'dead') 
#########################################################################################

##PHOSPHORESCENCE RATE###################################################################
class Phosph:
    def __init__(self,**kwargs):
        self.kind = 'phosph'
        self.lifetime = kwargs['life']

    def rate(self,**kwargs):
        return 1/self.lifetime[kwargs['material']]
     
    def action(self,particle,system,local):
        particle.kill(self.kind,system,system.t1,'dead')
#########################################################################################
 
##NONRADIATIVE DECAY RATE################################################################         
class Nonrad:
    def __init__(self,**kwargs):
        self.kind = 'nonrad'
        self.taxa = kwargs['rate']

    def rate(self,**kwargs):
        return self.taxa[kwargs['material']]
     
    def action(self,particle,system,local):
        particle.kill(self.kind,system,system.s1,'dead')        
#########################################################################################

##INTERSYSTEM CROSSING RATE##############################################################
class ISC:
    def __init__(self,**kwargs):
        self.kind = 'isc'
        self.taxa = kwargs['rate']
        self.map  = {'singlet':'triplet', 'triplet':'singlet'}

    def rate(self,**kwargs):
        material = kwargs['material']
        return self.taxa[material]
     
    def action(self,particle,system,local):
        if particle.species == 'singlet':
            system.set_particles([Triplet(particle.position)])
            particle.kill(self.kind,system,system.s1,'converted')
        elif particle.species == 'triplet':    
            system.set_particles([Singlet(particle.position)])
            particle.kill('r'+self.kind,system,system.s1,'converted')
#########################################################################################


