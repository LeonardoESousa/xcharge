
import numpy as np
import random
from kmc.particles import *

epsilon_vaccum = 8.854187e-12    #Permitivity in C/Vm
e              = -1.60217662e-19 #Electron charge    
kb             = 8.617e-5        #Boltzmann constant
hbar           = 6.582e-16       #Reduced Planck's constant



###RATES#################################################################################    

##FUNCTION FOR SETTING RADII#############################################################
def raios(num,Rf,mat,lifetime,mats):
    Raios = np.empty(num)
    Raios.fill(Rf[(mat,mat)])
    materiais = [i for i in lifetime.keys() if i != mat]
    for m in materiais:
        Raios[mats == m] =  Rf[(mat,m)]
    return Raios

def raios_dist(num,Rf,mat,lifetime,mats):
    Raios = np.array(random.choices(Rf[(mat,mat)][:,0],Rf[(mat,mat)][:,1],k=num))
    materiais = [i for i in lifetime.keys() if i != mat]
    for m in materiais:
        R2 = np.array(random.choices(Rf[(mat,m)][:,0],Rf[(mat,m)][:,1],k=num))
        Raios[mats == m] = R2[mats == m]
    return Raios


#########################################################################################

##STANDARD FORSTER TRANSFER RATE#########################################################
class Forster:
    def __init__(self,**kwargs):
        self.kind = 'jump'
        self.Rf = kwargs['Rf']
        self.lifetime = kwargs['life']
        self.mu = kwargs['mu']
        self.alpha = 1.15*0.53


    def rate(self,**kwargs):
        r      = kwargs['r']
        system = kwargs['system']
        ex     = kwargs['particle']
        mats   = system.mats    
        local  = ex.position    
        mat = mats[local]
        
        Rf = raios(len(mats),self.Rf,mat,self.lifetime,mats)
        
        x = (Rf/(self.alpha*self.mu[mat] + r))
        taxa = (1/self.lifetime[mat])*x*x*x*x*x*x
        taxa[r == 0] = 0
        return taxa


    def action(self,particle,system,local):
        particle.move(local)
#########################################################################################

##TRIPLET TO SINGLET FORSTER TRANSFER####################################################
class ForsterT:
    def __init__(self,**kwargs):
        self.kind = 'jump'
        self.Rf = kwargs['Rf']
        self.lifetime = kwargs['life']
        self.mu = kwargs['mu']
        self.alpha = 1.15*0.53

    def rate(self,**kwargs):
        r      = kwargs['r']
        system = kwargs['system']
        ex     = kwargs['particle']
        mats   = system.mats 
        local  = ex.position 
        mat = mats[local]
        
        Rf = raios(len(mats),self.Rf,mat,self.lifetime,mats)
        x = (Rf/(self.alpha*self.mu[mat] + r))
        taxa = (1/self.lifetime[mat])*x*x*x*x*x*x
        taxa[r == 0] = 0
        return taxa

    def action(self,particle,system,local):
        system.set_particles([Singlet(local)])
        particle.kill('tts',system,system.t1,'converted')
#########################################################################################

##FORSTER TRANSFER WITH ORIENTATION FACTORS ON THE FLY###################################
class ForsterKappa:
    def __init__(self,**kwargs):
        self.kind = 'jump'
        self.Rf = kwargs['Rf']
        self.lifetime = kwargs['life']
        self.mu = kwargs['mu']
        self.alpha = 1.15*0.53

    def rate(self,**kwargs):
        r      = kwargs['r']
        system = kwargs['system']
        ex     = kwargs['particle']
        mats   = system.mats    
        local  = ex.position    
        mat = mats[local]
        mus = np.copy(system.mu)

        R = np.copy(system.R) 
        dR = R - R[local,:]
        modulo = np.sqrt(np.sum(dR*dR,axis=1))[:,np.newaxis]
        dR /= modulo

        kappa = np.inner(mus[local,:],mus) -  3*(np.inner(mus[local,:],dR)*(np.sum(mus*dR,axis=1)))  
        Rf = raios(len(mats),self.Rf,mat,self.lifetime,mats)
        
        x = (Rf/(self.alpha*system.norma_mu[local] + r))
        taxa = (1/self.lifetime[mat])*(kappa*kappa)*x*x*x*x*x*x
        taxa[r == 0] = 0
        return taxa

    def action(self,particle,system,local):
        particle.move(local)
#########################################################################################

##STANDARD FORSTER TRANSFER RATE#########################################################
class ForsterRedShift:
    def __init__(self,**kwargs):
        self.kind = 'jump'
        self.Rf = kwargs['Rf']
        self.lifetime = kwargs['life']
        self.mu = kwargs['mu']
        self.T  = kwargs['T']
        self.alpha = 1.15*0.53


    def rate(self,**kwargs):
        r      = kwargs['r']
        system = kwargs['system']
        ex     = kwargs['particle']
        mats   = system.mats    
        local  = ex.position    
        mat = mats[local]
        num = len(mats)
        
        Rfs = raios_dist(num,self.Rf,mat,self.lifetime,mats)
        
        s1s   = np.copy(system.s1)
        s1s   = (s1s - s1s[local]) + abs(s1s - s1s[local]) 
        boltz = np.exp(-1*s1s/(2*kb*self.T)) 
        x = Rfs/(self.alpha*self.mu[mat] + r)
        taxa  = (1/self.lifetime[mat])*x*x*x*x*x*x*boltz
        taxa[r == 0] = 0
        return taxa


    def action(self,particle,system,local):
        particle.move(local)
#########################################################################################

##STANDARD DEXTER TRANSFER RATE##########################################################
class Dexter:
    def __init__(self,**kwargs):
        self.kind = 'jump'
        self.Rd = kwargs['Rd']
        self.lifetime = kwargs['life']
        self.L = kwargs['L']
       
    def rate(self,**kwargs):
        r      = kwargs['r']
        system = kwargs['system']
        ex     = kwargs['particle']
        mats   = system.mats  
        local  = ex.position  
        mat = mats[local]
        
        Rd = raios(len(mats),self.Rd,mat,self.lifetime,mats)
        
        taxa = (1/self.lifetime[mat])*np.exp((2*Rd/self.L[mat])*(1-r/Rd))
        taxa[r == 0] = 0
        return taxa

    def action(self,particle,system,local):
        particle.move(local)
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
        local    = particle.position 
        mats     = system.mats        
        mat      = mats[local]
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
        taxae = (1e-12)*(AtH)*np.exp(-2*in_loc_rad*r)*np.exp(-(DEe+abs(DEe))/(2*kb*self.T))
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
        mats     = system.mats        
        mat      = mats[local]
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
        taxah = (1e-12)*(AtH)*np.exp(-2*in_loc_rad*r)*np.exp(-(DEh+abs(DEh))/(2*kb*self.T))
        taxah[r == 0] = 0
        return taxah
                 
    def action(self,particle,system,local):
        e = Electron(particle.position)
        h = Hole(local) 
        h.identity = -1*e.identity   
        system.set_particles([e,h])
        particle.kill(self.kind,system,system.s1,'converted')
#########################################################################################


def corrected_energies(system,s,r):
    potential = np.copy(system.electrostatic())
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
        particle  = kwargs['particle']
        mats      = system.mats        
        mat       = mats[particle.position]
        
        AtH        = raios(len(r),self.AtH,mat,self.inv,mats)
        in_loc_rad = self.inv[mat]

        DE = corrected_energies(system,particle,r) 
        taxa = (1e-12)*(AtH)*np.exp(
                               -(in_loc_rad*r+in_loc_rad*r))*np.exp(-DE/((kb*self.T+kb*self.T)))                           	               
        taxa[r == 0] = 0
        return taxa
 
    def action(self,particle,system,local):
        particle.move(local)

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


