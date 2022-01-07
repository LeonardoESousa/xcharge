
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
        particle.move(local)
        particle.convert(system,system.t1,self.kind,'singlet')
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

##EXCITON DISSOCIATION RATE##############################################################
class Dissociation:
    def __init__(self,**kwargs):
        self.kind = 'dissociation'
        self.AtH = kwargs['AtH']
        self.inv      = kwargs['invrad']
        self.T        = kwargs['T']
        self.Map      = {}

    def rate(self,**kwargs):
        system   = kwargs['system']
        r        = kwargs['r']
        particle = kwargs['particle']
        local    = particle.position 
        mats     = system.mats        
        mat      = mats[local]
        num      = len(mats)

        lumos = np.copy(system.LUMO)
        homos = np.copy(system.HOMO)
        if particle.species   == 'singlet':
            s1s   = np.copy(system.s1) 
        elif particle.species == 'triplet':
            s1s   = np.copy(system.t1)

        AtH        = raios(num,self.AtH,mat,self.inv,mats)
        in_loc_rad = self.inv[mat]
        
        
        DEe = lumos - (homos[local] + s1s[local])
        DEh = (lumos[local] - s1s[local]) - homos  
        

        taxae = (1e-12)*(AtH)*np.exp(-2*in_loc_rad*r)*np.exp(-1*(DEe+abs(DEe))/(2*kb*self.T))
        taxah = (1e-12)*(AtH)*np.exp(-2*in_loc_rad*r)*np.exp(-1*(DEh+abs(DEh))/(2*kb*self.T))
        taxae[r == 0] = 0
        taxah[r == 0] = 0
        taxae = np.nan_to_num(taxae)
        taxah = np.nan_to_num(taxah)
        
        TE = np.sum(taxae)
        TH = np.sum(taxah)
        if random.uniform(0,1) <= TE/(TE+TH):
            self.Map[particle.identity] = 'electron'
            return taxae
        else:
            self.Map[particle.identity] = 'hole'
            return taxah    
            

    def label(self):
        return self.kind
     
    def action(self,particle,system,local):
        if self.Map[particle.identity] == 'electron':
            e = Electron(local)
            h = Hole(particle.position)
        else:
           e = Electron(particle.position)
           h = Hole(local) 

        h.identity = -1*e.identity   
        system.add_particle(e)
        system.add_particle(h)
        del self.Map[particle.identity]
        particle.kill('dissociation',system,system.s1)
#########################################################################################

##MILLER-ABRHAMS RATE####################################################################
class MillerAbrahams:
    def __init__(self,**kwargs):
        self.kind = 'miller-abrahams'
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
        
        charge   = particle.charge 

        AtH        = raios(num,self.AtH,mat,self.inv,mats)
        in_loc_rad = self.inv[mat]
        
        potential = np.copy(system.electrostatic())
        
        #case where two electrons or holes overlap
        duplicate  = [ x.charge for x in system.particles if x.position == local]
        for charge in duplicate:
            potential -= charge*abs(e)/(4*np.pi*system.epsilon*r*1e-10)
         
        indices_e  = [ x.position for x in system.particles if x.charge == -1 and x.position != local ]
        indices_h  = [ x.position for x in system.particles if x.charge == 1  and x.position != local]

        if particle.species == 'electron':
            engs  = np.copy(system.LUMO)
            homos = np.copy(system.HOMO)
            
            for m in indices_h:
                potential[m] = 0
                engs[m] = homos[m]
                
            for m in indices_e:
                potential[m] = 0
                engs[m] = -np.inf             
                
            engs += -1*potential
            DE = (engs - engs[local]) + abs(engs - engs[local])
        elif particle.species == 'hole':
            engs  = np.copy(system.HOMO)
            lumos = np.copy(system.LUMO)
          
            for m in indices_e:
                potential[m] = 0
                engs[m] = lumos[m]
                
            for m in indices_h:
                potential[m] = 0
                engs[m] = -np.inf 
                 
            engs += -1*potential
            DE = (engs[local] - engs) + abs(engs[local] - engs)
            

        taxa = (1e-12)*(AtH)*np.exp(
                               -2*in_loc_rad*r)*np.exp(-1*DE/(2*kb*self.T))
                                   	               
        taxa[r == 0] = 0
        return taxa

    def label(self):
        return self.kind
     
    def action(self,particle,system,local):
        
        indices_e  = [ x.position for x in system.particles if x.charge == -1 ]    	
        if particle.species == 'hole' and local in indices_e:
            pass
        else:
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
        particle.kill(self.kind,system,system.s1) 
#########################################################################################

##PHOSPHORESCENCE RATE###################################################################
class Phosph:
    def __init__(self,**kwargs):
        self.kind = 'phosph'
        self.lifetime = kwargs['life']

    def rate(self,**kwargs):
        return 1/self.lifetime[kwargs['material']]
     
    def action(self,particle,system,local):
        particle.kill(self.kind,system,system.t1)
#########################################################################################
 
##NONRADIATIVE DECAY RATE################################################################         
class Nonrad:
    def __init__(self,**kwargs):
        self.kind = 'nonrad'
        self.taxa = kwargs['rate']

    def rate(self,**kwargs):
        return self.taxa[kwargs['material']]
     
    def action(self,particle,system,local):
        particle.kill(self.kind,system,system.s1)        
#########################################################################################

##INTERSYSTEM CROSSING RATE##############################################################
class ISC:
    def __init__(self,**kwargs):
        self.kind = 'isc'
        self.taxa = kwargs['rate']

    def rate(self,**kwargs):
        material = kwargs['material']
        return self.taxa[material]
     
    def action(self,particle,system,local):
        if particle.species == 'singlet':
            particle.convert(system,system.s1,self.kind,'triplet')
        elif particle.species == 'triplet':    
            particle.convert(system,system.t1,self.kind,'singlet')
#########################################################################################


