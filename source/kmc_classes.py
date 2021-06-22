import numpy as np
import random

epsilon_vaccum = (1)*8.85*10**(-22) #Permitivity in C/VAngstrom
e              = -1.60217662*(10**(-19)) #Electron charge    
kb             = 8.617*(10**(-5))    # Boltzmann constant
hbar           = 6.582*(10**(-16)) #Reduced Planck's constant
        


class System:
    def __init__(self,X,Y,Z,Mats):
        self.X = X
        self.Y = Y
        self.Z = Z
        self.R = np.hstack((X[:,np.newaxis], Y[:,np.newaxis], Z[:,np.newaxis]))
        self.mats = Mats
        self.dead = []
        self.time = 0
        self.potential_time = -1
    
    def set_particles(self,Ss):
        self.particles = Ss            
    
    def set_dipoles(self,mus):
        self.mu = mus
        self.norma_mu = np.sqrt(np.sum(mus**2,axis=1))
        self.mu /= self.norma_mu[:,np.newaxis]

    def add_particle(self,s):
        self.particles.append(s)
    
    def count_particles(self):
        return len(self.particles)
    
    def set_medium(self,eps_rel):
        self.epsilon = eps_rel*epsilon_vaccum
        
    def get_num(self):
        return len(self.X)
           
    def set_energies(self,energy_dic):
    	
        s1 = energy_dic.get("s1")
        t1 = energy_dic.get("t1")
        HOMO = energy_dic.get("HOMO")
        LUMO = energy_dic.get("LUMO")      	
        
        
        self.s1   = s1 
        self.t1   = t1
        self.HOMO = HOMO
        self.LUMO = LUMO
        
    def remove(self,particle):
        self.particles.remove(particle)
        self.dead.append(particle) 

    def electrostatic(self):
        if self.time > self.potential_time:
            potential = np.zeros(len(self.X))
            for s in self.particles:
                if s.charge != 0:
                    dx = np.nan_to_num(self.X - self.X[s.position])
                    dy = np.nan_to_num(self.Y - self.Y[s.position])
                    dz = np.nan_to_num(self.Z - self.Z[s.position])
                    r  = np.sqrt(dx**2+dy**2+dz**2)
                    r[r == 0] = np.inf
                    potential += s.charge*abs(e)/(4*np.pi*self.epsilon*r)
            self.potential = potential
            self.potential_time = self.time
        else:
            pass   
        #print(self.potential)
        #input()
        return self.potential
            
   
class Particles:
    def __init__(self,species,initial):
        self.species = species
        self.initial = initial
        self.position = initial
        self.status = 'alive'
        self.identity = random.uniform(0,5)
        self.report = ''
    
    def move(self,local):
        self.position = local

    
    def make_text(self,system,energies,causamortis):
        time = system.time   #clock()
        X,Y,Z = system.X, system.Y, system.Z   #get_XYZ()
        Mats  = system.mats    #get_mats()
        x0,y0,z0 = X[self.initial],Y[self.initial],Z[self.initial]
        x, y, z  = X[self.position],Y[self.position],Z[self.position]
        dx = np.nan_to_num(x-x0)
        dy = np.nan_to_num(y-y0)
        dz = np.nan_to_num(z-z0)
        mat = Mats[self.position]
        energy = energies[self.position]
        texto = '{0:5.6e}    {1: 6.2f} {2: 6.2f} {3: 6.2f} {4:4} {5: 4.4f} {6:9} {7: 6.2f} {8: 6.2f} {9: 6.2f} {10:4}'.format(
                       time,dx,dy,dz,self.species,energy,mat,x,y,z,causamortis)
        self.report += texto+'\n'
    
    def kill(self,causamortis,system,energies):
        self.status = 'dead'
        self.make_text(system,energies,causamortis)
        system.remove(self)
    
    def convert(self,system,energies,causamortis,newkind):
        self.make_text(system,energies,causamortis)
        self.species = newkind

    def write(self):
        return self.report
        
 
class Electron(Particles):
    def __init__(self,initial):
        Particles.__init__(self,'electron',initial) 
        self.charge = -1

class Exciton(Particles):
    def __init__(self,kind,initial):
        Particles.__init__(self,kind,initial) 
        self.charge = 0

class Hole(Particles):
    def __init__(self,initial):
        Particles.__init__(self,'hole',initial) 
        self.charge = 1

###TAXAS#################################################################################    
def raios(num,Rf,mat,lifetime,mats):
    Raios = np.zeros(num) + Rf[(mat,mat)]
    materiais = [i for i in lifetime.keys() if i != mat]
    for m in materiais:
        R2 = np.zeros(num) + Rf[(mat,m)]
        Raios[mats == m] = R2[mats == m]
    return Raios


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
        num = len(mats)
        
        Rf = raios(num,self.Rf,mat,self.lifetime,mats)
        
        lifetime = self.lifetime[mat]
        mu       = self.mu[mat]
        
        taxa = (1/lifetime)*((Rf/(self.alpha*mu + r))**6)
        taxa = np.nan_to_num(taxa)
        return taxa


    def action(self,particle,system,local):
        particle.move(local)

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
        num = len(mats)
        
        Rf = raios(num,self.Rf,mat,self.lifetime,mats)
        lifetime = self.lifetime[mat]
        mu       = self.mu[mat]
        
        taxa = (1/lifetime)*((Rf/(self.alpha*mu + r))**6)
        taxa = np.nan_to_num(taxa)
        return taxa

    def action(self,particle,system,local):
        particle.move(local)
        energies = system.t1 
        particle.convert(system,energies,self.kind,'singlet')

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
        num = len(mats)
        mus = np.copy(system.mu)

        R = np.copy(system.R) 
        dR = R - R[local,:]
        modulo = np.sqrt(np.sum(dR**2,axis=1))[:,np.newaxis]
        dR /= modulo

        kappa = np.inner(mus[local,:],mus) -  3*(np.inner(mus[local,:],dR)*(np.sum(mus*dR,axis=1)))  
        #print(np.nanmean(kappa**2),max(kappa**2),min(kappa**2))
        Rf = raios(num,self.Rf,mat,self.lifetime,mats)
        
        lifetime = self.lifetime[mat]
        mu       = system.norma_mu[local]

        taxa = (1/lifetime)*(kappa**2)*((Rf/(self.alpha*mu + r))**6)
        #taxa = np.nan_to_num(taxa)
        return taxa

    def action(self,particle,system,local):
        particle.move(local)

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
        num = len(mats)
        
        Rd = raios(num,self.Rd,mat,self.lifetime,mats)
        
        lifetime = self.lifetime[mat]
        L        = self.L[mat]
        taxa = (1/lifetime)*np.exp((2*Rd/L)*(1-r/Rd))
        return taxa

    def action(self,particle,system,local):
        particle.move(local)


class Fluor:
    def __init__(self,**kwargs):
        self.kind = 'fluor'
        self.lifetime = kwargs['life']

    def rate(self,**kwargs):
        material = kwargs['material']
        lifetime = self.lifetime.get(material)
        taxa = 1/lifetime
        return taxa
     
    def action(self,particle,system,local):
        particle.kill(self.kind,system,system.s1) 

class Phosph:
    def __init__(self,**kwargs):
        self.kind = 'phosph'
        self.lifetime = kwargs['life']

    def rate(self,**kwargs):
        material = kwargs['material']
        lifetime = self.lifetime.get(material)
        taxa = 1/lifetime
        return taxa
     
    def action(self,particle,system,local):
        particle.kill(self.kind,system,system.t1)
        
class Nonrad:
    def __init__(self,**kwargs):
        self.kind = 'nonrad'
        self.taxa = kwargs['rate']

    def rate(self,**kwargs):
        material = kwargs['material']
        taxa = self.taxa.get(material)
        return taxa
     
    def action(self,particle,system,local):
        particle.kill(self.kind,system,system.s1)        

class ISC:
    def __init__(self,**kwargs):
        self.kind = 'isc'
        self.taxa = kwargs['rate']

    def rate(self,**kwargs):
        material = kwargs['material']
        taxa = self.taxa.get(material)
        return taxa
     
    def action(self,particle,system,local):
        if particle.species == 'singlet':
            energies = system.s1 
            newkind  = 'triplet'
        elif particle.species == 'triplet': 
            energies = system.t1  
            newkind  = 'singlet'   
        particle.convert(system,energies,self.kind,newkind)


class ForsterDye:
    def __init__(self,**kwargs):
        self.kind = 'fluor'
        self.lcc = kwargs['lcc'] #C-C distances in Angstrom
        self.t   = kwargs['t']   #Hopping integral in graphene in eV
        self.mu  = kwargs['mu']
        self.eps = kwargs['eps']
        self.switch = kwargs['switch']
    
    def rate(self,**kwargs):
        r      = kwargs['r']
        system = kwargs['system']
        ex     = kwargs['particle']
        mats   = system.mats  
        local  = ex.position  
        hw     =  system.s1[local] 
        mat = mats[local]
        num = len(mats)
        
        lcc  = self.lcc #C-C distances in Angstrom
        t    = self.t   #Hopping integral in graphene in eV
        mu   = self.mu[mat]
        eps  = self.eps
        
        switch = raios(num,self.switch,mat,self.mu,mats)
        
        epsilon = system.epsilon #Permitivity in C/VAngstrom    
        vf = t*(3/2)*lcc
        q  = np.linspace(0,hw/vf,1000)
        q  = q[:-1]
        funcao =  np.array([np.trapz(np.exp(-2*q*R)*(q**3)/np.sqrt(hw**2-(q**2)*(vf**2))) for R in r])  
        taxa = switch*10**(-12)*((mu*0.53)**2)*(1/48)*(e**2)/((2*np.pi*hbar)*(epsilon**2))*funcao
        taxa = np.nan_to_num(taxa)
        return taxa

    def action(self,particle,system,local):
        particle.move(local)        
        
        

class MillerAbrahams:
    def __init__(self,**kwargs):
        self.kind = 'miller-abrahams'
        self.AtH = kwargs['AtH']
        self.inv      = kwargs['invrad']
        self.T        = kwargs['T']

    def rate(self,**kwargs):
        system = kwargs['system']
        r      = kwargs['r']
        particle = kwargs['particle']
        local = particle.position 
        mats = system.mats        
        mat = mats[local]
        num = len(mats)
        
        charge = particle.charge 

        AtH      = raios(num,self.AtH,mat,self.inv,mats)
        in_loc_rad  = self.inv[mat]
        
        np.set_printoptions(suppress=True,formatter={'float': '{: 6.3f}'.format})        
        potential = np.copy(system.electrostatic())
        
        #case where two electrons or holes overlap
        duplicate  = [ x.charge for x in system.particles if x.position == local]
        for charge in duplicate:
            potential -= charge*abs(e)/(4*np.pi*system.epsilon*r)
         
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
                engs[m] = +np.inf 
                 
            engs += 1*potential
            DE = (engs - engs[local]) + abs(engs - engs[local])
            DE *= -1
            
        

        #taxa = (10**-12)*(2*np.pi/hbar)*(H**2)*(1/np.sqrt(4*np.pi*reorg*kb*self.T))*np.exp(
        #                       -2*in_loc_rad*r)*np.exp(-((DE + reorg)**2)/(4*reorg*kb*self.T))
        #
        taxa = (10E-12)*(AtH)*np.exp(
                               -2*in_loc_rad*r)*np.exp(charge*DE/(2*kb*self.T))
                               
                               
        '''                 
        # experimental                       
        if particle.species == 'hole':
            for i in indices_e:
            	taxa[i] = 1
        
        if particle.species == 'electron':
            for i in indices_h:
            	taxa[i] = 1
        '''   	
            	
            	               
        taxa[r == 0] = 0

        taxa = np.nan_to_num(taxa)
        return taxa

    def label(self):
        return self.kind
     
    def action(self,particle,system,local):
        
        indices_e  = [ x.position for x in system.particles if x.charge == -1 ]    	
        if particle.species == 'hole' and local in indices_e:
            pass
        else:
            particle.move(local)
        

class Dissociation:
    def __init__(self,**kwargs):
        self.kind = 'dissociation'
        self.AtH = kwargs['AtH']
        self.inv      = kwargs['invrad']
        self.T        = kwargs['T']
        self.Map      = {}

    def rate(self,**kwargs):
        system = kwargs['system']
        r      = kwargs['r']
        particle = kwargs['particle']
        local = particle.position 
        mats = system.mats        
        mat = mats[local]
        num = len(mats)

        lumos = np.copy(system.LUMO)
        homos = np.copy(system.HOMO)
        s1s   = np.copy(system.s1) 

        AtH      = raios(num,self.AtH,mat,self.inv,mats)
        in_loc_rad  = self.inv[mat]
        
        np.set_printoptions(suppress=True,formatter={'float': '{: 6.3f}'.format})        
        r[r == np.inf] = 0 
        
        DEe = lumos - (homos[local] + s1s[local])
        DEh = (lumos[local] - s1s[local]) - homos  
        

        taxae = (10E-12)*(AtH)*np.exp(-2*in_loc_rad*r)*np.exp((DEe+abs(DEe))/(2*kb*self.T))
        taxah = (10E-12)*(AtH)*np.exp(-2*in_loc_rad*r)*np.exp((DEh+abs(DEh))/(2*kb*self.T))
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


