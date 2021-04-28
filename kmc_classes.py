import numpy as np

epsilon = (3.5)*8.85*10**(-22) #Permitivity in C/VAngstrom
e       = -1.60217662*(10**(-19)) #Electron charge    
kb      = 8.617*(10**(-5))    # Boltzmann constant
hbar    = 6.582*(10**(-16)) #Reduced Planck's constant
        


class System:
    def __init__(self,X,Y,Z,Mats):
        self.X = X
        self.Y = Y
        self.Z = Z
        self.mats = Mats
        self.dead = []
        self.time = 0
        self.potential_time = -1
        
    def set_clock(self,time):
        self.time = time
    
    def clock(self):
        return self.time
    
    def set_particles(self,Ss):
        self.particles = Ss            
    
    def add_particle(self,s):
        self.particles.append(s)
    
    def count_particles(self):
        return len(self.particles)
    
    def get_particles(self):
        lista = self.particles
        return lista
    
    def get_XYZ(self):
        return self.X, self.Y, self.Z
    
    def get_num(self):
        return len(self.X)
    
    def set_orbital(self,h,l):
        self.HOMO = h
        self.LUMO = l
    
    def set_s1(self, energy):
        self.s1 = energy
    
    def set_t1(self,energy):
        self.t1 = energy
        
    def get_HOMO(self):
        return self.HOMO

    def get_LUMO(self):
        return self.LUMO
    
    def get_s1(self):
        return self.s1
    
    def get_t1(self):
        return self.t1
    
    def get_mats(self):
        return self.mats
        
    def remove(self,particle):
        self.particles.remove(particle)
        self.dead.append(particle) 

    def get_dead(self):
        return self.dead

    def electrostatic(self):
        if self.time > self.potential_time:
            potential = np.zeros(len(self.X))
            for s in self.particles:
                if s.get_charge() != 0:
                    dx = np.nan_to_num(self.X - self.X[s.location()])
                    dy = np.nan_to_num(self.Y - self.Y[s.location()])
                    dz = np.nan_to_num(self.Z - self.Z[s.location()])
                    r  = np.sqrt(dx**2+dy**2+dz**2)
                    #r[r == 0] = np.inf
                    potential += s.get_charge()*abs(e)/(4*np.pi*epsilon*r)
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
        self.report = ''
    
    def move(self,local):
        self.position = local

    def kind(self):
        return self.species
    
    def make_text(self,system,energies,causamortis):
        time = system.clock()
        X,Y,Z = system.get_XYZ()
        Mats  = system.get_mats()
        x0,y0,z0 = X[self.initial],Y[self.initial],Z[self.initial]
        x, y, z  = X[self.position],Y[self.position],Z[self.position]
        dx = np.nan_to_num(x-x0)
        dy = np.nan_to_num(y-y0)
        dz = np.nan_to_num(z-z0)
        mat = Mats[self.position]
        energy = energies[self.position]
        texto = '{0:10.3f}    {1: 6.2f} {2: 6.2f} {3: 6.2f} {4:4} {5: 4.4f} {6:9} {7: 6.2f} {8: 6.2f} {9: 6.2f} {10:4}'.format(
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
        
    def location(self):
        return self.position
 
class Electron(Particles):
    def __init__(self,initial):
        Particles.__init__(self,'electron',initial) 
        self.charge = -1

    def get_charge(self):
        return self.charge

class Exciton(Particles):
    def __init__(self,kind,initial):
        Particles.__init__(self,kind,initial) 
        self.charge = 0

    def get_charge(self):
        return self.charge

class Hole(Particles):
    def __init__(self,initial):
        Particles.__init__(self,'hole',initial) 
        self.charge = 1

    def get_charge(self):
        return self.charge

###TAXAS#################################################################################    
class Forster:
    def __init__(self,**kwargs):
        self.kind = 'jump'
        self.Rf = kwargs['Rf']
        self.lifetime = kwargs['life']
        self.mu = kwargs['mu']
       
    def rate(self,**kwargs):
        r      = kwargs['r']
        system = kwargs['system']
        ex     = kwargs['particle']
        mats   = system.get_mats()
        local  = ex.location()
        
        Rf       = np.array([self.Rf.get((mats[local],i)) for i in mats])
        lifetime = self.lifetime.get(mats[local])
        mu       = self.mu.get(mats[local])
        
        alpha = 1.15*0.53
        taxa = (1/lifetime)*((Rf/(alpha*mu + r))**6)
        taxa = np.nan_to_num(taxa)
        return taxa

    def label(self):
        return self.kind

    def action(self,particle,system,local):
        particle.move(local)

class ForsterT:
    def __init__(self,**kwargs):
        self.kind = 'jump'
        self.Rf = kwargs['Rf']
        self.lifetime = kwargs['life']
        self.mu = kwargs['mu']
       
    def rate(self,**kwargs):
        r      = kwargs['r']
        system = kwargs['system']
        ex     = kwargs['particle']
        mats   = system.get_mats()
        local  = ex.location()
        
        Rf       = np.array([self.Rf.get((mats[local],i)) for i in mats])
        lifetime = self.lifetime.get(mats[local])
        mu       = self.mu.get(mats[local])
        
        alpha = 1.15*0.53
        taxa = (1/lifetime)*((Rf/(alpha*mu + r))**6)
        taxa = np.nan_to_num(taxa)
        return taxa

    def label(self):
        return self.kind

    def action(self,particle,system,local):
        particle.move(local)
        energies = system.get_t1()
        particle.convert(system,energies,self.kind,'singlet')

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
        mats   = system.get_mats()
        local  = ex.location()
        
        Rd       = np.array([self.Rd.get((mats[local],i)) for i in mats])
        lifetime = self.lifetime.get(mats[local])
        L        = self.L.get(mats[local])
        taxa = (1/lifetime)*np.exp((2*Rd/L)*(1-r/Rd))
        return taxa

    def label(self):
        return self.kind

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

    def label(self):
        return self.kind
     
    def action(self,particle,system,local):
        particle.kill(self.kind,system,system.get_s1()) 

class Phosph:
    def __init__(self,**kwargs):
        self.kind = 'phosph'
        self.lifetime = kwargs['life']

    def rate(self,**kwargs):
        material = kwargs['material']
        lifetime = self.lifetime.get(material)
        taxa = 1/lifetime
        return taxa

    def label(self):
        return self.kind
     
    def action(self,particle,system,local):
        particle.kill(self.kind,system,system.get_t1())
        
class Nonrad:
    def __init__(self,**kwargs):
        self.kind = 'nonrad'
        self.taxa = kwargs['rate']

    def rate(self,**kwargs):
        material = kwargs['material']
        taxa = self.taxa.get(material)
        return taxa

    def label(self):
        return self.kind
     
    def action(self,particle,system,local):
        particle.kill(self.kind,system,system.get_s1())        

class ISC:
    def __init__(self,**kwargs):
        self.kind = 'isc'
        self.taxa = kwargs['rate']

    def rate(self,**kwargs):
        material = kwargs['material']
        taxa = self.taxa.get(material)
        return taxa

    def label(self):
        return self.kind
     
    def action(self,particle,system,local):
        if particle.kind() == 'singlet':
            energies = system.get_s1()
            newkind  = 'triplet'
        elif particle.kind() == 'triplet': 
            energies = system.get_t1()
            newkind  = 'singlet'   
        particle.convert(system,energies,self.kind,newkind)


class ForsterDye:
    def __init__(self,**kwargs):
        self.kind = 'fluor'
        self.lcc = kwargs['lcc'] #C-C distances in Angstrom
        self.t   = kwargs['t']   #Hopping integral in graphene in eV
        self.mu  = kwargs['mu']
        self.eps = kwargs['eps']
    
    def rate(self,**kwargs):
        r = kwargs['r']
        system = kwargs['system']
        local  = kwargs['location']
        hw   =  system.get_s1()[local]
        e    = -1.60217662*(10**(-19)) #Electron charge
        hbar = 6.582*(10**(-16)) #Reduced Planck's constant
        lcc  = self.lcc #C-C distances in Angstrom
        t    = self.t #Hopping integral in graphene in eV
        mu   = self.mu
        eps  = self.eps
        epsilon = eps*8.85*10**(-22) #Permitivity in C/VAngstrom    
        vf = t*(3/2)*lcc
        q  = np.linspace(0,hw/vf,1000)
        q  = q[:-1]
        funcao = np.exp(-2*q*r)*(q**3)/np.sqrt(hw**2-(q**2)*(vf**2))
        taxa = 10**(-12)*((mu*0.53)**2)*(1/48)*(e**2)/((2*np.pi*hbar)*(epsilon**2))*np.trapz(funcao,q)
        return taxa

    def label(self):
        return self.kind

    def action(self,particle,system,local):
        particle.move(local)        
        
class MillerAbrahams:
    def __init__(self,**kwargs):
        self.kind = 'miller-abrahams'
        self.H = kwargs['H']
        self.inv      = kwargs['invrad']
        self.T        = kwargs['T']

    def rate(self,**kwargs):
        
        system = kwargs['system']
        r      = kwargs['r']
        particle = kwargs['particle']
        local = particle.location()
        mats = system.get_mats()
        
        H      = np.array([self.H.get((mats[local],i)) for i in mats])
        in_loc_rad  = self.inv.get(mats[local])
        
                
        r[r == np.inf] = 0 
        potential = 1*system.electrostatic()
        potential -= particle.get_charge()*abs(e)/(4*np.pi*epsilon*r)
        potential = np.nan_to_num(potential,posinf=np.inf,neginf=-np.inf)  
      
        if particle.kind() == 'electron':
            engs  = system.get_LUMO()
            homos = system.get_HOMO()
            indices = np.where(potential == -np.inf)
            for m in indices[0]:
                potential[m] = 0
                engs[m] = homos[m]
            engs += -1*potential
            DE = engs - engs[local]

        elif particle.kind() == 'hole':
            engs  = system.get_HOMO()
            engs += 1*potential
            DE = (engs[local] - engs)
    
            
        

        #taxa = (10**-12)*(2*np.pi/hbar)*(H**2)*(1/np.sqrt(4*np.pi*reorg*kb*self.T))*np.exp(
        #                       -2*in_loc_rad*r)*np.exp(-((DE + reorg)**2)/(4*reorg*kb*self.T))
        #
        taxa = (10**-12)*(2*np.pi/hbar)*(H**2)*np.exp(
                               -2*in_loc_rad*r)*np.exp(particle.get_charge()*(DE +abs(DE))/(2*kb*self.T))
        taxa[r == 0] = 0

        taxa = np.nan_to_num(taxa)
        #print(particle.kind())
        #print(taxa)
        #print(DE)
        #input()
        print(particle.kind(),DE[np.where(taxa == max(taxa))[0][0]])
        return taxa

    def label(self):
        return self.kind
     
    def action(self,particle,system,local):
        particle.move(local)    
