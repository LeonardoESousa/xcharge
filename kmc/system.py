import numpy as np
import kmc.utils
epsilon_vaccum =  8.854187e-12   #Permitivity in C/Vm
e              = -1.60217662e-19 #Electron charge    
kb             = 8.617e-5        #Boltzmann constant
hbar           = 6.582e-16       #Reduced Planck's constant
        


class System:
    def __init__(self):
        self.dead = []
        self.time = 0
        self.potential_time = -1
        self.IT = 0 #number of steps
        self.annihi_radius = {}
        self.positions_by_species = {}
        self.neighbor_cache = None
                
    def set_morph(self,X,Y,Z,Mats):
        self.X = X
        self.Y = Y
        self.Z = Z
        self.Lx = max(X) - min(X)
        self.Ly = max(Y) - min(Y)
        self.Lz = max(Z) - min(Z)
        self.R = np.hstack((X[:,np.newaxis], Y[:,np.newaxis], Z[:,np.newaxis]))
        Mats = np.array(Mats,dtype=np.int32)
        self.mats = Mats
        self.uniq = np.unique(Mats)       


    def set_basic_info(self,monomolecular,processes,identifier,animation_mode,time_limit,pause,anni,distance):
        self.processes           = processes
        self.monomolecular       = monomolecular
        self.identifier          = identifier
        self.animation_mode      = animation_mode
        self.time_limit          = time_limit
        self.pause               = pause
        self.bimolec             = anni
        self.distance            = distance
    
    def set_particles(self,Ss):
        try:
            self.particles += Ss
        except:
            self.particles = Ss
        for p in Ss:
            self.positions_by_species.setdefault(p.species, set()).add(p.position)
               
    def reset_particles(self):
        self.particles = None 
        self.positions_by_species = {}
        
    def set_dipoles(self,mus):
        self.mu = mus
        self.norma_mu = np.sqrt(np.sum(mus*mus,axis=1))
        self.mu /= self.norma_mu[:,np.newaxis]
    
    def count_particles(self):
        return len(self.particles)
        
    def set_medium(self,eps_rel):
        self.eps_rel = eps_rel
        self.epsilon = eps_rel*epsilon_vaccum
        
    def get_num(self):
        return len(self.X)
           
    def set_energies(self,energy, kind):
        setattr(self, kind.lower(), energy)
        
    def remove(self,particle):
        self.particles.remove(particle)
        positions = self.positions_by_species.get(particle.species)
        if positions is not None:
            positions.discard(particle.position)
        self.dead.append(particle) 

    def set_electric_field(self, field):
        self.field = np.array(field)/(1e8)
        

    def electrostatic(self):
        if self.time > self.potential_time:
            potential = 0
            for s in self.particles:
                if s.charge != 0:
                    dx, dy, dz = self.distance(self, s.position)
                    r = kmc.utils.distances(dx,dy,dz,len(dx))*(1e-10) #np.sqrt(dx*dx + dy*dy + dz*dz)*(1e-10)
                    r[r == 0] = np.inf
                    potential += s.charge*abs(e)/(4*np.pi*self.epsilon*r)
            self.potential = potential
            self.potential_time = self.time   
        return self.potential
    
    def update_position(self, particle, old_position, new_position):
        positions = self.positions_by_species.setdefault(particle.species, set())
        positions.discard(old_position)
        positions.add(new_position)
    
    def build_neighbor_cache(self, distance_fn, cutoff):
        cache = []
        num_sites = len(self.X)
        for site in range(num_sites):
            dx, dy, dz = distance_fn(self, site)
            r = kmc.utils.distances(dx, dy, dz, len(dx))
            mask = np.where(r < cutoff)[0]
            cache.append(
                (
                    mask,
                    dx[mask],
                    dy[mask],
                    dz[mask],
                    r[mask],
                )
            )
        self.neighbor_cache = cache
            
