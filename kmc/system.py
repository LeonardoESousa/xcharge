import numpy as np

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
                
    def set_morph(self,X,Y,Z,Mats):
        self.X = X
        self.Y = Y
        self.Z = Z
        self.Lx = max(X) - min(X)
        self.Ly = max(Y) - min(Y)
        self.Lz = max(Z) - min(Z)
        self.R = np.hstack((X[:,np.newaxis], Y[:,np.newaxis], Z[:,np.newaxis]))
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
               
    def reset_particles(self):
        self.particles = None 
        
    def set_dipoles(self,mus):
        self.mu = mus
        self.norma_mu = np.sqrt(np.sum(mus*mus,axis=1))
        self.mu /= self.norma_mu[:,np.newaxis]
    
    def count_particles(self):
        return len(self.particles)
        
    def append_annihi_radius(self,new_radius):
        #updating annihi dic for each instance called
        self.annihi_radius.update(new_radius) 
        key  = list(new_radius.keys())[0]
        val  = new_radius[key]
        self.annihi_radius.update({(key[1],key[0]):val}) #adding the exchange dictionary

    def set_medium(self,eps_rel):
        self.eps_rel = eps_rel
        self.epsilon = eps_rel*epsilon_vaccum
        
    def get_num(self):
        return len(self.X)
           
    def set_energies(self,energy, kind):
        setattr(self, kind.lower(), energy)
        
    def remove(self,particle):
        self.particles.remove(particle)
        self.dead.append(particle) 

    def set_electric_field(self, field):
        self.field = np.array(field)/(1e8)
        

    def electrostatic(self):
        if self.time > self.potential_time:
            potential = 0
            for s in self.particles:
                if s.charge != 0:
                    dx, dy, dz = self.distance(self, s.position)
                    r = np.sqrt(dx*dx + dy*dy + dz*dz)*(1e-10)
                    r[r == 0] = np.inf
                    potential += s.charge*abs(e)/(4*np.pi*self.epsilon*r)
            self.potential = potential
            self.potential_time = self.time   
        return self.potential
            
