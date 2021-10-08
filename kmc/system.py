import numpy as np

epsilon_vaccum = 8.85e-22        #Permitivity in C/VAngstrom
e              = -1.60217662e-19 #Electron charge    
kb             = 8.617e-5        #Boltzmann constant
hbar           = 6.582e-16       #Reduced Planck's constant
        


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
   
  
    def set_basic_info(self,monomolecular,processes,identifier,animation_mode,time_limit,pause,anni,anni_funcs_array):
    	self.processes        = processes
    	self.monomolecular    = monomolecular
    	self.identifier       = identifier
    	self.animation_mode   = animation_mode
    	self.time_limit       = time_limit
    	self.pause            =   pause
    	self.anni_funcs_array = anni_funcs_array
    	self.anni             = anni
    
    def set_particles(self,Ss):
        self.particles = Ss            
    
    def set_dipoles(self,mus):
        if mus: #in case you dont have the dipoles
            self.mu = mus
            self.norma_mu = np.sqrt(np.sum(mus**2,axis=1))
            self.mu /= self.norma_mu[:,np.newaxis]

    def add_particle(self,s):
        self.particles.append(s)
    
    def count_particles(self):
        return len(self.particles)
    
    def set_medium(self,eps_rel):
    	self.eps_rel = eps_rel
    	self.epsilon = eps_rel*epsilon_vaccum
        
    def get_num(self):
        return len(self.X)
           
    def set_energies(self,energy, kind):
        if kind.lower() == 's1':
            self.s1   = energy 
        elif kind.lower() == 't1':
            self.t1   = energy
        elif kind.lower() == 'homo':
            self.HOMO = energy
        elif kind.lower() == 'lumo':
            self.LUMO = energy
        
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
        return self.potential
            
