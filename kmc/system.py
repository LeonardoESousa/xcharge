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
        
    def set_morph(self,X,Y,Z,Mats):
        self.X = X
        self.Y = Y
        self.Z = Z
        self.R = np.hstack((X[:,np.newaxis], Y[:,np.newaxis], Z[:,np.newaxis]))
        self.mats = Mats
        self.uniq = np.unique(Mats)       


    def set_basic_info(self,monomolecular,processes,identifier,animation_mode,time_limit,pause,anni,anni_funcs_array):
        self.processes           = processes
        self.monomolecular       = monomolecular
        self.identifier          = identifier
        self.animation_mode      = animation_mode
        self.time_limit          = time_limit
        self.pause               = pause
        self.bimolec_funcs_array = anni_funcs_array
        self.bimolec             = anni
    
    def set_particles(self,Ss):
        
        try:
            self.particles = self.particles + Ss
        except:    
            self.particles = Ss            

    def set_dipoles(self,mus):
        self.mu = mus
        self.norma_mu = np.sqrt(np.sum(mus*mus,axis=1))
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

    def set_electric_field(self, field):
        field = field/(1e8)
        comp_x = -1*field[0]*self.X
        comp_y = -1*field[1]*self.Y
        comp_z = -1*field[2]*self.Z
        self.electric_potential = comp_x + comp_y + comp_z
        

    def electrostatic(self):
        if self.time > self.potential_time:
            potential = np.copy(self.electric_potential)
            for s in self.particles:
                if s.charge != 0:
                    dx = np.nan_to_num(self.X - self.X[s.position])
                    dy = np.nan_to_num(self.Y - self.Y[s.position])
                    dz = np.nan_to_num(self.Z - self.Z[s.position])
                    r  = np.sqrt(dx*dx+dy*dy+dz*dz)*(1e-10)
                    r[r == 0] = np.inf
                    potential += s.charge*abs(e)/(4*np.pi*self.epsilon*r)
            self.potential = potential
            self.potential_time = self.time
        else:
            pass   
        return self.potential
            
