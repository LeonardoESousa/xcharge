import numpy as np
import random


DIFFUSION_LOG_COLUMNS = (
    'AveragedDiffusionTime',
    'D_est_xx', 'D_est_yy', 'D_est_zz', 'D_est_xy', 'D_est_xz', 'D_est_yz',
    'MSD_est_x', 'MSD_est_y', 'MSD_est_z',
    'Cov_est_xy', 'Cov_est_xz', 'Cov_est_yz',
    'Drift_est_x', 'Drift_est_y', 'Drift_est_z',
    'MeanDelta_est_x', 'MeanDelta_est_y', 'MeanDelta_est_z',
    'Var_est_x', 'Var_est_y', 'Var_est_z',
)


def diffusion_log_values(particle, system):
    """Lifetime-averaged transport and integrated unwrapped covariance.

    Particles never averaged report NaN, including retained FP products.
    """
    duration = particle.averaged_diffusion_time
    if duration <= 0:
        return [float('nan')] * len(DIFFUSION_LOG_COLUMNS)
    integral = particle.averaged_diffusion_integral
    tensor = integral / duration
    covariance = 2.0 * integral
    mean = particle.averaged_drift_integral
    second_moment = covariance + np.outer(mean, mean)
    indices = ((0,0), (1,1), (2,2), (0,1), (0,2), (1,2))
    return ([duration] + [tensor[i,j] for i,j in indices]
            + list(np.diag(second_moment))
            + [covariance[i,j] for i,j in indices[3:]]
            + list(mean/duration) + list(mean) + list(np.diag(covariance)))

class Particles:
    def __init__(self,species,initial):
        self.species  = species
        self.initial  = initial
        self.position = initial
        self.status   = 'alive'
        self.identity = random.uniform(0,5)
        self._report_lines = []
        self.texto    = ''
        self.process  = None
        self.destination  = None
        self.Dx = 0
        self.Dy = 0
        self.Dz = 0
        self.averaged_diffusion_time = 0.0
        self.averaged_diffusion_integral = np.zeros((3, 3))
        self.averaged_drift_integral = np.zeros(3)
        
    def move(self,local,system,displacement=None):
        old_position = self.position
        Dx, Dy, Dz = (system.edge_displacement(self.position, local)
                      if displacement is None else displacement)
        self.Dx += Dx
        self.Dy += Dy
        self.Dz += Dz
        self.position = local
        system.update_position(self, old_position, local)
        

    def make_text(self,system,energies,causamortis):
        X,Y,Z = system.X, system.Y, system.Z  
        Mats  = system.mats 
        x, y, z  = X[self.position],Y[self.position],Z[self.position]
        mat = Mats[self.position]
        energy = energies[self.position]
        self.texto = f'TEMPO,{self.Dx:.0f},{self.Dy:.0f},{self.Dz:.0f},{self.species},{energy:.2f},{mat:.0f},{x:.0f},{y:.0f},{z:.0f},{causamortis},{self.status}'
        if system.diffusion_averager is not None:
            self.texto += ''.join(f',{value:.12e}' for value in diffusion_log_values(self, system))
        
    def stamp_time(self,system):
        if self.texto != '':
            texto = self.texto
            texto = texto.replace('TEMPO',f'{system.time:.12e}')
            self._report_lines.append(texto + '\n')
            self.texto = ''
        
    def kill(self,causamortis,system,energies,result):
        #print(f'killing a {self.species}')
        self.status = result
        self.make_text(system,energies,causamortis)
        system.remove(self)
    
    def write(self):
        return ''.join(self._report_lines)

    @property
    def report(self):
        """Backward-compatible text view of the buffered particle report."""
        return ''.join(self._report_lines)

    @report.setter
    def report(self, value):
        self._report_lines = [value] if value else []


        
 
class Electron(Particles):
    def __init__(self,initial):
        Particles.__init__(self,'electron',initial) 
        self.charge = -1
        self.color  = "red"
        self.marker = "$e^-$"

class Singlet(Particles):
    def __init__(self,initial):
        Particles.__init__(self,'singlet',initial) 
        self.charge = 0
        self.color  = "orange"
        self.marker = "$S_1$"
        
class Triplet(Particles):
    def __init__(self,initial):
        Particles.__init__(self,'triplet',initial) 
        self.charge = 0
        self.color  = "green"
        self.marker = "$T_1$"

class Hole(Particles):
    def __init__(self,initial):
        Particles.__init__(self,'hole',initial) 
        self.charge = 1
        self.color  = "blue"
        self.marker = "$h^+$"

class Vacancy(Particles):
    def __init__(self,initial):
        Particles.__init__(self,'vacancy',initial) 
        self.charge = 0
        self.color  = "black"
        self.marker = "$V$"

class Interstitial(Particles):
    def __init__(self,initial):
        Particles.__init__(self,'interstitial',initial) 
        self.charge = 0
        self.color  = "purple"
        self.marker = "$I$"        

class Frenkelpair(Particles):
    def __init__(self,initial):
        Particles.__init__(self,'frenkelpair',initial) 
        self.charge = 0
        self.color  = "brown"
        self.marker = "$FP$"
        self.ghost_site = None
        self.origin_site = None

class Frenkelpair0(Particles):
    def __init__(self,initial):
        Particles.__init__(self,'frenkelpair0',initial) 
        self.charge = 0
        self.color  = "yellow"
        self.marker = "$FP0$"
        self.ghost_site = None
        self.origin_site = None

class Frenkelpair2(Particles):
    def __init__(self,initial):
        Particles.__init__(self,'frenkelpair2',initial) 
        self.charge = 0
        self.color  = "gray"
        self.marker = "$FP2$"
        self.ghost_site = None
        self.origin_site = None

class I2(Particles): #two interstitial forming I2
    def __init__(self,initial):
        Particles.__init__(self,'I2',initial) 
        self.charge = 0
        self.color  = "magenta"
        self.marker = "$I2$"
        self.ghost_site = None
        self.origin_site = None

class Ghost(Particles):
    def __init__(self,initial):
        Particles.__init__(self,'ghost',initial) 
        self.charge = 0
        self.color  = None
        self.marker = None
