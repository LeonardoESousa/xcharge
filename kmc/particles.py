import numpy as np
import random


epsilon_vaccum = 8.85e-22        #Permitivity in C/VAngstrom
e              = -1.60217662e-19 #Electron charge    
kb             = 8.617e-5        # Boltzmann constant
hbar           = 6.582e-16       #Reduced Planck's constant

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
        time = system.time   
        X,Y,Z = system.X, system.Y, system.Z  
        Mats  = system.mats 
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
