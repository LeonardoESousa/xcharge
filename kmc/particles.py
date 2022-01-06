import numpy as np
import random

class Particles:
    def __init__(self,species,initial):
        self.species  = species
        self.initial  = initial
        self.position = initial
        self.status   = 'alive'
        self.identity = random.uniform(0,5)
        self.report   = ''
    
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
        texto = '{0:<15.6e}  {1:<10.2f}  {2:<10.2f}  {3:<10.2f}  {4:<10}  {5:<10.4f}  {6:<10}  {7:<10.2f}  {8:<10.2f}  {9:<10.2f}  {10:<10}'.format(
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
        self.color  = "red"
        self.marker = "$e^-$"

class Exciton(Particles):
    def __init__(self,kind,initial):
        Particles.__init__(self,kind,initial) 
        self.charge = 0
        if kind == 'singlet':
            self.color  = "orange"
            self.marker = "$S_1$"
        elif kind == 'triplet':
            self.color  = "green"
            self.marker = "$T_1$"

class Hole(Particles):
    def __init__(self,initial):
        Particles.__init__(self,'hole',initial) 
        self.charge = 1
        self.color  = "blue"
        self.marker = "$h^+$"