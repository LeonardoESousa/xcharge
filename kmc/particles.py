import numpy as np
import random
import sys

class Particles:
    def __init__(self,species,initial):
        self.species  = species
        self.initial  = initial
        self.position = initial
        self.status   = 'alive'
        self.identity = random.uniform(0,5)
        self.report   = ''
        self.process  = None
        self.destination  = None
        
    def move(self,local):
        self.position = local
        
    
    def make_text(self,system,energies,causamortis):
        time = system.time   
        X,Y,Z = system.X, system.Y, system.Z  
        Mats  = system.mats 
        x0,y0,z0 = X[self.initial],Y[self.initial],Z[self.initial]
        x, y, z  = X[self.position],Y[self.position],Z[self.position]
        dx = x-x0 
        dy = y-y0 
        dz = z-z0 
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
        self.status = 'dead'
        self.make_text(system,energies,causamortis)
        Particula = getattr(sys.modules[__name__], newkind.title())
        system.set_particles([Particula(self.position)])
        system.remove(self)

    def write(self):
        return self.report


        
 
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