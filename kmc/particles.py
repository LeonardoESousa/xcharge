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
        self.Dx = 0
        self.Dy = 0
        self.Dz = 0
        
    def move(self,local,system):
        Dx, Dy, Dz = system.distance(system,self.position)
        self.Dx += Dx[local]
        self.Dy += Dy[local]
        self.Dz += Dz[local]
        self.position = local
        

    def make_text(self,system,energies,causamortis):
        time = system.time   
        X,Y,Z = system.X, system.Y, system.Z  
        Mats  = system.mats 
        x, y, z  = X[self.position],Y[self.position],Z[self.position]
        mat = Mats[self.position]
        energy = energies[self.position]
        texto = '{0:<15.6e}  {1:<10.2f}  {2:<10.2f}  {3:<10.2f}  {4:<10}  {5:<10.4f}  {6:<10}  {7:<10.2f}  {8:<10.2f}  {9:<10.2f}  {10:<11}  {11:<10}'.format(
                       time,self.Dx,self.Dy,self.Dz,self.species,energy,mat,x,y,z,causamortis,self.status)
        self.report += texto+'\n'
    
    def kill(self,causamortis,system,energies,result):
        self.status = result
        self.make_text(system,energies,causamortis)
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