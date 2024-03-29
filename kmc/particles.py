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
        self.texto    = ''
        self.process  = None
        self.destination  = None
        self.Dx = 0
        self.Dy = 0
        self.Dz = 0
        
    def move(self,local,system):
        Dx, Dy, Dz = system.distance(system,self.position,local)
        self.Dx += Dx
        self.Dy += Dy
        self.Dz += Dz
        self.position = local
        

    def make_text(self,system,energies,causamortis):
        X,Y,Z = system.X, system.Y, system.Z  
        Mats  = system.mats 
        x, y, z  = X[self.position],Y[self.position],Z[self.position]
        mat = Mats[self.position]
        energy = energies[self.position]
        self.texto = f'TEMPO,{self.Dx:.0f},{self.Dy:.0f},{self.Dz:.0f},{self.species},{energy:.2f},{mat:.0f},{x:.0f},{y:.0f},{z:.0f},{causamortis},{self.status}'
        
    def stamp_time(self,system):
        if self.texto != '':
            texto = self.texto
            texto = texto.replace('TEMPO',f'{system.time:.0f}')
            self.report += texto+'\n'
            self.texto = ''
        
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