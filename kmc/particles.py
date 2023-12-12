import random

class Particles:
    def __init__(self, species, initial):
        self.species = species
        self.initial = initial
        self.position = initial
        self.status = "alive"
        self.identity = random.uniform(0, 5)
        self.report = ""
        self.process = None
        self.destination = None
        self.Dx = 0
        self.Dy = 0
        self.Dz = 0

    def move(self, local, system):
        Dx, Dy, Dz = system.distance(system, self.position, local)
        self.Dx += Dx
        self.Dy += Dy
        self.Dz += Dz
        self.position = local

    def make_text(self, system, energy, causamortis):
        X, Y, Z = system.X, system.Y, system.Z
        Mats = system.mats
        x, y, z = X[self.position], Y[self.position], Z[self.position]
        mat = Mats[self.position]
        self.report += f"{system.time:.0f},{self.Dx:.0f},{self.Dy:.0f},{self.Dz:.0f},{self.species},{energy:.2f},{mat:.0f},{x:.0f},{y:.0f},{z:.0f},{causamortis},{self.status}"
        self.report += "\n"
        
    def kill(self, causamortis, system, energy, result):
        self.status = result
        self.make_text(system, energy, causamortis)
        system.remove(self)

    def write(self):
        return self.report


class Singlet(Particles):
    def __init__(self, initial):
        Particles.__init__(self, "singlet", initial)
        self.charge = 0
        self.color = "orange"
        self.marker = "$S_1$"
        self.conformer = None


class Triplet(Particles):
    def __init__(self, initial):
        Particles.__init__(self, "triplet", initial)
        self.charge = 0
        self.color = "green"
        self.marker = "$T_1$"
        self.conformer = None


class Electron(Particles):
    def __init__(self, initial):
        Particles.__init__(self, "electron", initial)
        self.charge = -1
        self.color = "red"
        self.marker = "$e^-$"


class Hole(Particles):
    def __init__(self, initial):
        Particles.__init__(self, "hole", initial)
        self.charge = 1
        self.color = "blue"
        self.marker = "$h^+$"

class Ghost(Particles):
    def __init__(self, initial):
        Particles.__init__(self, "ghost", initial)
        self.charge = 0
        self.color = "black"
        self.marker = "$G$"
        