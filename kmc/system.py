# pylint: disable=no-member
# pylint: disable=import-error
# pylint: disable=no-name-in-module
import numpy as np
import kmc.utils
import kmc.variables

EPSILON_0 = kmc.variables.EPSILON_0  # Permitivity in C/Vm
E = kmc.variables.E  # Electron charge
KB = kmc.variables.KB  # Boltzmann constant
HBAR = kmc.variables.HBAR  # Reduced Planck's constant


class System:
    def __init__(self):
        self.dead = []
        self.time = 0
        self.potential_time = -1
        self.IT = 0  # number of steps

    def set_morph(self, X, Y, Z, Mats):
        self.X = X
        self.Y = Y
        self.Z = Z
        self.Lx = max(X) - min(X)
        self.Ly = max(Y) - min(Y)
        self.Lz = max(Z) - min(Z)
        self.R = np.hstack((X[:, np.newaxis], Y[:, np.newaxis], Z[:, np.newaxis]))
        Mats = np.array(Mats, dtype=np.int32)
        self.mats = Mats
        self.uniq = np.unique(Mats)

    def set_basic_info(
        self,
        monomolecular,
        processes,
        identifier,
        animation_mode,
        time_limit,
        pause,
        anni,
        distance,
    ):
        self.processes = processes
        self.monomolecular = monomolecular
        self.identifier = identifier
        self.animation_mode = animation_mode
        self.time_limit = time_limit
        self.pause = pause
        self.bimolec = anni
        self.distance = distance

    def set_particles(self, Ss):
        try:
            self.particles += Ss
        except (TypeError,AttributeError):
            self.particles = Ss

    def reset_particles(self):
        self.particles = None

    def count_particles(self):
        return len(self.particles)

    def get_num(self):
        return len(self.X)

    def set_energies(self, energy, kind):
        setattr(self, kind.lower(), energy)

    def remove(self, particle):
        self.particles.remove(particle)
        self.dead.append(particle)
