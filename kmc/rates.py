# pylint: disable=no-member
# pylint: disable=import-error
# pylint: disable=no-name-in-module

import numpy as np
import kmc.particles
import kmc.utils
import kmc.variables
import random

EPSILON_0 = kmc.variables.EPSILON_0  # Permitivity in C/Vm
E = kmc.variables.E  # Electron charge
KB = kmc.variables.KB  # Boltzmann constant
HBAR = kmc.variables.HBAR  # Reduced Planck's constant


###RATES#################################################################################

##FUNCTION FOR SETTING RADII#############################################################


def raios(num, Rf, mat, lifetime, mats):
    # Initialize the Raios array with the value of Rf[(mat,mat)]
    Raios = np.empty(num)
    Raios.fill(Rf[(mat, mat)])

    # Use NumPy's where function to set the values of Raios for the other materials
    for m in lifetime.keys():
        if m != mat:
            Raios = np.where(mats == m, Rf[(mat, m)], Raios)

    return Raios


# function to convert dictionary with (i,j) keys to ixj array
def dict_to_array(d):
    keys = d.keys()
    num_keys = len(set(key[0] for key in keys))
    radius = np.empty((num_keys, num_keys))
    for key in keys:
        radius[key[0], key[1]] = d[key]
    return radius


#########################################################################################


##STANDARD FORSTER TRANSFER RATE#########################################################
class Forster:
    def __init__(self, **kwargs):
        self.kind = "jump"
        self.Rf = dict_to_array(kwargs["Rf"])
        self.lifetime = kwargs["life"]
        self.mu = kwargs["mu"]
        self.alpha = 1.15 * 0.53

    def rate(self, **kwargs):
        r = kwargs["r"]
        mats = kwargs["mats"]
        mat = kwargs["matlocal"]
        num = len(mats)
        taxa = kmc.utils.forster(
            self.Rf[mat, :],
            mats,
            num,
            self.alpha * self.mu[mat],
            r,
            1 / self.lifetime[mat],
        )
        return taxa

    def action(self, particle, system, local):
        particle.move(local, system)


#########################################################################################


##TRIPLET TO SINGLET FORSTER TRANSFER####################################################
class ForsterT:
    def __init__(self, **kwargs):
        self.kind = "jump"
        self.Rf = dict_to_array(kwargs["Rf"])
        self.lifetime = kwargs["life"]
        self.mu = kwargs["mu"]
        self.alpha = 1.15 * 0.53

    def rate(self, **kwargs):
        r = kwargs["r"]
        mats = kwargs["mats"]
        mat = kwargs["matlocal"]
        num = len(mats)
        taxa = kmc.utils.forster(
            self.Rf[mat, :],
            mats,
            num,
            self.alpha * self.mu[mat],
            r,
            1 / self.lifetime[mat],
        )
        return taxa

    def action(self, particle, system, local):
        particle.move(local, system)
        particle.kill("tts", system, system.t1, "converted")
        system.set_particles([kmc.particles.Singlet(local)])


#########################################################################################


##FORSTER ANNIHILATION RADIUS#########################################################
class ForsterAnniRad:
    def __init__(self, **kwargs):
        self.kind = "jump"
        self.Rf = dict_to_array(kwargs["Rf"])
        self.lifetime = kwargs["life"]
        self.mu = kwargs["mu"]
        self.alpha = 1.15 * 0.53
        self.anni_rad = kwargs["anni_rad"]

    def rate(self, **kwargs):
        r = kwargs["r"]
        system = kwargs["system"]
        ex = kwargs["particle"]
        mats = kwargs["mats"]
        cut = kwargs["cut"]
        mat = kwargs["matlocal"]
        num = len(mats)
        relevant_particles = [
            p
            for p in system.particles
            if p.identity != ex.identity and p.position in cut
        ]
        ss = [
            (
                np.where(cut == p.position)[0][0],
                self.anni_rad[(mat, system.mats[p.position])][p.species],
            )
            for p in relevant_particles
        ]
        replace_pos = np.array([ele[0] for ele in ss], dtype=np.int32)
        replace_raios = np.array([ele[1] for ele in ss], dtype=np.double)
        mum = len(replace_pos)
        taxa = kmc.utils.forster_anni(
            self.Rf[mat, :],
            mats,
            num,
            self.alpha * self.mu[mat],
            r,
            1 / self.lifetime[mat],
            replace_pos,
            replace_raios,
            mum,
        )
        return taxa

    def action(self, particle, system, local):
        particle.move(local, system)


#########################################################################################


##STANDARD DEXTER TRANSFER RATE##########################################################
class Dexter:
    def __init__(self, **kwargs):
        self.kind = "jump"
        self.Rd = dict_to_array(kwargs["Rd"])
        self.lifetime = kwargs["life"]
        self.L = kwargs["L"]

    def rate(self, **kwargs):
        r = kwargs["r"]
        mats = kwargs["mats"]
        mat = kwargs["matlocal"]
        num = len(mats)
        taxa = kmc.utils.dexter(
            self.Rd[mat, :], 1 / self.L[mat], 1 / self.lifetime[mat], mats, num, r
        )
        return taxa

    def action(self, particle, system, local):
        particle.move(local, system)


#########################################################################################


##STANDARD DEXTER TRANSFER RATE##########################################################
class Marcus:
    def __init__(self, **kwargs):
        self.kind = "jump"
        self.coupling = dict_to_array(kwargs["coupling"])
        self.reorg = kwargs["reorg"]
        self.level = kwargs["level"]
        self.kbt = kwargs["temperature"] * KB
        self.decay = kwargs["decay"]
        self.prefactor = 1e-12 * 2 * np.pi / HBAR
        
    def rate(self, **kwargs):
        system = kwargs["system"]
        cut = kwargs["cut"]
        r = kwargs["r"]
        particle = kwargs["particle"]
        mats = kwargs["mats"]
        mat = kwargs["matlocal"]
        num = len(mats)
        energy = getattr(system, self.level, None)
        site_energy = energy[particle.position]
        taxa = kmc.utils.marcus(
            self.coupling[mat, :],
            energy[cut],
            self.reorg[mat],
            self.prefactor,
            mats,
            num,
            site_energy,
            self.kbt,
            self.decay,
            r,
        )
        return taxa

    def action(self, particle, system, local):
        particle.move(local, system)

#########################################################################################
def gauss(x_value, mean, std):
    y_value = (1 / (np.sqrt(2 * np.pi) * std)) * np.exp(-0.5 * ((x_value - mean) / std) ** 2)
    return y_value

# Speed of light
C = 299792458 # m/s
C *= 1e10 # Å/s
CONST = (HBAR**3) * (9 * (C**4)) / (8 * np.pi)
CONST *= 1e-12

def radius(X, YD, YA,kappa2):
    # Calculates the overlap
    X4 = X*X
    X4 = X4*X4
    Overlap = YA * YD / (X4)

    # Integrates overlap
    IntOver = np.trapz(Overlap, X)

    # Calculates radius sixth power
    radius6 = kappa2 * CONST * IntOver
    return radius6


#########################################################################################

class DynamicForster:
    def __init__(self, **kwargs):
        self.kind = "jump"
        self.excited = kwargs["excited"]
        self.keys = list(self.excited.keys())
        self.ground = kwargs["ground"]
        self.kappa = kwargs["kappa"]
        self.x = np.linspace(0.01, 10, 500)

    def rate(self, **kwargs):
        r = kwargs["r"]
        system = kwargs["system"]
        static = system.static
        particle = kwargs["particle"]
        if particle.conformer is None:
            particle.conformer = random.choice(self.keys)
        engs_s1, sigma_s1, diff_rate  = self.excited[particle.conformer]
        engs_s1 += static[particle.position]
        emission = diff_rate*gauss(self.x, engs_s1, sigma_s1)
        cut = kwargs["cut"]
        acceptors = system.s0[cut]
        static = static[cut]
        taxa = np.zeros(len(acceptors))
        for j, dist in enumerate(r):
            if dist != 0:
                engs_s0, sigma_s0, cross_section = self.ground[acceptors[j]]
                engs_s0 += static[j]
                absorption = cross_section*gauss(self.x, engs_s0, sigma_s0)
                dist6 = dist*dist*dist
                dist6 = dist6*dist6
                taxa[j] = radius(self.x, absorption, emission ,self.kappa)/dist6
        return taxa

    def action(self, particle, system, local):
        particle.move(local, system)


# MONOMOLECULAR RATES#####################################################################

class DynamicFluor:
    def __init__(self, **kwargs):
        self.kind = "fluor"
        self.excited = kwargs["excited"]
        self.keys = list(self.excited.keys())

    def rate(self, **kwargs):
        particle = kwargs["particle"]
        if particle.conformer is None:
            particle.conformer = random.choice(self.keys)
        particle = kwargs["particle"]
        _, _, diff_rate  = self.excited[particle.conformer]
        return 1e-12*diff_rate/HBAR

    def action(self, particle, system, local):
        site_energy = system.static[local]
        mean_energy = site_energy + self.excited[particle.conformer][0]
        #sample from gaussian distribution with np.random.normal(mean, std)
        particle_energy = np.random.normal(mean_energy, self.excited[particle.conformer][1])
        particle.kill(self.kind, system, particle_energy, "dead")

##INTERSYSTEM CROSSING RATE##############################################################
class DynamicISC:
    def __init__(self, **kwargs):
        self.kind = "isc"
        self.map = {"singlet": "triplet", "triplet": "singlet"}
        self.excited = kwargs["excited"]
        self.keys = list(self.excited.keys())

    def rate(self, **kwargs):
        particle = kwargs["particle"]
        if particle.conformer is None:
            particle.conformer = random.choice(self.keys)
        particle = kwargs["particle"]
        _, _, diff_rate  = self.excited[particle.conformer]
        return diff_rate

    def action(self, particle, system, local):
        site_energy = system.static[local]
        mean_energy = site_energy + self.excited[particle.conformer][0]
        #sample from gaussian distribution with np.random.normal(mean, std)
        particle_energy = np.random.normal(mean_energy, self.excited[particle.conformer][1])
        if particle.species == "singlet":
            system.set_particles([kmc.particles.Triplet(particle.position)])
            particle.kill(self.kind, system, particle_energy, "converted")
        elif particle.species == "triplet":
            system.set_particles([kmc.particles.Singlet(particle.position)])
            particle.kill("r" + self.kind, system, particle_energy, "converted")


#########################################################################################

class NonAdiabatic:
    def __init__(self, **kwargs):
        self.kind = 'nonadiabatic'
        self.frequency = kwargs['frequency']
        self.conformers = kwargs['conformers']

    def rate(self, **kwargs):
        return self.frequency

    def action(self, particle, system, local):
        particle.conformer = random.choice(self.conformers)

class NonAdiabaticGround:
    def __init__(self, **kwargs):
        self.kind = 'nonadiabatic'
        self.frequency = kwargs['frequency']
        self.conformers = kwargs['conformers']

    def rate(self, **kwargs):
        return self.frequency

    def action(self, particle, system, local):
        if len(system.particles) == 1:
            system.remove(particle)
        else:
            self.conformers.assign_to_system(system)

##FLUORESCENCE RATE######################################################################
class Fluor:
    def __init__(self, **kwargs):
        self.kind = "fluor"
        self.lifetime = kwargs["life"]

    def rate(self, **kwargs):
        return 1 / self.lifetime[kwargs["material"]]

    def action(self, particle, system, local):
        particle.kill(self.kind, system, system.s1[local], "dead")


#########################################################################################


##PHOSPHORESCENCE RATE###################################################################
class Phosph:
    def __init__(self, **kwargs):
        self.kind = "phosph"
        self.lifetime = kwargs["life"]

    def rate(self, **kwargs):
        return 1 / self.lifetime[kwargs["material"]]

    def action(self, particle, system, local):
        particle.kill(self.kind, system, system.t1[local], "dead")


#########################################################################################


##NONRADIATIVE DECAY RATE################################################################
class Nonrad:
    def __init__(self, **kwargs):
        self.kind = "nonrad"
        self.taxa = kwargs["rate"]

    def rate(self, **kwargs):
        return self.taxa[kwargs["material"]]

    def action(self, particle, system, local):
        particle.kill(self.kind, system, system.s1[local], "dead")


#########################################################################################


##INTERSYSTEM CROSSING RATE##############################################################
class ISC:
    def __init__(self, **kwargs):
        self.kind = "isc"
        self.taxa = kwargs["rate"]
        self.map = {"singlet": "triplet", "triplet": "singlet"}

    def rate(self, **kwargs):
        material = kwargs["material"]
        return self.taxa[material]

    def action(self, particle, system, local):
        if particle.species == "singlet":
            system.set_particles([kmc.particles.Triplet(particle.position)])
            particle.kill(self.kind, system, system.s1[local], "converted")
        elif particle.species == "triplet":
            system.set_particles([kmc.particles.Singlet(particle.position)])
            particle.kill("r" + self.kind, system, system.s1[local], "converted")


#########################################################################################
