# pylint: disable=no-member
# pylint: disable=import-error
# pylint: disable=no-name-in-module

import numpy as np
import kmc.particles
import kmc.utils
import kmc.variables

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

# MONOMOLECULAR RATES#####################################################################


##FLUORESCENCE RATE######################################################################
class Fluor:
    def __init__(self, **kwargs):
        self.kind = "fluor"
        self.lifetime = kwargs["life"]

    def rate(self, **kwargs):
        return 1 / self.lifetime[kwargs["material"]]

    def action(self, particle, system, local):
        particle.kill(self.kind, system, system.s1, "dead")


#########################################################################################


##PHOSPHORESCENCE RATE###################################################################
class Phosph:
    def __init__(self, **kwargs):
        self.kind = "phosph"
        self.lifetime = kwargs["life"]

    def rate(self, **kwargs):
        return 1 / self.lifetime[kwargs["material"]]

    def action(self, particle, system, local):
        particle.kill(self.kind, system, system.t1, "dead")


#########################################################################################


##NONRADIATIVE DECAY RATE################################################################
class Nonrad:
    def __init__(self, **kwargs):
        self.kind = "nonrad"
        self.taxa = kwargs["rate"]

    def rate(self, **kwargs):
        return self.taxa[kwargs["material"]]

    def action(self, particle, system, local):
        particle.kill(self.kind, system, system.s1, "dead")


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
            particle.kill(self.kind, system, system.s1, "converted")
        elif particle.species == "triplet":
            system.set_particles([kmc.particles.Singlet(particle.position)])
            particle.kill("r" + self.kind, system, system.s1, "converted")


#########################################################################################
