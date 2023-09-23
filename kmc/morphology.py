import sys
import random
import numpy as np
from kmc.particles import Singlet, Triplet, Electron, Hole

# SET OF FUNCS THAT GENERATE THE MORPHOLOGY OF THE SYSTEM
# note: always define the function by list (param) that contains the things needed
#### CHOOSE A FUNC TO GENERATE PARTICLES


def randomized(available_sites, num_sites, system, kwargs):
    acceptable_materials = kwargs["mat"]
    materials = system.mats
    # get indices of materials that are contained in list of acceptable materials
    acceptable_indices = np.where(np.isin(materials, acceptable_materials))[0]
    available_sites = np.array(available_sites)[acceptable_indices]
    # Initially select random sites
    selected_sites = random.choices(available_sites, k=num_sites)
    # Convert the set of selected sites back to a list before returning
    return selected_sites


##CLASS FOR GENERATING PARTICLES IN THE SYSTEM###########################################
class CreateParticles:
    def __init__(self, kind, prob, num, method, **kwargs):
        self.kind = kind
        prob = np.array(prob)
        prob = prob / np.sum(prob)
        prob = np.cumsum(prob)
        self.prob = prob
        self.num = num
        self.method = method
        self.argv = kwargs

    def assign_to_system(self, system):
        selected = self.method(range(len(system.X)), self.num, system, self.argv)
        for number in selected:
            kind = self.kind[np.where(random.uniform(0, 1) <= self.prob)[0][0]]
            Particula = getattr(sys.modules[__name__], kind.title())
            particle = Particula(number)
            system.set_particles([particle])


#########################################################################################


##CLASSES TO ASSIGN ENERGIES TO LATTICE##################################################
class GaussianEnergy:
    def __init__(self, s1s):
        self.s1s = s1s

    def assign_to_system(self, system):
        type_en = self.s1s["level"]
        uniq = system.uniq
        N = len(system.mats)

        means = np.empty(N)
        stds = np.empty(N)
        means.fill(self.s1s[uniq[0]][0])
        stds.fill(self.s1s[uniq[0]][1])
        for m in uniq[1:]:
            mask = system.mats == m
            means[mask] = self.s1s[m][0]
            stds[mask] = self.s1s[m][1]
        s1 = np.random.normal(means, stds, N)
        system.set_energies(s1, type_en)


#########################################################################################


class Lattice:
    """Arguments:
    num_sites:   number of sites in the lattice (integer)
    vector:      3-element list with distance between sites in x,y and z directions in Å (list)
    disorder:    3-element list with standard deviation of distances between sites in x,y and z directions in Å (list)
    composition: n-element list with the proportion of each of the n different materials in the lattice  (list)
    """

    def __init__(
        self, num_sites, vector, disorder, composition
    ):  # initializing the lattice class with some basic info given by the user
        self.num_sites = num_sites
        self.vector = vector
        self.disorder = disorder
        self.composition = np.cumsum([i / np.sum(composition) for i in composition])

    def make(self):  # Generating the set X,Y,Z,Mats
        dim = []
        for elem in self.vector:
            if elem != 0:
                dim.append(1)
            else:
                dim.append(0)

        num = int(self.num_sites ** (1 / np.sum(dim)))
        numx = max(dim[0] * num, 1)
        numy = max(dim[1] * num, 1)
        numz = max(dim[2] * num, 1)
        Numx = np.array(range(numx)) * self.vector[0]
        Numy = np.array(range(numy)) * self.vector[1]
        Numz = np.array(range(numz)) * self.vector[2]

        total = numx * numy * numz
        X, Y, Z = np.meshgrid(Numx, Numy, Numz, indexing="ij")
        X, Y, Z = (
            X.reshape(
                total,
            ),
            Y.reshape(
                total,
            ),
            Z.reshape(
                total,
            ),
        )
        # add noise to the lattice
        X = X + np.random.normal(0, self.disorder[0], total)
        Y = Y + np.random.normal(0, self.disorder[1], total)
        Z = Z + np.random.normal(0, self.disorder[2], total)

        luck = np.random.uniform(0, 1, total)
        Mats = np.zeros(total)
        for i in reversed(range(len(self.composition))):
            Mats[luck < self.composition[i]] = i
        return X, Y, Z, Mats

    def assign_to_system(self, system):  # adding the X,Y,Z,Mats to the system
        X, Y, Z, Mats = self.make()
        system.set_morph(X, Y, Z, Mats)
