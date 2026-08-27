from collections import Counter
from itertools import product

import numpy as np
from scipy.spatial import cKDTree

import kmc.utils


EPSILON_VACUUM = 8.854187e-12
ELEMENTARY_CHARGE = 1.60217662e-19


class System:
    def __init__(self):
        self.dead = []
        self.particles = []
        self.time = 0.0
        self.potential_time = -1.0
        self.IT = 0
        self.positions_by_species = {}
        self.occupancy = Counter()
        self.cell = None
        self._cutoff_neighbors = None
        self._first_neighbors = None
        self._reverse_neighbors = None
        self.rate_cache = {}
        self.rng = np.random.default_rng()

    def set_morph(self, X, Y, Z, Mats, cell=None):
        self.X = np.asarray(X, dtype=float)
        self.Y = np.asarray(Y, dtype=float)
        self.Z = np.asarray(Z, dtype=float)
        self.R = np.column_stack((self.X, self.Y, self.Z))
        self.mats = np.asarray(Mats, dtype=np.int32)
        self.uniq = np.unique(self.mats)
        if cell is not None:
            self.set_cell(cell)
        else:
            self.Lx = float(np.ptp(self.X))
            self.Ly = float(np.ptp(self.Y))
            self.Lz = float(np.ptp(self.Z))

    def set_cell(self, cell):
        cell = np.asarray(cell, dtype=float)
        if cell.shape != (3, 3):
            raise ValueError("The periodic cell must be a 3 x 3 matrix.")
        self.cell = cell
        self.Lx, self.Ly, self.Lz = np.linalg.norm(cell, axis=1)

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
        generation,
        cutoff,
        periodic,
        first_neighbor_rtol,
        particle_condition,
        creation_threshold,
        time_units,
    ):
        self.processes = processes
        self.monomolecular = monomolecular
        self.identifier = identifier
        self.animation_mode = animation_mode
        self.time_limit = float(time_limit)
        self.pause = pause
        self.bimolec = anni
        self.distance = distance
        self.generation = generation
        self.cutoff = float(cutoff)
        self.periodic = bool(periodic)
        self.first_neighbor_rtol = float(first_neighbor_rtol)
        self.particle_condition = bool(particle_condition)
        self.creation_threshold = float(creation_threshold)
        self.time_units = time_units

    def set_rng(self, seed):
        self.rng = np.random.default_rng(seed)

    def _translation_vectors(self):
        if not self.periodic:
            return np.zeros((1, 3))
        if self.cell is None:
            raise ValueError(
                "Periodic simulations require the full lattice cell. "
                "Use ReadCIF or provide a cell-aware morphology."
            )
        active = [row for row in self.cell if np.linalg.norm(row) > 0]
        shifts = []
        for coefficients in product((-1, 0, 1), repeat=len(active)):
            shift = np.zeros(3)
            for coefficient, vector in zip(coefficients, active):
                shift += coefficient * vector
            shifts.append(shift)
        return np.asarray(shifts)

    @staticmethod
    def _deduplicate_neighbors(
        origin,
        image_indices,
        image_sites,
        image_coordinates,
        coordinates,
        exclude_origin=True,
    ):
        selected = {}
        for image_index in image_indices:
            destination = int(image_sites[image_index])
            if exclude_origin and destination == origin:
                continue
            displacement = image_coordinates[image_index] - coordinates[origin]
            distance = float(np.linalg.norm(displacement))
            previous = selected.get(destination)
            if previous is None or distance < previous[0]:
                selected[destination] = (distance, displacement)
        ordered = sorted(selected.items())
        if not ordered:
            empty_i = np.empty(0, dtype=np.int32)
            empty_f = np.empty(0, dtype=float)
            return empty_i, empty_f, empty_f.copy(), empty_f.copy(), empty_f.copy()
        destinations = np.fromiter((item[0] for item in ordered), dtype=np.int32)
        distances = np.fromiter((item[1][0] for item in ordered), dtype=float)
        displacements = np.asarray([item[1][1] for item in ordered])
        return (
            destinations,
            displacements[:, 0],
            displacements[:, 1],
            displacements[:, 2],
            distances,
        )

    def build_neighbor_topology(self):
        """Precompute cutoff neighbors and material-resolved first shells."""
        translations = self._translation_vectors()
        n_sites = len(self.R)
        image_coordinates = np.concatenate(
            [self.R + translation for translation in translations], axis=0
        )
        image_sites = np.tile(np.arange(n_sites, dtype=np.int32), len(translations))
        all_tree = cKDTree(image_coordinates)

        self._cutoff_neighbors = []
        for origin, coordinate in enumerate(self.R):
            image_indices = all_tree.query_ball_point(coordinate, self.cutoff)
            self._cutoff_neighbors.append(
                self._deduplicate_neighbors(
                    origin,
                    image_indices,
                    image_sites,
                    image_coordinates,
                    self.R,
                )
            )

        first_by_site = [dict() for _ in range(n_sites)]
        for material in self.uniq:
            material_sites = np.flatnonzero(self.mats == material).astype(np.int32)
            material_images = np.concatenate(
                [self.R[material_sites] + translation for translation in translations],
                axis=0,
            )
            material_image_sites = np.tile(material_sites, len(translations))
            tree = cKDTree(material_images)
            k = min(max(8, 2 * len(translations)), len(material_images))

            for origin, coordinate in enumerate(self.R):
                distances, indices = tree.query(coordinate, k=k)
                distances = np.atleast_1d(distances)
                indices = np.atleast_1d(indices)
                nearest = None
                for distance, image_index in zip(distances, indices):
                    destination = int(material_image_sites[int(image_index)])
                    if destination != origin and distance > 1e-12:
                        nearest = float(distance)
                        break
                if nearest is None:
                    continue
                radius = nearest * (1.0 + self.first_neighbor_rtol) + 1e-10
                image_indices = tree.query_ball_point(coordinate, radius)
                data = self._deduplicate_neighbors(
                    origin,
                    image_indices,
                    material_image_sites,
                    material_images,
                    self.R,
                )
                if data[0].size:
                    first_by_site[origin][int(material)] = data

        self._first_neighbors = []
        for material_data in first_by_site:
            parts = list(material_data.values())
            if not parts:
                empty_i = np.empty(0, dtype=np.int32)
                empty_f = np.empty(0, dtype=float)
                self._first_neighbors.append(
                    (empty_i, empty_f, empty_f.copy(), empty_f.copy(), empty_f.copy())
                )
                continue
            self._first_neighbors.append(
                tuple(np.concatenate([part[i] for part in parts]) for i in range(5))
            )

        reverse = [set() for _ in range(n_sites)]
        for neighborhood in (self._cutoff_neighbors, self._first_neighbors):
            for origin, data in enumerate(neighborhood):
                for destination in data[0]:
                    reverse[int(destination)].add(origin)
        self._reverse_neighbors = reverse

    def neighbor_data(self, site, mode="cutoff"):
        if self._cutoff_neighbors is None or self._first_neighbors is None:
            raise RuntimeError("Neighbor topology has not been initialized.")
        if mode == "first":
            return self._first_neighbors[site]
        if mode == "cutoff":
            return self._cutoff_neighbors[site]
        raise ValueError(f"Unknown neighbor mode: {mode}")

    def first_neighbors_of_material(self, site, material):
        data = self._first_neighbors[site]
        destinations = data[0]
        return destinations[self.mats[destinations] == material]

    def affected_origins(self, changed_sites):
        affected = set(int(site) for site in changed_sites)
        for site in changed_sites:
            affected.update(self._reverse_neighbors[int(site)])
        return affected

    def set_particles(self, particles):
        for particle in particles:
            self.particles.append(particle)
            species_positions = self.positions_by_species.setdefault(
                particle.species, Counter()
            )
            species_positions[particle.position] += 1
            self.occupancy[particle.position] += 1

    def reset_particles(self):
        self.particles = []
        self.dead = []
        self.positions_by_species = {}
        self.occupancy = Counter()
        self.time = 0.0
        self.potential_time = -1.0
        self.IT = 0
        self.rate_cache = {}

    def species_positions(self, species):
        return self.positions_by_species.get(species, {}).keys()

    def set_dipoles(self, mus):
        self.mu = mus
        self.norma_mu = np.sqrt(np.sum(mus * mus, axis=1))
        self.mu /= self.norma_mu[:, np.newaxis]

    def count_particles(self):
        return len(self.particles)

    def set_medium(self, eps_rel):
        self.eps_rel = eps_rel
        self.epsilon = eps_rel * EPSILON_VACUUM

    def get_num(self):
        return len(self.X)

    def set_energies(self, energy, kind):
        setattr(self, kind.lower(), energy)

    @staticmethod
    def _decrement(counter, key):
        counter[key] -= 1
        if counter[key] <= 0:
            del counter[key]

    def remove(self, particle):
        self.particles.remove(particle)
        species_positions = self.positions_by_species.get(particle.species)
        if species_positions is not None:
            self._decrement(species_positions, particle.position)
        self._decrement(self.occupancy, particle.position)
        self.dead.append(particle)

    def set_electric_field(self, field):
        self.field = np.asarray(field, dtype=float) / 1e8

    def minimum_image_displacement(self, local, destination=None):
        if destination is None:
            displacement = self.R - self.R[local]
        else:
            displacement = np.asarray(self.R[destination] - self.R[local])
        if not self.periodic:
            return tuple(np.moveaxis(displacement, -1, 0))
        if self.cell is None:
            raise ValueError("Periodic displacement requested without a lattice cell.")
        best = np.array(displacement, copy=True)
        best_norm = np.sum(best * best, axis=-1)
        for translation in self._translation_vectors():
            candidate = displacement + translation
            norm = np.sum(candidate * candidate, axis=-1)
            mask = norm < best_norm
            best = np.where(np.expand_dims(mask, axis=-1), candidate, best)
            best_norm = np.where(mask, norm, best_norm)
        return tuple(np.moveaxis(best, -1, 0))

    def electrostatic(self):
        if self.time > self.potential_time:
            potential = 0.0
            for particle in self.particles:
                if particle.charge != 0:
                    dx, dy, dz = self.minimum_image_displacement(particle.position)
                    r = kmc.utils.distances(dx, dy, dz, len(dx)) * 1e-10
                    r[r == 0] = np.inf
                    potential += (
                        particle.charge
                        * ELEMENTARY_CHARGE
                        / (4 * np.pi * self.epsilon * r)
                    )
            self.potential = potential
            self.potential_time = self.time
        return self.potential

    def update_position(self, particle, old_position, new_position):
        species_positions = self.positions_by_species[particle.species]
        self._decrement(species_positions, old_position)
        species_positions[new_position] += 1
        self._decrement(self.occupancy, old_position)
        self.occupancy[new_position] += 1
