"""Indexed event storage for the exact KMC direct method.

The classes in this module do not know about particular particle species or
rate models.  They only require a callback that builds the event catalogue of
one particle, which keeps the scheduler usable by future mechanisms.
"""

from dataclasses import dataclass

import numpy as np


@dataclass(slots=True)
class ParticleEvents:
    """All currently available events for one particle."""

    cumulative_rates: np.ndarray
    processes: tuple
    destinations: np.ndarray
    total: float

    @classmethod
    def empty(cls):
        return cls(
            cumulative_rates=np.empty(0, dtype=float),
            processes=(),
            destinations=np.empty(0, dtype=np.int32),
            total=0.0,
        )

    def select(self, rng):
        if self.total <= 0.0:
            raise RuntimeError("Cannot select from an empty event catalogue.")
        target = rng.random() * self.total
        index = int(np.searchsorted(self.cumulative_rates, target, side="right"))
        return self.processes[index], int(self.destinations[index])


class FenwickTree:
    """Dynamically sized Fenwick tree of non-negative propensities."""

    def __init__(self, capacity=16):
        self.capacity = max(1, int(capacity))
        self.values = np.zeros(self.capacity, dtype=float)
        self.tree = np.zeros(self.capacity + 1, dtype=float)
        self.total = 0.0

    def _rebuild(self):
        self.tree.fill(0.0)
        for index, value in enumerate(self.values):
            if value:
                tree_index = index + 1
                while tree_index <= self.capacity:
                    self.tree[tree_index] += value
                    tree_index += tree_index & -tree_index
        self.total = float(self.values.sum())

    def grow(self, capacity):
        capacity = int(capacity)
        if capacity <= self.capacity:
            return
        values = np.zeros(capacity, dtype=float)
        values[: self.capacity] = self.values
        self.capacity = capacity
        self.values = values
        self.tree = np.zeros(capacity + 1, dtype=float)
        self._rebuild()

    def set(self, index, value):
        value = float(value)
        if not np.isfinite(value) or value < 0.0:
            raise ValueError("KMC propensities must be finite and non-negative.")
        delta = value - self.values[index]
        if delta == 0.0:
            return
        self.values[index] = value
        self.total += delta
        tree_index = index + 1
        while tree_index <= self.capacity:
            self.tree[tree_index] += delta
            tree_index += tree_index & -tree_index

    def find(self, target):
        """Return the index containing ``target`` in ``[0, total)``."""
        if not 0.0 <= target < self.total:
            raise ValueError("Fenwick-tree target is outside the total rate.")
        index = 0
        accumulated = 0.0
        bit = 1 << (self.capacity.bit_length() - 1)
        while bit:
            candidate = index + bit
            if (
                candidate <= self.capacity
                and accumulated + self.tree[candidate] <= target
            ):
                index = candidate
                accumulated += self.tree[candidate]
            bit >>= 1
        return index


class EventScheduler:
    """Exact particle-event scheduler with local propensity updates."""

    def __init__(self, system, event_builder):
        self.system = system
        self.event_builder = event_builder
        capacity = max(16, 2 * len(system.particles))
        self.tree = FenwickTree(capacity)
        self.slots = [None] * capacity
        self.free_slots = list(range(capacity - 1, -1, -1))
        self.particle_slots = {}
        self.catalogues = {}
        for particle in system.particles:
            self.register(particle)

    @property
    def total_rate(self):
        return self.tree.total

    def _grow(self):
        old_capacity = len(self.slots)
        new_capacity = max(16, 2 * old_capacity)
        self.tree.grow(new_capacity)
        self.slots.extend([None] * (new_capacity - old_capacity))
        self.free_slots.extend(range(new_capacity - 1, old_capacity - 1, -1))

    def register(self, particle):
        if particle in self.particle_slots:
            self.refresh_particle(particle)
            return
        if not self.free_slots:
            self._grow()
        slot = self.free_slots.pop()
        self.slots[slot] = particle
        self.particle_slots[particle] = slot
        self.refresh_particle(particle)

    def unregister(self, particle):
        slot = self.particle_slots.pop(particle, None)
        self.catalogues.pop(particle, None)
        if slot is None:
            return
        self.tree.set(slot, 0.0)
        self.slots[slot] = None
        self.free_slots.append(slot)

    def refresh_particle(self, particle):
        slot = self.particle_slots.get(particle)
        if slot is None:
            return
        catalogue = self.event_builder(particle, self.system)
        self.catalogues[particle] = catalogue
        self.tree.set(slot, catalogue.total)

    def refresh_all(self):
        for particle in tuple(self.particle_slots):
            self.refresh_particle(particle)

    def apply_changes(self, changes):
        removed, added, changed_sites, dirty_particles, charge_changed = changes
        for particle in removed:
            self.unregister(particle)
        for particle in added:
            self.register(particle)

        if charge_changed:
            self.refresh_all()
            return

        dirty = set(dirty_particles)
        if changed_sites:
            affected_origins = self.system.affected_origins(changed_sites)
            for site in affected_origins:
                dirty.update(self.system.particles_at_site(site))
        for particle in dirty:
            if particle in self.particle_slots:
                self.refresh_particle(particle)

    def select(self, target):
        slot = self.tree.find(target)
        particle = self.slots[slot]
        if particle is None:  # pragma: no cover - guards internal consistency
            raise RuntimeError("Selected an empty KMC event slot.")
        catalogue = self.catalogues[particle]
        process, destination = catalogue.select(self.system.rng)
        return particle, process, destination
