"""Vacancy and interstitial migration toward periodic grain boundaries.

Run from the repository root with

    python -m kmc input_examples/migration_barrier.py

The barrier function returns a destination-dependent activation-barrier
change in eV. Negative values increase the migration rate.
"""

import numpy as np

import kmc.morphology as morphology
from kmc.potentials import MigrationBarrier
from kmc.rates import Migration


# BASIC PARAMETERS
identifier = "migration_barrier"
time_units = "s"
time_limit = 2.0e-8
rounds = 4
n_proc = 1
frozen_lattice = True
periodic = True
bimolec = False
particle_condition = True
random_seed = 12345
append_output = False

# Migration uses the first-neighbor graph, not this general cutoff graph.
cutoff = 0.1
first_neighbor_rtol = 0.05

animation_mode = False
save_animation = False


# PERIODIC MIGRATION-BARRIER PROFILE
temperature = 298.0


def grain_boundary_barrier(x, period):
    """Return ΔE‡(x) in eV for boundaries at x=0 and x=period."""
    distance_to_boundary = np.minimum(x, period - x)
    depth = 0.03  # eV; negative correction at the boundary
    width = 3.0   # Angstrom
    return -depth * np.exp(-(distance_to_boundary / width) ** 2)


migration_barrier = MigrationBarrier(
    function=grain_boundary_barrier,
    temperature=temperature,
    axis="x",
)


# MIGRATION PROCESSES
# V and I can retain different unmodified Eyring rates. The one registered
# migration_barrier profile is applied automatically to both processes.
vacancy_migration = Migration(k={(0, 0): 2.0e8})      # s^-1
interstitial_migration = Migration(k={(0, 0): 5.0e8})  # s^-1

processes = {
    "vacancy": [vacancy_migration],
    "interstitial": [interstitial_migration],
}
monomolecular = {
    "vacancy": [],
    "interstitial": [],
}
generation = []


# PERIODIC ORTHORHOMBIC LATTICE
# 4 x 4 x 4 sites with a cell length of 12 Å along each direction.
lattice_func = morphology.Lattice(
    num_sites=64,
    vector=[3.0, 3.0, 3.0],
    disorder=[0.0, 0.0, 0.0],
    composition=[1.0],
)

# Particle reporting expects an S1 energy array even for defect-only runs.
s1_energy = morphology.Gaussian_energy(
    {0: (0.0, 0.0), "level": "s1"}
)


# INITIAL PARTICLES
vacancy = morphology.Create_Particles(
    "vacancy",
    1,
    morphology.randomized,
    mat=[0],
)
interstitial = morphology.Create_Particles(
    "interstitial",
    1,
    morphology.randomized,
    mat=[0],
)
