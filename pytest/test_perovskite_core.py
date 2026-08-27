import numpy as np

from kmc.particles import Interstitial, Vacancy
from kmc.rates import DissociationFP, FP_generation, Formation, Migration
from kmc.system import System


def make_system(coordinates, materials, *, periodic=False, cell=None, rtol=0.05):
    system = System()
    coordinates = np.asarray(coordinates, dtype=float)
    system.set_morph(
        coordinates[:, 0],
        coordinates[:, 1],
        coordinates[:, 2],
        materials,
        cell=cell,
    )
    system.periodic = periodic
    system.cutoff = 100.0
    system.first_neighbor_rtol = rtol
    system.build_neighbor_topology()
    system.s1 = np.zeros(len(coordinates))
    return system


def test_migration_uses_material_resolved_first_shell_only():
    system = make_system(
        [[0, 0, 0], [5, 0, 0], [2, 0, 0], [7, 0, 0]],
        [0, 0, 1, 1],
        periodic=True,
        cell=np.diag([10.0, 10.0, 10.0]),
    )
    particle = Interstitial(0)
    system.set_particles([particle])
    migration = Migration(
        k={(0, 0): 3.0, (0, 1): 0.0, (1, 0): 0.0, (1, 1): 3.0}
    )
    cut, dx, dy, dz, r = system.neighbor_data(0, "first")
    rates = migration.rate(
        r=r,
        dx=dx,
        dy=dy,
        dz=dz,
        system=system,
        particle=particle,
        mats=system.mats[cut],
        matlocal=system.mats[0],
        cut=cut,
    )
    assert set(cut[rates > 0]) == {1}


def test_formation_restores_the_original_ghost_material():
    system = make_system([[0, 0, 0], [1, 0, 0]], [0, 1])
    interstitial = Interstitial(0)
    vacancy = Vacancy(1)
    system.set_particles([interstitial, vacancy])
    formation = Formation(
        k={(0, 0): 0.0, (0, 1): 2.0, (1, 0): 2.0, (1, 1): 0.0}
    )

    cut, dx, dy, dz, r = system.neighbor_data(0, "first")
    rates = formation.rate(
        r=r,
        dx=dx,
        dy=dy,
        dz=dz,
        system=system,
        particle=interstitial,
        mats=system.mats[cut],
        matlocal=system.mats[0],
        cut=cut,
    )
    assert list(cut[rates > 0]) == [1]

    formation.action(interstitial, system, 1)
    pair = next(p for p in system.particles if p.species == "frenkelpair0")
    assert pair.origin_site == 1
    assert system.mats[1] == 999

    DissociationFP(k={0: 1.0, 1: 1.0}).action(pair, system, pair.position)
    assert system.mats[1] == 1
    assert sorted(p.species for p in system.particles) == ["interstitial", "vacancy"]


def test_vacancy_initiated_formation_uses_the_vacancy_as_ghost():
    system = make_system([[0, 0, 0], [1, 0, 0]], [0, 1])
    interstitial = Interstitial(0)
    vacancy = Vacancy(1)
    system.set_particles([interstitial, vacancy])
    formation = Formation(
        k={(0, 0): 0.0, (0, 1): 2.0, (1, 0): 2.0, (1, 1): 0.0}
    )
    formation.action(vacancy, system, 0)
    pair = next(p for p in system.particles if p.species == "frenkelpair0")
    assert pair.position == 0
    assert pair.ghost_site == 1
    assert pair.origin_site == 1


def test_generation_uses_first_neighbors_and_disables_when_exhausted():
    system = make_system(
        [[0, 0, 0], [1, 0, 0], [4, 0, 0]],
        [0, 1, 1],
    )
    generator = FP_generation(k=5.0, pairs=[[1, 0]])
    assert generator.available_pairs(system) == [(0, 1)]
    assert generator.rate(system=system) == 5.0

    generator.action(system, num=1)
    assert generator.rate(system=system) == 0.0


def test_triclinic_minimum_image_uses_full_cell_matrix():
    cell = np.asarray([[4.0, 0.0, 0.0], [3.0, 3.0, 0.0], [0.0, 0.0, 5.0]])
    system = make_system(
        [[0, 0, 0], [3.8, 0.1, 0]],
        [0, 0],
        periodic=True,
        cell=cell,
    )
    dx, dy, dz = system.minimum_image_displacement(0, 1)
    assert np.isclose(np.sqrt(dx * dx + dy * dy + dz * dz), np.sqrt(0.05))
