import numpy as np

from kmc.particles import Interstitial, Vacancy
from kmc.potentials import BOLTZMANN_EV, MigrationBarrier
from kmc.rates import Migration
from kmc.system import System


def make_system():
    system = System()
    system.set_morph(
        [0.0, 1.0, 2.0],
        [0.0, 0.0, 0.0],
        [0.0, 0.0, 0.0],
        [0, 0, 0],
    )
    return system


def migration_rates(system, particle):
    migration = Migration(k={(0, 0): 10.0})
    cut = np.asarray([1, 2], dtype=np.int32)
    return migration.rate(
        system=system,
        particle=particle,
        mats=system.mats[cut],
        matlocal=0,
        cut=cut,
    )


def test_migration_is_unchanged_without_barrier_profile():
    system = make_system()
    particle = Interstitial(0)
    system.set_particles([particle])

    np.testing.assert_allclose(migration_rates(system, particle), [10.0, 10.0])


def test_same_barrier_profile_biases_interstitial_and_vacancy_migration():
    temperature = 300.0
    system = make_system()
    system.migration_barrier_change = np.asarray([0.0, -0.05, 0.05])
    system.migration_kbt = BOLTZMANN_EV * temperature

    expected = 10.0 * np.exp(
        -system.migration_barrier_change[[1, 2]] / system.migration_kbt
    )

    for particle_class in (Interstitial, Vacancy):
        particle = particle_class(0)
        system.set_particles([particle])
        np.testing.assert_allclose(migration_rates(system, particle), expected)
        system.remove(particle)


def test_periodic_profile_is_evaluated_once_at_wrapped_site_coordinates():
    system = System()
    system.set_morph(
        [0.0, 0.5, 9.5],
        [0.0, 0.0, 0.0],
        [0.0, 0.0, 0.0],
        [0, 0, 0],
        cell=np.diag([10.0, 1.0, 1.0]),
    )
    system.periodic = True

    def boundary_profile(x, period):
        distance = np.minimum(x, period - x)
        return -0.1 * np.exp(-(distance / 1.0) ** 2)

    MigrationBarrier(
        boundary_profile,
        temperature=300.0,
        axis="x",
    ).assign_to_system(system)

    assert system.migration_barrier_change[0] == -0.1
    assert np.isclose(
        system.migration_barrier_change[1],
        system.migration_barrier_change[2],
    )


def test_directional_profile_uses_smooth_angular_weight():
    temperature = 300.0
    system = System()
    system.set_morph(
        [0.0, 0.0, 1.0, 1.0],
        [0.0, 1.0, 0.0, 1.0],
        [0.0, 0.0, 0.0, 0.0],
        [0, 0, 0, 0],
        cell=np.diag([10.0, 10.0, 10.0]),
    )
    system.periodic = True
    system.migration_barrier_change = np.asarray([0.0, 0.06, 0.06, 0.06])
    system.migration_kbt = BOLTZMANN_EV * temperature
    system.migration_barrier_axis = 1
    system.migration_barrier_directional = True
    system.migration_barrier_angular_exponent = 2.0

    particle = Interstitial(0)
    system.set_particles([particle])
    migration = Migration(k={(0, 0): 10.0})
    cut = np.asarray([1, 2, 3], dtype=np.int32)
    actual = migration.rate(
        system=system,
        particle=particle,
        mats=system.mats[cut],
        matlocal=0,
        cut=cut,
    )

    expected_barriers = np.asarray([0.06, 0.0, 0.03])
    expected = 10.0 * np.exp(
        -expected_barriers / system.migration_kbt
    )
    np.testing.assert_allclose(actual, expected)


def test_nondirectional_profile_preserves_previous_behavior():
    temperature = 300.0
    system = make_system()
    system.migration_barrier_change = np.asarray([0.0, 0.05, 0.05])
    system.migration_kbt = BOLTZMANN_EV * temperature
    system.migration_barrier_axis = 1
    system.migration_barrier_directional = False

    particle = Vacancy(0)
    system.set_particles([particle])

    expected = 10.0 * np.exp(-0.05 / system.migration_kbt)
    np.testing.assert_allclose(
        migration_rates(system, particle),
        [expected, expected],
    )


def test_symmetric_coupling_uses_site_energy_difference():
    temperature = 300.0
    system = make_system()
    system.migration_barrier_change = np.asarray([0.04, 0.02, 0.06])
    system.migration_kbt = BOLTZMANN_EV * temperature
    system.migration_barrier_coupling = "symmetric"

    particle = Interstitial(0)
    system.set_particles([particle])
    actual = migration_rates(system, particle)

    expected_delta_barriers = 0.5 * np.asarray([-0.02, 0.02])
    expected = 10.0 * np.exp(
        -expected_delta_barriers / system.migration_kbt
    )
    np.testing.assert_allclose(actual, expected)
