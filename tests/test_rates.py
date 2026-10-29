import numpy as np

from kmc.particles import Interstitial, Singlet
from kmc.rates import Fluor, Forster, Migration
from kmc.system import System


def test_fluorescence_rate_is_inverse_lifetime():
    fluorescence = Fluor(life={0: 4.0, 1: 10.0})

    assert fluorescence.rate(material=0) == 0.25
    assert fluorescence.rate(material=1) == 0.1


def test_forster_rate_matches_implemented_expression():
    transfer = Forster(
        Rf={(0, 0): 6.0, (0, 1): 8.0, (1, 0): 7.0, (1, 1): 9.0},
        life={0: 4.0, 1: 5.0},
        mu={0: 2.0, 1: 3.0},
    )
    distances = np.asarray([0.0, 5.0, 10.0])
    materials = np.asarray([0, 0, 1], dtype=np.int32)

    rates = transfer.rate(
        r=distances,
        system=None,
        particle=Singlet(0),
        mats=materials,
        matlocal=0,
    )

    alpha_mu = 1.15 * 0.53 * 2.0
    expected = np.asarray(
        [
            0.0,
            (6.0 / (alpha_mu + 5.0)) ** 6 / 4.0,
            (8.0 / (alpha_mu + 10.0)) ** 6 / 4.0,
        ]
    )
    np.testing.assert_allclose(rates, expected)


def test_migration_rates_give_exact_cubic_lattice_diffusion_coefficient():
    coordinates = np.asarray(
        [
            [0, 0, 0],
            [1, 0, 0],
            [-1, 0, 0],
            [0, 1, 0],
            [0, -1, 0],
            [0, 0, 1],
            [0, 0, -1],
        ],
        dtype=float,
    )
    system = System()
    system.set_morph(
        coordinates[:, 0],
        coordinates[:, 1],
        coordinates[:, 2],
        np.zeros(len(coordinates), dtype=int),
    )
    system.periodic = False
    system.cutoff = 0.1
    system.first_neighbor_rtol = 0.01
    system.build_neighbor_topology()
    particle = Interstitial(0)
    system.set_particles([particle])
    hop_rate = 2.0e5
    migration = Migration(k={(0, 0): hop_rate})

    cut, dx, dy, dz, distances = system.neighbor_data(0, "first")
    rates = migration.rate(
        r=distances,
        dx=dx,
        dy=dy,
        dz=dz,
        system=system,
        particle=particle,
        mats=system.mats[cut],
        matlocal=0,
        cut=cut,
    )
    diffusion = np.sum(rates * distances**2) / 6.0

    assert len(cut) == 6
    assert diffusion == hop_rate
