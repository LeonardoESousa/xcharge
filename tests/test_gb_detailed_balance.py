"""An axis-dependent GB energy must not drive equilibrium currents."""

import numpy as np

from kmc.potentials import MigrationBarrier
from kmc.system import System


def test_tilted_hops_obey_detailed_balance_and_zero_cycle_affinity():
    system = System()
    system.set_morph([0, 1, 2], [0, 1, 0], [0, 0, 1], [0, 0, 0])
    system.periodic = False
    system.cutoff = 3.0
    system.first_neighbor_rtol = 2.0
    MigrationBarrier(function=lambda y, period: 0.08*y,
                     temperature=298, axis="y", coupling="symmetric",
                     directional=True, angular_exponent=2).assign_to_system(system)
    system.build_neighbor_topology()
    energies = system.migration_barrier_change
    probability = np.exp(-energies/system.migration_kbt)
    probability /= probability.sum()
    rates = {}
    drift = np.zeros(3)
    for origin in range(3):
        destinations, dx, dy, dz, distances = system.neighbor_data(origin, "first")
        factors = system.migration_rate_factors(origin, destinations)
        for destination, factor, x, y, z in zip(destinations, factors, dx, dy, dz):
            rates[origin, int(destination)] = factor
            drift += probability[origin]*factor*np.array([x, y, z])
    for (i, j), rate in rates.items():
        np.testing.assert_allclose(probability[i]*rate, probability[j]*rates[j, i])
    forward = rates[0, 1]*rates[1, 2]*rates[2, 0]
    backward = rates[1, 0]*rates[2, 1]*rates[0, 2]
    np.testing.assert_allclose(forward, backward)
    np.testing.assert_allclose(drift, 0.0, atol=1e-14)
    # Both sites lie in the same GB-parallel plane.
    assert rates[0, 2] == rates[2, 0] == 1.0
