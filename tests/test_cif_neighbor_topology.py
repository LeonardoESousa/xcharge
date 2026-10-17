"""Regression checks for geometry shared by bulk and GB material labels."""

import numpy as np

from kmc.morphology import ReadCIF
from kmc.particles import Interstitial
from kmc.rates import Migration
from kmc.system import System


def cif_system(tmp_path, *, labelled=True):
    rows = []
    for element, offset in (("Cs", 0), ("Br", 1)):
        for number in range(6):
            suffix = "_GB" if labelled and number == 5 else ""
            rows.append(
                f"{element} {element}{number + 1}{suffix} "
                f"{(2 * number + offset) / 12:.12f} 0 0"
            )
    path = tmp_path / "chain.cif"
    path.write_text(
        "data_chain\n_cell_length_a 12\n_cell_length_b 20\n"
        "_cell_length_c 20\n_cell_angle_alpha 90\n_cell_angle_beta 90\n"
        "_cell_angle_gamma 90\n_space_group_name_H-M_alt 'P 1'\n"
        "loop_\n_space_group_symop_operation_xyz\n'x,y,z'\n"
        "loop_\n_atom_site_type_symbol\n_atom_site_label\n"
        "_atom_site_fract_x\n_atom_site_fract_y\n_atom_site_fract_z\n"
        + "\n".join(rows) + "\n"
    )
    system = System()
    ReadCIF(str(path), [1, 1, 1],
            {"Cs": 0, "Br": 1, "Cs_GB": 2, "Br_GB": 3}).assign_to_system(system)
    system.periodic = True
    system.cutoff = 0.1
    system.first_neighbor_rtol = 0.05
    system.build_neighbor_topology()
    system.distance = lambda s, i, j: s.minimum_image_displacement(i, j)
    return system


def edges(system):
    return {
        (i, int(j)): np.array([x, y, z])
        for i, (destinations, dx, dy, dz, _) in enumerate(system._first_neighbors)
        for j, x, y, z in zip(destinations, dx, dy, dz)
    }


def test_gb_labels_do_not_change_physical_neighbors(tmp_path):
    labelled = cif_system(tmp_path)
    bulk = cif_system(tmp_path, labelled=False)
    assert set(labelled.mats) == {0, 1, 2, 3}
    assert set(labelled.neighbor_species) == {"Cs", "Br"}
    labelled_edges, bulk_edges = edges(labelled), edges(bulk)
    assert labelled_edges.keys() == bulk_edges.keys()
    for edge, vector in labelled_edges.items():
        np.testing.assert_allclose(vector, bulk_edges[edge], atol=1e-9)
        i, j = edge
        if labelled.neighbor_species[i] == labelled.neighbor_species[j]:
            assert np.isclose(np.linalg.norm(vector), 2.0)
            np.testing.assert_allclose(vector, -labelled_edges[j, i], atol=1e-9)


def test_gb_rates_and_occupancy_still_apply(tmp_path):
    system = cif_system(tmp_path)
    origin = int(np.flatnonzero(np.isclose(system.X, 9))[0])
    particle = Interstitial(origin)
    system.set_particles([particle])
    migration = Migration(k={(a, b): (7.0 if a >= 2 or b >= 2 else 3.0)
                             for a in (1, 3) for b in (1, 3)})
    cut, dx, dy, dz, r = system.neighbor_data(origin, "first")
    def rates():
        return migration.rate(r=r, dx=dx, dy=dy, dz=dz, system=system,
                              particle=particle, mats=system.mats[cut],
                              matlocal=system.mats[origin], cut=cut)
    values = rates()
    assert sorted(values[values > 0]) == [3.0, 7.0]
    gb_site = int(cut[values == 7.0][0])
    system.set_particles([Interstitial(gb_site)])
    assert rates()[cut == gb_site][0] == 0.0


def test_periodic_hop_and_reverse_cancel(tmp_path):
    system = cif_system(tmp_path)
    i = int(np.flatnonzero(np.isclose(system.X, 1))[0])
    j = int(np.flatnonzero(np.isclose(system.X, 11))[0])
    particle = Interstitial(i)
    system.set_particles([particle])
    particle.move(j, system)
    np.testing.assert_allclose([particle.Dx, particle.Dy, particle.Dz],
                               [-2, 0, 0], atol=1e-9)
    particle.move(i, system)
    np.testing.assert_allclose([particle.Dx, particle.Dy, particle.Dz],
                               [0, 0, 0], atol=1e-9)


def test_distorted_sublattice_has_reverse_hops():
    system = System()
    system.set_morph([0, 1, 2.2], [0]*3, [0]*3, [0]*3)
    system.periodic = False
    system.cutoff = 0.1
    system.first_neighbor_rtol = 0.05
    system.build_neighbor_topology()
    graph = edges(system)
    assert set(graph) == {(0, 1), (1, 0), (1, 2), (2, 1)}
    for (i, j), vector in graph.items():
        np.testing.assert_allclose(vector, -graph[j, i])


def test_homogeneous_periodic_walk_has_no_drift(tmp_path):
    system = cif_system(tmp_path)
    graph = edges(system)
    starts = np.flatnonzero(system.mats % 2 == 1)
    neighbors = {}
    for i in starts:
        neighbors[i] = [(j, v) for (origin, j), v in graph.items()
                        if origin == i and system.mats[j] % 2 == 1]
        assert len(neighbors[i]) == 2
    rng = np.random.default_rng(294)
    displacements = []
    for _ in range(1000):
        position = int(rng.choice(starts))
        delta = np.zeros(3)
        for _ in range(200):
            position, vector = neighbors[position][rng.integers(2)]
            delta += vector
        displacements.append(delta[0])
    displacements = np.asarray(displacements)
    expected_msd = 200 * 2.0**2
    assert abs(displacements.mean()) < 4*np.sqrt(expected_msd/len(displacements))
    assert np.isclose(np.mean(displacements**2), expected_msd, rtol=0.15)
