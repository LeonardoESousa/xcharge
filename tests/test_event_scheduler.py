import numpy as np

from kmc.events import EventScheduler, FenwickTree, ParticleEvents
from kmc.particles import Vacancy
from kmc.system import System


def _catalogue(total):
    return ParticleEvents(
        cumulative_rates=np.asarray([total], dtype=float),
        processes=("move",),
        destinations=np.asarray([0], dtype=np.int32),
        total=float(total),
    )


def test_fenwick_tree_updates_and_selects_exact_intervals():
    tree = FenwickTree(4)
    for index, value in enumerate((2.0, 0.0, 3.0, 5.0)):
        tree.set(index, value)

    assert tree.total == 10.0
    assert [tree.find(value) for value in (0.0, 1.999, 2.0, 4.999, 5.0, 9.999)] == [
        0,
        0,
        2,
        2,
        3,
        3,
    ]

    tree.grow(9)
    tree.set(2, 1.0)
    tree.set(8, 4.0)
    assert tree.total == 12.0
    assert tree.find(11.999) == 8


def test_system_particle_indexes_follow_add_move_and_remove():
    system = System()
    system.set_morph([0.0, 1.0, 2.0], [0.0] * 3, [0.0] * 3, [0, 0, 0])
    first = Vacancy(0)
    second = Vacancy(2)
    system.set_particles([first, second])
    system._clear_event_changes()

    first.move(1, system)
    assert system.site_occupancy.tolist() == [0, 1, 1]
    assert system.species_occupancy("vacancy", [0, 1, 2]).tolist() == [
        False,
        True,
        True,
    ]
    assert system.particles_at_site(1) == (first,)

    system.remove(first)
    assert system.particles == [second]
    assert system.particles_at_site(1) == ()
    assert system.site_occupancy.tolist() == [0, 0, 1]


def test_scheduler_refreshes_only_particles_affected_by_local_change():
    system = System()
    system.set_morph([0.0, 1.0, 5.0], [0.0] * 3, [0.0] * 3, [0, 0, 0])
    moving = Vacancy(0)
    nearby = Vacancy(1)
    distant = Vacancy(2)
    system.set_particles([moving, nearby, distant])
    system._reverse_neighbors = [{0, 1}, {0, 1}, {2}]
    system._clear_event_changes()

    evaluations = {moving: 0, nearby: 0, distant: 0}

    def build(particle, _system):
        evaluations[particle] += 1
        return _catalogue(particle.position + 1.0)

    scheduler = EventScheduler(system, build)
    evaluations.update({moving: 0, nearby: 0, distant: 0})

    system.begin_event(moving, 1)
    moving.move(1, system)
    scheduler.apply_changes(system.consume_event_changes())

    assert evaluations[moving] == 1
    assert evaluations[nearby] == 1
    assert evaluations[distant] == 0
    assert scheduler.total_rate == 2.0 + 2.0 + 3.0


def test_scheduler_registers_and_unregisters_particles():
    system = System()
    system.set_morph([0.0, 1.0], [0.0, 0.0], [0.0, 0.0], [0, 0])
    original = Vacancy(0)
    system.set_particles([original])
    system._clear_event_changes()
    scheduler = EventScheduler(system, lambda particle, _: _catalogue(2.0))

    system.begin_event(original, 1)
    replacement = Vacancy(1)
    system.set_particles([replacement])
    system.remove(original)
    scheduler.apply_changes(system.consume_event_changes())

    assert original not in scheduler.particle_slots
    assert replacement in scheduler.particle_slots
    assert scheduler.total_rate == 2.0
