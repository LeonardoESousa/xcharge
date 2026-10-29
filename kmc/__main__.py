import importlib
import inspect
import multiprocessing
import os
import random
import sys

import matplotlib.pyplot as plt
import numpy as np
from matplotlib import animation

try:
    import tqdm
except ImportError:  # pragma: no cover - optional progress display
    tqdm = None

import kmc.bimolecular
import kmc.variables
import kmc.vis
from kmc.events import EventScheduler, ParticleEvents
from kmc.system import System

try:
    from importlib import metadata
except ImportError:  # pragma: no cover - Python < 3.8
    import importlib_metadata as metadata


def _load_parameters(path):
    module_name = os.path.splitext(os.path.basename(path))[0]
    spec = importlib.util.spec_from_file_location(module_name, os.path.abspath(path))
    if spec is None or spec.loader is None:
        raise ValueError(f"Could not load input file: {path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


if len(sys.argv) < 2:
    raise SystemExit("Usage: kmc INPUT.py [input arguments ...]")

param = _load_parameters(sys.argv[1])


def package_version():
    try:
        return metadata.version("KMC")
    except metadata.PackageNotFoundError:
        from kmc.__version__ import __version__

        return __version__


def set_variable(name):
    if hasattr(param, name):
        return getattr(param, name)
    return getattr(kmc.variables, name)


assigners = [
    value
    for value in vars(param).values()
    if hasattr(value, "assign_to_system") and not inspect.isclass(value)
]

monomolecular = param.monomolecular
processes = param.processes
generation = set_variable("generation")
identifier = set_variable("identifier")
animation_mode = set_variable("animation_mode")
save_animation = set_variable("save_animation")
animation_exten = set_variable("animation_exten")
time_limit = set_variable("time_limit")
pause = set_variable("pause")
marker_type = set_variable("marker_type")
rotate = set_variable("rotate")
frozen_lattice = set_variable("frozen_lattice")
bimolec = set_variable("bimolec")
periodic = set_variable("periodic")
n_proc = int(set_variable("n_proc"))
rounds = int(set_variable("rounds"))
cutoff = set_variable("cutoff")
colors_dic = set_variable("colors_dic")
sizes_dic = set_variable("sizes_dic")
square_ratio = set_variable("square_ratio")
scatter_alpha = set_variable("scatter_alpha")
clean_vis = set_variable("clean_vis")
material_label = set_variable("material_label")
material_leg = set_variable("material_leg")
time_units = set_variable("time_units")
particle_condition = set_variable("particle_condition")
creation_threshold = set_variable("creation_threshold")
print_site_position = set_variable("print_site_position")
first_neighbor_rtol = set_variable("first_neighbor_rtol")
random_seed = set_variable("random_seed")
append_output = set_variable("append_output")

PARTICLE_ASSIGNERS = {
    "Create_Particles",
    "Create_Particles_PROB",
    "Create_ParticlesFP",
}


def regular_distance(system, local, destination=None):
    if destination is None:
        displacement = system.R - system.R[local]
    else:
        displacement = system.R[destination] - system.R[local]
    return tuple(np.moveaxis(displacement, -1, 0))


def periodic_distance(system, local, destination=None):
    return system.minimum_image_displacement(local, destination)


distance = periodic_distance if periodic else regular_distance


def _seed_legacy_rngs(seed):
    if seed is None:
        return
    seed = int(seed)
    random.seed(seed)
    np.random.seed(seed % (2**32))


def make_system(seed=None, include_particles=True):
    _seed_legacy_rngs(seed)
    system = System()
    system.set_rng(seed)
    system.set_basic_info(
        monomolecular=monomolecular,
        processes=processes,
        identifier=identifier,
        animation_mode=animation_mode,
        time_limit=time_limit,
        pause=pause,
        anni=bimolec,
        distance=distance,
        generation=generation,
        cutoff=cutoff,
        periodic=periodic,
        first_neighbor_rtol=first_neighbor_rtol,
        particle_condition=particle_condition,
        creation_threshold=creation_threshold,
        time_units=time_units,
    )

    for assigner in assigners:
        if (
            assigner.__class__.__name__ not in PARTICLE_ASSIGNERS
            and not getattr(assigner, "requires_morphology", False)
        ):
            assigner.assign_to_system(system)
    if not hasattr(system, "R"):
        raise ValueError("No morphology was assigned to the system.")

    for assigner in assigners:
        if getattr(assigner, "requires_morphology", False):
            assigner.assign_to_system(system)

    system.build_neighbor_topology()

    if include_particles:
        for assigner in assigners:
            if assigner.__class__.__name__ in PARTICLE_ASSIGNERS:
                assigner.assign_to_system(system)
        system.event_scheduler = EventScheduler(system, _build_particle_events)
        system._clear_event_changes()
    return system


def anni_general(system, reactions, local):
    particles = list(system.particles_at_site(local))
    if len(particles) < 2:
        return
    pair = [0, 1]
    key = tuple(sorted(particles[index].species for index in pair))
    reaction = reactions.get(key)
    if reaction is not None:
        reaction(particles, system, pair)


def _checked_rates(values, process):
    values = np.asarray(values, dtype=float).reshape(-1)
    if np.any(~np.isfinite(values)):
        raise ValueError(f"Non-finite rate produced by {process!r}.")
    if np.any(values < 0):
        raise ValueError(f"Negative rate produced by {process!r}.")
    return values


def _build_particle_events(particle, system):
    local = particle.position
    rate_blocks = []
    destinations = []
    labels = []

    for transfer in system.processes[particle.species]:
        mode = getattr(transfer, "neighbor_mode", "cutoff")
        cut, dx, dy, dz, r = system.neighbor_data(local, mode)
        mats = system.mats[cut]
        rates = transfer.rate(
            r=r,
            dx=dx,
            dy=dy,
            dz=dz,
            system=system,
            particle=particle,
            mats=mats,
            matlocal=system.mats[local],
            cut=cut,
        )
        rates = _checked_rates(rates, transfer)
        if len(rates) != len(cut):
            raise ValueError(
                f"{transfer!r} returned {len(rates)} rates for {len(cut)} destinations."
            )
        positive = rates > 0.0
        if np.any(positive):
            selected_rates = rates[positive]
            rate_blocks.append(selected_rates)
            destinations.append(np.asarray(cut[positive], dtype=np.int32))
            labels.extend([transfer] * len(selected_rates))

    for process in system.monomolecular[particle.species]:
        rate = _checked_rates(
            [process.rate(material=system.mats[local])], process
        )
        if rate[0] > 0.0:
            rate_blocks.append(rate)
            destinations.append(np.asarray([local], dtype=np.int32))
            labels.append(process)

    if not rate_blocks:
        return ParticleEvents.empty()

    rates = np.concatenate(rate_blocks)
    return ParticleEvents(
        cumulative_rates=np.cumsum(rates),
        processes=tuple(labels),
        destinations=np.concatenate(destinations),
        total=float(rates.sum()),
    )


def decision(particle, system):
    """Compatibility wrapper returning one particle's total propensity."""
    catalogue = _build_particle_events(particle, system)
    if catalogue.total <= 0.0:
        particle.process = None
        particle.destination = None
        return 0.0
    particle.process, particle.destination = catalogue.select(system.rng)
    return catalogue.total


def choose_generation(system):
    rates = _checked_rates(
        [generator.rate(system=system) for generator in system.generation],
        "generation",
    )
    total = float(rates.sum())
    if total <= 0:
        return None, 0.0
    index = int(
        np.searchsorted(np.cumsum(rates), system.rng.random() * total, side="right")
    )
    return system.generation[index], total


def _active(system):
    if system.time >= system.time_limit:
        return False
    if system.particle_condition and system.count_particles() == 0:
        return False
    return True


def _advance_to_next_event(system):
    """Execute one exact KMC event, respecting deterministic time boundaries."""
    while _active(system):
        generation_enabled = system.time < system.creation_threshold
        if generation_enabled:
            generation_event, generation_rate = choose_generation(system)
        else:
            generation_event, generation_rate = None, 0.0

        if system.event_scheduler is None:
            system.event_scheduler = EventScheduler(system, _build_particle_events)
            system._clear_event_changes()
        particle_total = system.event_scheduler.total_rate
        total_rate = generation_rate + particle_total

        boundaries = [system.time_limit]
        if generation_enabled and np.isfinite(system.creation_threshold):
            boundaries.append(system.creation_threshold)
        boundary = min(value for value in boundaries if value >= system.time)

        if total_rate <= 0:
            if np.isfinite(boundary) and boundary > system.time:
                system.time = boundary
                continue
            return False

        dt = -np.log(max(system.rng.random(), np.finfo(float).tiny)) / total_rate
        if system.time + dt >= boundary:
            system.time = boundary
            continue

        system.time += dt
        system.IT += 1
        target = system.rng.random() * total_rate
        if target < generation_rate:
            system.begin_event()
            generation_event.action(system, num=1)
            system.event_scheduler.apply_changes(system.consume_event_changes())
            return True

        target -= generation_rate
        # Guard against the last-bit roundoff that can otherwise put a target
        # exactly on the open upper edge of the Fenwick tree.
        target = min(target, np.nextafter(particle_total, -np.inf))
        particle, process, destination = system.event_scheduler.select(target)
        particle.process = process
        particle.destination = destination
        system.begin_event(particle, destination)
        process.action(particle, system, destination)
        if bimolec:
            anni_general(
                system,
                kmc.bimolecular.bimolec_funcs_array,
                destination,
            )
        particle.stamp_time(system)
        system.event_scheduler.apply_changes(
            system.consume_event_changes()
        )
        return True
    return False


def _energy_array(system, particle):
    if particle.species == "electron" and hasattr(system, "lumo"):
        return system.lumo
    if particle.species == "hole" and hasattr(system, "homo"):
        return system.homo
    if particle.species == "triplet" and hasattr(system, "t1"):
        return system.t1
    return system.s1


def _finalize(system):
    for particle in system.particles.copy():
        particle.kill("alive", system, _energy_array(system, particle), "alive")
        particle.stamp_time(system)


def step_nonani(system):
    while _advance_to_next_event(system):
        pass
    _finalize(system)


def step_ani(system):
    if _advance_to_next_event(system):
        return system.particles.copy()
    return []


step = step_ani if animation_mode else step_nonani


def open_log():
    filename = f"Simulation_{identifier}.txt"
    if not os.path.isfile(filename) or not append_output:
        version = package_version()
        with open(filename, "w", encoding="utf-8") as output:
            output.write(f"# Version: {version}\n")
            output.write(f"# Time units: {time_units}\n")
            output.write(
                "Time,DeltaX,DeltaY,DeltaZ,Type,Energy,Location,"
                "FinalX,FinalY,FinalZ,CausaMortis,Status\n"
            )
    return filename


def spectra_text(system):
    return "".join(particle.write() for particle in system.dead) + "END\n"


def animate(_frame, viewer):
    return viewer.update()


def run_animation():
    system = make_system(seed=random_seed)
    figure = plt.figure()
    axis = figure.add_subplot(111, projection="3d")
    viewer = kmc.vis.KMCViewer(
        system=system,
        step_function=step,
        fig=figure,
        ax=axis,
        colors_dic=colors_dic,
        sizes_dic=sizes_dic,
        material_label=material_label,
        marker_option=marker_type,
        rotate=rotate,
        square_ratio=square_ratio,
        clean_vis=clean_vis,
        scatter_alpha=scatter_alpha,
        print_site_position=print_site_position,
        material_leg=material_leg,
        time_units=time_units,
    )
    ani = animation.FuncAnimation(
        figure,
        animate,
        fargs=[viewer],
        interval=1,
        blit=False,
        repeat=False,
        cache_frame_data=False,
        save_count=1000,
    )
    return ani


_frozen_template = None


def _init_worker(template):
    global _frozen_template
    _frozen_template = template


def _particle_assigners(system):
    for assigner in assigners:
        if assigner.__class__.__name__ in PARTICLE_ASSIGNERS:
            assigner.assign_to_system(system)


def reroll_system(system, seed):
    _seed_legacy_rngs(seed)
    system.reset_particles()
    system.set_rng(seed)
    _particle_assigners(system)
    system.event_scheduler = EventScheduler(system, _build_particle_events)
    system._clear_event_changes()
    return system


def run_dynamic(task):
    _, seed = task
    system = make_system(seed=seed)
    step_nonani(system)
    return spectra_text(system)


def run_frozen(task):
    _, seed = task
    if _frozen_template is None:
        raise RuntimeError("Frozen-lattice worker was not initialized.")
    system = reroll_system(_frozen_template, seed)
    step_nonani(system)
    return spectra_text(system)


def _round_seeds():
    sequence = np.random.SeedSequence(random_seed)
    return [int(child.generate_state(1)[0]) for child in sequence.spawn(rounds)]


def main():
    version = package_version()
    print("####################################################################")
    print("Xcharge: A Kinetic Monte Carlo Model for Exciton and Charge Dynamics")
    print(f"Version  : {version}\n")

    if animation_mode:
        ani = run_animation()
        path = identifier + "_animation." + animation_exten
        if save_animation and animation_exten == "gif":
            ani.save(path, writer="pillow", fps=10)
        elif save_animation and animation_exten == "mp4":
            writer = animation.FFMpegWriter(
                fps=10,
                extra_args=["-crf", "10", "-preset", "slow", "-pix_fmt", "yuv420p"],
            )
            ani.save(path, writer=writer)
        plt.show()
        return 0

    filename = open_log()
    seeds = _round_seeds()
    tasks = list(enumerate(seeds))
    template = None
    runner = run_dynamic
    if frozen_lattice:
        template_seed = int(np.random.SeedSequence(random_seed).generate_state(1)[0])
        template = make_system(seed=template_seed, include_particles=False)
        runner = run_frozen

    with open(filename, "a", encoding="utf-8") as output:
        if n_proc == 1:
            _init_worker(template)
            results = map(runner, tasks)
            iterator = tqdm.tqdm(results, total=rounds) if tqdm else results
            for result in iterator:
                output.write(result)
        else:
            chunksize = max(1, rounds // (4 * n_proc))
            with multiprocessing.Pool(
                n_proc,
                initializer=_init_worker,
                initargs=(template,),
            ) as pool:
                results = pool.imap(runner, tasks, chunksize=chunksize)
                iterator = tqdm.tqdm(results, total=rounds) if tqdm else results
                for result in iterator:
                    output.write(result)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
