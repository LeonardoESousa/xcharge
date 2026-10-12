"""Independent-particle, spatial slow-scale averaging of defect diffusion.

Migration is a virtual fast process. Its stationary distributions weight the
original slow channels, including their origins and destinations. This is an
approximation requiring fast spatial mixing and negligible free-defect
correlations/blocking; it is not a first-encounter solver.
"""

import copy
from dataclasses import dataclass

import numpy as np
from scipy.sparse import coo_matrix, diags
from scipy.sparse.csgraph import connected_components
from scipy.sparse.linalg import spsolve, splu, cg

from kmc.rates import FP_generation, Formation, Migration


FAST_SPECIES = frozenset(("interstitial", "vacancy"))


def _rates(values):
    values = np.asarray(values, dtype=float).reshape(-1)
    if np.any(~np.isfinite(values)) or np.any(values < 0):
        raise ValueError("Diffusion averaging requires finite, non-negative rates.")
    return values


class _FastView:
    """Read-only geometry with only retained particles blocking diffusion."""

    def __init__(self, system):
        self.system = system
        self.site_occupancy = np.zeros(len(system.mats), dtype=np.int32)
        for particle in system.particles:
            if particle.species not in FAST_SPECIES:
                self.site_occupancy[particle.position] += 1
        self.site_occupancy[system.mats == 999] = 1
        self.occupancy = {
            int(i): int(self.site_occupancy[i])
            for i in np.flatnonzero(self.site_occupancy)
        }

    def __getattr__(self, name):
        return getattr(self.system, name)

    def species_occupancy(self, _species, sites):
        # Formation.rate supplies the intrinsic spatial kernel; actual partner
        # probabilities are multiplied in separately, rather than counted twice.
        return self.site_occupancy[np.asarray(sites, dtype=int)] == 0


@dataclass
class MixedGroup:
    species: str
    particles: tuple
    probability: np.ndarray
    diffusion: np.ndarray
    drift: np.ndarray


@dataclass
class AveragedChannel:
    process: object
    initiators: tuple
    counterparts: tuple
    origins: np.ndarray
    destinations: np.ndarray
    cumulative_rates: np.ndarray
    total: float
    generation: bool = False

    def select(self, rng):
        target = min(rng.random() * self.total, np.nextafter(self.total, -np.inf))
        index = int(np.searchsorted(self.cumulative_rates, target, side="right"))
        return int(self.origins[index]), int(self.destinations[index])


@dataclass
class AveragedModel:
    groups: tuple
    channels: tuple
    cumulative_rates: np.ndarray
    total_rate: float
    eligible_pairs: float


def stationary_distribution(generator, tolerance=1e-10):
    """Solve pi Q=0 on one irreducible class, with a sparse residual check.

    A spanning-tree detailed-balance construction avoids LU fill-in for the
    usual reversible migration graphs. Nonreversible graphs use sparse LU.
    """
    n = generator.shape[0]
    if n == 1:
        return np.ones(1), 0.0
    scale = float(np.max(-generator.diagonal()))
    if not np.isfinite(scale) or scale <= 0:
        raise ValueError("Invalid diffusion generator.")
    q = (generator / scale).tocsr()
    q.sort_indices()
    log_weight = np.full(n, np.nan)
    log_weight[0] = 0.0
    stack = [0]
    reversible = True
    while stack:
        i = stack.pop()
        for offset in range(q.indptr[i], q.indptr[i + 1]):
            j, forward = int(q.indices[offset]), float(q.data[offset])
            if i == j or forward <= 0 or not np.isnan(log_weight[j]):
                continue
            start, stop = q.indptr[j], q.indptr[j + 1]
            reverse_offset = int(np.searchsorted(q.indices[start:stop], i)) + start
            if reverse_offset >= stop or q.indices[reverse_offset] != i:
                reversible = False
                break
            reverse = float(q.data[reverse_offset])
            if reverse <= 0:
                reversible = False
                break
            log_weight[j] = log_weight[i] + np.log(forward) - np.log(reverse)
            stack.append(j)
        if not reversible:
            break
    if reversible and np.all(np.isfinite(log_weight)):
        pi = np.exp(log_weight - log_weight.max())
        pi /= pi.sum()
        residual = float(np.linalg.norm(q.T @ pi, ord=1))
        if residual <= tolerance:
            return pi, residual
    matrix = q.T.tolil()
    matrix[-1, :] = np.ones(n)
    rhs = np.zeros(n)
    rhs[-1] = 1.0
    pi = np.asarray(spsolve(matrix.tocsc(), rhs))
    if np.any(~np.isfinite(pi)) or np.min(pi) < -tolerance:
        raise ValueError("Could not solve the stationary diffusion distribution.")
    pi = np.maximum(pi, 0.0)
    if pi.sum() <= 0:
        raise ValueError("Empty stationary diffusion distribution.")
    pi /= pi.sum()
    residual = float(np.linalg.norm(q.T @ pi, ord=1))
    if residual > tolerance:
        raise ValueError(f"Stationary diffusion residual is too large: {residual:g}.")
    return pi, residual


def effective_transport(graph, probability, displacements, tolerance=1e-9, channels=None):
    """Long-time centered diffusion tensor and drift of a lattice walk.

    Solve Q chi = v-b with one fixed gauge, v = sum(pi_i b_i).
    Corrected edge increments are
    d_ij + chi_j - chi_i; D = sum(pi_i w_ij increment increment.T)/2.
    Bounded directions cancel through chi rather than spuriously diffusing.
    Without channels, edge vectors follow graph.tocoo() order. With channels,
    unmerged (row, column, rate) arrays retain distinct periodic-image jumps.
    This describes virtual independent migration, not reaction-conditioned paths.
    """
    graph = graph.astype(float)
    edges = (graph.tocoo() if channels is None else
             coo_matrix((channels[2], (channels[0], channels[1])), shape=graph.shape))
    n = graph.shape[0]
    vectors = np.asarray(displacements, dtype=float).reshape(-1, 3)
    pi = np.asarray(probability, dtype=float)
    if vectors.shape != (len(edges.data), 3) or np.any(~np.isfinite(vectors)):
        raise ValueError("Invalid diffusion edge displacements.")
    if not len(edges.data):
        return np.zeros((3, 3)), np.zeros(3)
    exit_rates = np.asarray(graph.sum(axis=1)).ravel()
    scale = float(exit_rates.max())
    q = (graph - diags(exit_rates)) / scale
    weights = edges.data / scale
    b = np.zeros((n, 3))
    np.add.at(b, edges.row, weights[:, None] * vectors)
    velocity = pi @ b
    characteristic = float(np.sum(pi[edges.row] * weights * np.linalg.norm(vectors, axis=1)))
    centered_b = b - velocity
    chi = np.zeros((n, 3))
    # The unweighted periodic cell can be large. Reversible generators permit
    # a symmetric positive definite grounded solve, avoiding sparse LU fill-in.
    laplacian = (-diags(pi) @ q).tocsr()
    asymmetry = laplacian - laplacian.T
    reversible = (not asymmetry.nnz or np.max(np.abs(asymmetry.data)) <=
                  1e-10 * max(float(np.max(np.abs(laplacian.data))), np.finfo(float).tiny))
    if n > 5000 and reversible:
        matrix = laplacian[1:, 1:].tocsr()
        preconditioner = diags(1.0 / matrix.diagonal())
        for axis in range(3):
            rhs = pi[1:] * centered_b[1:, axis]
            if not np.any(rhs):
                continue
            solution, info = cg(matrix, rhs, M=preconditioner,
                                rtol=1e-11, atol=0.0, maxiter=max(1000, 5*n))
            if info != 0:
                raise ValueError("Diffusion corrector did not converge.")
            chi[1:, axis] = solution
    elif n > 1:
        chi[1:] = splu((-q[1:, 1:]).tocsc()).solve(centered_b[1:])
    residual = np.max(np.abs(q @ chi + centered_b))
    reference = max(float(np.max(np.abs(centered_b))), characteristic, np.finfo(float).tiny)
    if np.any(~np.isfinite(chi)) or residual > 10*tolerance*reference:
        raise ValueError("Diffusion corrector residual is too large.")
    corrected = vectors + chi[edges.col] - chi[edges.row]
    # Suppress roundoff in directions whose displacement is a bounded state
    # function; this is essential when physical hopping rates are enormous.
    corrected[np.abs(corrected) < 1e-10 * max(float(np.max(np.abs(vectors))), 1e-300)] = 0
    tensor = 0.5 * np.einsum('e,ei,ej->ij', pi[edges.row]*edges.data,
                             corrected, corrected)
    if np.any(~np.isfinite(tensor)):
        raise ValueError("Diffusion tensor overflowed.")
    drift = velocity * scale
    if np.any(~np.isfinite(drift)):
        raise ValueError("Transport drift overflowed.")
    return tensor, drift


def effective_diffusion(graph, probability, displacements, tolerance=1e-9):
    """Backward-compatible tensor-only interface; tensor is centered on drift."""
    return effective_transport(graph, probability, displacements, tolerance)[0]


class DiffusionAverager:
    def __init__(self, system, event_builder):
        if system.animation_mode or system.bimolec:
            raise ValueError("diffusion_average requires animation_mode=False and bimolec=False.")
        self.system = system
        self.event_builder = event_builder
        self.model = None
        self.signature = None
        self.stats = dict(stationary_solves=0, model_builds=0, averaged_events=0,
                          boundary_steps=0, max_stationary_residual=0.0,
                          displacement_available=False)

    def _signature(self):
        s = self.system
        return (s.mats.tobytes(), tuple(
            (id(p), p.species, None if p.species in FAST_SPECIES else p.position)
            for p in s.particles), s.time < s.creation_threshold)

    def _graph(self, species, view):
        s = self.system
        transfers = s.processes.get(species, ())
        if any(type(p) not in (Migration, Formation) for p in transfers):
            raise ValueError("diffusion_average supports only Migration and Formation "
                             "as spatial processes of free I/V defects.")
        prototype = copy.copy(next(p for p in s.particles if p.species == species))
        rows, columns, values, vectors = [], [], [], []
        for origin in range(len(s.mats)):
            if view.site_occupancy[origin]:
                continue
            prototype.position = origin
            for process in transfers:
                if type(process) is not Migration:
                    continue
                cut, dx, dy, dz, r = s.neighbor_data(origin, process.neighbor_mode)
                rates = _rates(process.rate(
                    r=r, dx=dx, dy=dy, dz=dz, system=view, particle=prototype,
                    mats=s.mats[cut], matlocal=s.mats[origin], cut=cut))
                if len(rates) != len(cut):
                    raise ValueError("Migration returned an incorrect destination count.")
                positive = rates > 0
                rows.extend([origin] * int(positive.sum()))
                columns.extend(cut[positive])
                values.extend(rates[positive])
                vectors.extend(np.column_stack((dx, dy, dz))[positive])
        n = len(s.mats)
        graph = coo_matrix((values, (rows, columns)), shape=(n, n)).tocsr()
        graph.sum_duplicates()
        graph.sort_indices()
        _, labels = connected_components(graph, directed=True, connection="strong")
        channels = (np.asarray(rows, dtype=int), np.asarray(columns, dtype=int),
                    np.asarray(values, dtype=float), np.asarray(vectors, dtype=float).reshape(-1, 3))
        return graph, labels, channels

    def _groups(self, view):
        s = self.system
        groups = []
        for species in sorted(FAST_SPECIES):
            particles = [p for p in s.particles if p.species == species]
            if not particles:
                continue
            graph, labels, channels = self._graph(species, view)
            by_class = {}
            for particle in particles:
                by_class.setdefault(int(labels[particle.position]), []).append(particle)
            for label, members in by_class.items():
                sites = np.flatnonzero(labels == label)
                outside = labels[graph[sites].indices] != label
                if np.any(outside):
                    raise ValueError("Cannot average a transient diffusion class; "
                                     "the fast process must be equilibrated in a closed class.")
                restricted = graph[sites][:, sites]
                q = restricted - diags(np.asarray(restricted.sum(axis=1)).ravel())
                pi, residual = stationary_distribution(q)
                self.stats["stationary_solves"] += 1
                self.stats["max_stationary_residual"] = max(
                    self.stats["max_stationary_residual"], residual)
                probability = np.zeros(len(s.mats))
                probability[sites] = pi
                rows, columns, rates, vectors = channels
                selected = (labels[rows] == label) & (labels[columns] == label)
                indices = np.full(len(s.mats), -1, dtype=int)
                indices[sites] = np.arange(len(sites))
                local_channels = (indices[rows[selected]], indices[columns[selected]], rates[selected])
                tensor, drift = effective_transport(restricted, pi, vectors[selected], channels=local_channels)
                groups.append(MixedGroup(species, tuple(members), probability, tensor, drift))
        return tuple(groups)

    @staticmethod
    def _channel(process, initiators, counterparts, origins, destinations,
                 rates, generation=False):
        rates = _rates(rates)
        positive = rates > 0
        if not np.any(positive):
            return None
        rates = rates[positive]
        total = float(rates.sum())
        if not np.isfinite(total):
            raise ValueError("Averaged event rate overflowed.")
        return AveragedChannel(process, tuple(initiators), tuple(counterparts),
            np.asarray(origins, dtype=np.int32)[positive],
            np.asarray(destinations, dtype=np.int32)[positive],
            np.cumsum(rates), total, generation)

    def build_model(self):
        signature = self._signature()
        if signature == self.signature:
            return self.model
        s = self.system
        view = _FastView(s)
        groups = self._groups(view)
        channels = []
        eligible_pairs = 0.0
        for group in groups:
            prototype = copy.copy(group.particles[0])
            opposite = "vacancy" if group.species == "interstitial" else "interstitial"
            for process in s.processes.get(group.species, ()):
                if type(process) is not Formation:
                    continue
                for partner in groups:
                    if partner.species != opposite:
                        continue
                    origins, destinations, weights = [], [], []
                    for origin in np.flatnonzero(group.probability):
                        prototype.position = int(origin)
                        cut, dx, dy, dz, r = s.neighbor_data(origin, process.neighbor_mode)
                        intrinsic = _rates(process.rate(
                            r=r, dx=dx, dy=dy, dz=dz, system=view, particle=prototype,
                            mats=s.mats[cut], matlocal=s.mats[origin], cut=cut))
                        if len(intrinsic) != len(cut):
                            raise ValueError("Formation returned an incorrect destination count.")
                        occupancy = (len(group.particles) * len(partner.particles)
                                     * group.probability[origin] * partner.probability[cut])
                        positive = intrinsic > 0
                        eligible_pairs += float(occupancy[positive].sum())
                        origins.extend([int(origin)] * int(positive.sum()))
                        destinations.extend(cut[positive])
                        weights.extend((intrinsic * occupancy)[positive])
                    channel = self._channel(process, group.particles, partner.particles,
                                            origins, destinations, weights)
                    if channel is not None:
                        channels.append(channel)
            sites = np.flatnonzero(group.probability)
            for process in s.monomolecular.get(group.species, ()):
                intrinsic = _rates([process.rate(material=s.mats[i]) for i in sites])
                channel = self._channel(process, group.particles, (), sites, sites,
                    len(group.particles) * group.probability[sites] * intrinsic)
                if channel is not None:
                    channels.append(channel)
        for particle in s.particles:
            if particle.species in FAST_SPECIES:
                continue
            if s.processes.get(particle.species, ()):
                raise ValueError("Spatial processes of retained particles are not supported "
                                 "while diffusion_average is active.")
            events = self.event_builder(particle, s)
            rates = np.diff(np.r_[0.0, events.cumulative_rates])
            for process, destination, rate in zip(events.processes, events.destinations, rates):
                channel = self._channel(process, (particle,), (), [particle.position],
                                        [destination], [rate])
                if channel is not None:
                    channels.append(channel)
        if s.time < s.creation_threshold:
            # Fixed free particles are retained as obstacles to generation.
            for group in groups:
                if np.count_nonzero(group.probability) == 1:
                    view.occupancy[int(np.argmax(group.probability))] = 1
            for generator in s.generation:
                if type(generator) is not FP_generation:
                    raise ValueError("diffusion_average supports only FP_generation "
                                     "during active free-defect diffusion.")
                rate = _rates([generator.rate(system=view)])[0]
                if rate == 0:
                    continue
                pairs = generator.available_pairs(view)
                # Match FP_generation's constant global clock and uniform choice
                # of an available edge, neglecting dilute mobile-defect blocking.
                pairs = np.asarray(pairs, dtype=np.int32)
                channel = self._channel(generator, (), (), pairs[:, 0], pairs[:, 1],
                                        np.full(len(pairs), rate / len(pairs)), True)
                channels.append(channel)
        totals = np.asarray([c.total for c in channels])
        total = float(totals.sum())
        if not np.isfinite(total):
            raise ValueError("Total averaged event rate overflowed.")
        self.model = AveragedModel(groups, tuple(channels), np.cumsum(totals),
                                   total, eligible_pairs)
        self.signature = signature
        self.stats["model_builds"] += 1
        return self.model

    def _positions(self, model, pinned=None, reserved=()):
        """Stage a physical realization before changing the clock or particles."""
        s = self.system
        positions = dict(pinned or {})
        occupied = set(int(i) for i in np.flatnonzero(s.mats == 999))
        occupied.update(p.position for p in s.particles if p.species not in FAST_SPECIES)
        occupied.update(int(i) for i in reserved)
        for particle, site in positions.items():
            if particle.species in FAST_SPECIES:
                if site in occupied:
                    raise ValueError("Averaged event would occupy a blocked site.")
                occupied.add(site)
        remaining = [(p, group.probability) for group in model.groups
                     for p in group.particles if p not in positions]
        s.rng.shuffle(remaining)
        for particle, probability in remaining:
            for _ in range(64):
                site = int(s.rng.choice(len(probability), p=probability))
                if site not in occupied:
                    break
            else:
                available = probability.copy()
                available[list(occupied)] = 0
                if available.sum() <= 0:
                    raise ValueError("Cannot realize independent mixed defects without "
                                     "overlap; diffusion_average requires dilute occupancy.")
                available /= available.sum()
                site = int(s.rng.choice(len(available), p=available))
            positions[particle] = site
            occupied.add(site)
        return positions

    def _realize(self, positions):
        s = self.system
        for particle, position in positions.items():
            if particle.species not in FAST_SPECIES:
                continue
            old = particle.position
            if old != position:
                particle.position = position
                s.update_position(particle, old, position)
            # Paths and unwrapped displacements cannot be reconstructed from
            # equilibrium endpoint sampling. Never report fictitious MSD data.
            particle.Dx = particle.Dy = particle.Dz = float("nan")

    def _accumulate_diffusion(self, model, dt):
        if dt <= 0:
            return
        for group in model.groups:
            for particle in group.particles:
                particle.averaged_diffusion_time += dt
                particle.averaged_diffusion_integral += group.diffusion * dt
                particle.averaged_drift_integral += group.drift * dt

    def advance(self, direct_step):
        s = self.system
        while s.time < s.time_limit and (not s.particle_condition or s.particles):
            if not any(p.species in FAST_SPECIES for p in s.particles):
                return direct_step(s)
            model = self.build_model()
            boundaries = [s.time_limit]
            if s.time < s.creation_threshold and np.isfinite(s.creation_threshold):
                boundaries.append(s.creation_threshold)
            boundary = min(boundaries)
            if model.total_rate <= 0:
                if not np.isfinite(boundary):
                    return False
                positions = self._positions(model)
                self._accumulate_diffusion(model, boundary - s.time)
                s.time = boundary
                self._realize(positions)
                s.event_scheduler = None
                s._clear_event_changes()
                self.stats["boundary_steps"] += 1
                continue
            dt = -np.log(max(s.rng.random(), np.finfo(float).tiny)) / model.total_rate
            if s.time + dt >= boundary:
                positions = self._positions(model)
                self._accumulate_diffusion(model, boundary - s.time)
                s.time = boundary
                self._realize(positions)
                s.event_scheduler = None
                s._clear_event_changes()
                self.stats["boundary_steps"] += 1
                continue
            target = min(s.rng.random() * model.total_rate,
                         np.nextafter(model.total_rate, -np.inf))
            channel = model.channels[int(np.searchsorted(
                model.cumulative_rates, target, side="right"))]
            origin, destination = channel.select(s.rng)
            initiator = None
            pinned = {}
            if channel.initiators:
                initiator = channel.initiators[int(s.rng.integers(len(channel.initiators)))]
                pinned[initiator] = origin
            if channel.counterparts:
                partner = channel.counterparts[int(s.rng.integers(len(channel.counterparts)))]
                pinned[partner] = destination
            reserved = (origin, destination) if channel.generation else ()
            positions = self._positions(model, pinned, reserved)
            self._accumulate_diffusion(model, dt)
            s.time += dt
            s.IT += 1
            s.begin_event()
            self._realize(positions)
            before = tuple(s.particles)
            if channel.generation:
                channel.process.create(s, origin, destination)
            else:
                initiator.process = channel.process
                initiator.destination = destination
                channel.process.action(initiator, s, destination)
            for particle in before:
                particle.stamp_time(s)
            s.event_scheduler = None
            s._clear_event_changes()
            self.stats["averaged_events"] += 1
            return True
        return False
