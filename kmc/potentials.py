"""Optional spatial energy profiles used by KMC rate models."""

import numpy as np


BOLTZMANN_EV = 8.617333262e-5


class MigrationBarrier:
    """Precompute a one-dimensional migration-barrier correction.

    ``function(coordinate, period)`` must return an energy profile, in eV, at
    every supplied coordinate. For a periodic simulation, coordinates are
    wrapped into ``[0, period)`` before the function is called.

    With ``coupling="destination"``, the profile is interpreted as an absolute
    correction to the destination barrier. With ``coupling="symmetric"``, it
    is interpreted as a site-energy landscape and half of the destination
    minus origin energy is added to the barrier. The latter preserves the
    supplied rate in spatially flat regions and obeys detailed balance for
    symmetric unmodified hop rates.

    If ``directional`` is true in destination mode, angular weighting applies
    to the symmetric bond part (the mean endpoint correction), while the
    half energy difference remains unweighted. This preserves detailed balance
    on reciprocal graphs with symmetric base rates. In symmetric mode the
    profile supplies only a site-energy difference, so it is never angularly
    weighted. Scalar site profiles are nondirectional by default.

    Xcharge evaluates this assigner after the morphology has been constructed.
    """

    def __init__(
        self,
        function,
        temperature,
        axis="x",
        directional=False,
        angular_exponent=2.0,
        coupling="destination",
    ):
        if not callable(function):
            raise TypeError("MigrationBarrier function must be callable.")
        if temperature <= 0:
            raise ValueError("MigrationBarrier temperature must be positive.")
        if axis not in {"x", "y", "z"}:
            raise ValueError("MigrationBarrier axis must be 'x', 'y', or 'z'.")
        if not isinstance(directional, (bool, np.bool_)):
            raise TypeError("MigrationBarrier directional must be boolean.")
        if not np.isfinite(angular_exponent) or angular_exponent <= 0:
            raise ValueError(
                "MigrationBarrier angular_exponent must be positive."
            )
        if coupling not in {"destination", "symmetric"}:
            raise ValueError(
                "MigrationBarrier coupling must be 'destination' or "
                "'symmetric'."
            )

        self.function = function
        self.temperature = float(temperature)
        self.axis = axis
        self.directional = bool(directional)
        self.angular_exponent = float(angular_exponent)
        self.coupling = coupling
        self.requires_morphology = True

    def assign_to_system(self, system):
        if not hasattr(system, "R"):
            raise RuntimeError(
                "MigrationBarrier requires lattice coordinates."
            )

        axis_index = {"x": 0, "y": 1, "z": 2}[self.axis]
        coordinates = np.asarray(system.R[:, axis_index], dtype=float)

        if system.cell is not None:
            vector = np.asarray(system.cell[axis_index], dtype=float)
            other_components = np.delete(vector, axis_index)
            if np.any(np.abs(other_components) > 1.0e-10):
                raise ValueError(
                    "A Cartesian MigrationBarrier currently requires an "
                    "orthorhombic cell."
                )
            period = abs(float(vector[axis_index]))
        else:
            period = float(np.ptp(coordinates))

        if period <= 0:
            raise ValueError(
                f"The lattice has zero extent along the {self.axis!r} axis."
            )

        if getattr(system, "periodic", False):
            evaluation_coordinates = np.mod(coordinates, period)
        else:
            evaluation_coordinates = coordinates

        values = np.asarray(
            self.function(evaluation_coordinates, period),
            dtype=float,
        )
        if values.ndim == 0:
            values = np.full(coordinates.shape, float(values), dtype=float)
        if values.shape != coordinates.shape:
            raise ValueError(
                "MigrationBarrier function must return either one value or "
                "one value per lattice site."
            )
        if np.any(~np.isfinite(values)):
            raise ValueError("MigrationBarrier produced a non-finite value.")

        system.migration_barrier_change = values
        system.migration_kbt = BOLTZMANN_EV * self.temperature
        system.migration_barrier_axis = axis_index
        system.migration_barrier_directional = self.directional
        system.migration_barrier_angular_exponent = self.angular_exponent
        system.migration_barrier_coupling = self.coupling
