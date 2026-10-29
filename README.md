# Xcharge

**Kinetic Monte Carlo simulations of defects, excitons, and charge carriers.**

Xcharge models stochastic migration and reactions on spatial lattices. It supports
vacancies, interstitials, Frenkel pairs, singlet and triplet excitons, electrons,
and holes. Morphologies can be generated in Python or read from a crystal
structure in CIF format.

Version **1.0.0** is distributed as the Python package `KMC` and provides the
`kmc` command. Input files import the `kmc` module.

## Capabilities

- Defect migration, Frenkel-pair generation, formation, dissociation, and
  annihilation; interstitial-pair formation and dissociation.
- First-neighbor defect processes, with separate rates for bulk and grain-boundary
  environments.
- Periodic boundary conditions using the full lattice cell, including triclinic
  cells; cumulative displacements across periodic boundaries.
- Optional spatial energy profiles for migration, such as attractive GB wells.
- Exciton transfer through Förster and Dexter processes, fluorescence,
  phosphorescence, intersystem crossing, and annihilation.
- Charge hopping through Miller–Abrahams rates and charge dissociation.
- Independent simulation rounds across multiple CPU processes, reproducible
  random seeds, text logs, and animations.

Processes are enabled explicitly in the input. Defining a rate object alone does
not activate that mechanism.

## Installation
Python 3.10 or later is required. Xcharge supports NumPy 1.26 and 2.x; compatible dependencies are installed automatically by pip.

```bash
git clone https://github.com/LeonardoESousa/xcharge.git
cd xcharge
python -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip setuptools wheel
python -m pip install .
```

On Windows, activate the environment with `.venv\Scripts\Activate.ps1` in
PowerShell.

An optional C extension accelerates `kmc.utils`. Installation uses its Python
implementation if compilation is unavailable. Check which implementation is loaded:

```bash
python -c "import kmc.utils; print(kmc.utils.__file__)"
```

A compiled module ends in `.so` on Linux or `.pyd` on Windows. Build it with a
compiler suitable for the machines where jobs will run.

After updating the source, reinstall it in the environment used by your jobs:

```bash
python -m pip install --no-deps --force-reinstall .
```

## Quick start

Run the self-contained migration-barrier example from the repository directory:

```bash
kmc input_examples/migration_barrier.py
```

The equivalent module invocation is:

```bash
python -m kmc input_examples/migration_barrier.py
```

This example runs four rounds with a vacancy and an interstitial on a periodic
lattice and writes `Simulation_migration_barrier.txt` to the working directory.
The other [input examples](input_examples/) cover excitons, charges, and defects;
some require a structure file or additional command-line arguments.

Inputs are Python programs and can read extra arguments, for example:

```bash
kmc base.py 100
```

The meaning of `100` is defined by `base.py`.

## Input structure

An input defines the morphology, initial particles, enabled processes, and run
settings. These are the main controls:

| Setting | Meaning |
| --- | --- |
| `identifier` | Output filename suffix |
| `rounds` | Number of independent realizations |
| `n_proc` | Number of worker processes on one machine |
| `frozen_lattice` | Reuse the morphology between rounds |
| `periodic` | Enable periodic boundary conditions |
| `time_limit` | Maximum simulated time per round |
| `time_units` | Time label written to the log; default is `"ps"` |
| `random_seed` | Optional seed for reproducible rounds |
| `append_output` | Append new rounds to an existing log; default is `False` |
| `animation_mode` | Display the dynamics instead of running the batch workflow |

**Keep units consistent.** `time_units` labels the output; it does not convert
rates or lifetimes. With CIF coordinates, distances are in ångströms. For rates in
s⁻¹, use times in seconds and set `time_units = "s"`. An Eyring rate can be defined as:

```python
import numpy as np

T = 298.0
kB = 8.617333262e-5  # eV/K
h = 4.135667696e-15  # eV s


def eyring_rate(barrier_eV):
    return (kB * T / h) * np.exp(-barrier_eV / (kB * T))
```

Register spatial processes in `processes`, local reactions in `monomolecular`,
and ongoing particle generation in `generation`. For example:

```python
processes = {
    "vacancy": [vacancy_migration, fp_formation],
    "interstitial": [interstitial_migration],
    "frenkelpair0": [],
    "frenkelpair2": [],
}
monomolecular = {
    "vacancy": [],
    "interstitial": [],
    "frenkelpair0": [fp0_annihilation, fp0_dissociation],
    "frenkelpair2": [fp2_annihilation, fp2_dissociation],
}
generation = []
```

Here, FP formation is initiated only by the vacancy. `generation = []` disables
ongoing generation; initial particles are assigned separately. The process objects
in this snippet must be constructed with the rates appropriate to your model.

### CIF structures and GB labels

```python
from pathlib import Path
import kmc.morphology as morphology

material_label = {"Cs": 0, "Br": 1, "Cs_GB": 2, "Br_GB": 3}
structure_file = Path(__file__).resolve().with_name("grain_boundary.cif")
lattice_func = morphology.ReadCIF(
    str(structure_file), [1, 1, 1], material_label, ["Pb"]
)
```

CIF labels such as `Br123_GB` are normalized to `Br_GB` for material indexing.
These labels select bulk or GB rates. Chemical species define the neighbor
geometry: `Br` and `Br_GB` belong to the same physical sublattice.

`Migration`, `Formation`, `I2Formation`, and `FP_generation` use precomputed
first-neighbor shells. `first_neighbor_rtol` groups distorted shell distances
with a relative tolerance (default `0.15`). Same-sublattice migration edges have
matching reverse edges. The general `cutoff` applies to distance-dependent
processes; it does not set the range of these first-neighbor defect processes.

### Migration energy profiles

To represent a GB well, supply a vectorized site-energy function in eV:

```python
from kmc.potentials import MigrationBarrier


def gb_energy(y, cell_period):
    spacing = 164.597  # Å; set this for the actual structure
    shape = 0.5 * (1.0 + np.cos(2.0 * np.pi * y / spacing))
    return 0.152283 * (1.0 - shape**17.144)


migration_barrier = MigrationBarrier(
    function=gb_energy, temperature=T, axis="y", coupling="symmetric"
)
```

The coordinate origin must align with the GB planes, and the profile must match
the structure's periodicity. Spatial profiles currently require an orthorhombic
cell aligned with the chosen Cartesian axis.

With `coupling="symmetric"`, a hop from site *i* to *j* receives the barrier
correction `(U[j] - U[i]) / 2`. For symmetric base rates this gives detailed
balance. A profile depending only on `y` leaves hops with unchanged `y` unbiased.
No angular weighting is applied to site-energy differences.

The alternative `coupling="destination"` adds an absolute destination barrier
correction; its optional angular weighting has a different interpretation from
a conservative site-energy well. One registered profile applies to both vacancy
and interstitial `Migration` processes. The base material-dependent rates remain
in effect.

## Parallel rounds and SLURM

`n_proc` distributes independent rounds within one node. To use several nodes,
submit independent jobs or a SLURM array with distinct output identifiers.
For four jobs of 250 rounds each, include in the input:

```python
import os

chunk_id = int(os.environ.get("SLURM_ARRAY_TASK_ID", 0))
identifier = f"my_simulation_chunk{chunk_id}"
rounds = 250
n_proc = 40
random_seed = 12345 + chunk_id  # distinct reproducible seed for each job
```

A corresponding submission script is:

```bash
#!/bin/bash
#SBATCH --time=2-0
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=40
#SBATCH --array=0-3
#SBATCH --output=slurm-%A_%a.out

set -euo pipefail
source /path/to/venv/bin/activate
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1

kmc base.py 100
```

Set the partition and load the Python environment required by your cluster.
Ensure each task can access its input and CIF files. Give different parameter
sets distinct identifiers as well. Appending logs adds independent realizations;
it does **not** resume an interrupted trajectory.

## Output and analysis

Each batch writes `Simulation_<identifier>.txt`. It has comment headers,
CSV-formatted records, and an `END` separator after each round. Records include
the simulation time, cumulative displacement, particle type, site information,
process, and final status.

`DeltaX`, `DeltaY`, and `DeltaZ` are signed cumulative hop displacements, including
periodic crossings. `FinalX`, `FinalY`, and `FinalZ` are coordinates inside the
simulation cell. `Time` is measured from the start of the round, not necessarily
from the creation of an individual defect. The reported `Energy` is an assigned
site energy, not a migration activation barrier.

For diffusion analysis, use elapsed particle time and unwrapped displacements.
In GB structures, examine directional mean-squared displacements and their
linear-time regime before interpreting a single scalar diffusion coefficient.

To display an animation, set `animation_mode = True`. To save it, also set
`save_animation = True` and `animation_exten = "gif"` or `"mp4"`. MP4 output
requires an available video writer such as FFmpeg.

## Development and tests

```bash
python -m pip install -e ".[dev]"
python -m pytest -q
```

Custom rate models can provide `rate(...)` and `action(...)` methods, following
the implementations in [kmc/rates.py](kmc/rates.py). Keep new mechanisms covered by
tests of their physical behavior and particle/lattice bookkeeping.

## Authors and license

Developed by **Leonardo Evaristo de Sousa** and **Tiago de Sousa Araújo Cassiano**.
See the [project wiki](https://github.com/LeonardoESousa/xcharge/wiki) for additional
background and older workflows.

Licensed under the [GNU General Public License, version 3 or later](LICENSE).
