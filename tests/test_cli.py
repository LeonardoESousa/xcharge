import csv
import os
import subprocess
import sys
from pathlib import Path


CONFIG = """
import numpy as np

import kmc.morphology as morphology
from kmc.rates import Fluor

identifier = "smoke"
rounds = 6
n_proc = 1
frozen_lattice = True
periodic = False
cutoff = 1.1
time_limit = np.inf
random_seed = 2468
animation_mode = False
bimolec = False

lifetimes = {0: 3.0}
fluorescence = Fluor(life=lifetimes)
processes = {"singlet": []}
monomolecular = {"singlet": [fluorescence]}

lattice_func = morphology.Lattice(
    num_sites=4,
    vector=[1.0, 0.0, 0.0],
    disorder=[0.0, 0.0, 0.0],
    composition=[1.0],
)
s1 = morphology.Gaussian_energy({0: (2.5, 0.0), "level": "s1"})
exciton = morphology.Create_Particles(
    "singlet", 1, morphology.randomized, mat=[0]
)
"""


def run_simulation(workdir: Path):
    input_path = workdir / "input.py"
    input_path.write_text(CONFIG, encoding="utf-8")
    repo_root = Path(__file__).resolve().parents[1]
    env = os.environ.copy()
    old_pythonpath = env.get("PYTHONPATH")
    env["PYTHONPATH"] = (
        str(repo_root)
        if not old_pythonpath
        else os.pathsep.join((str(repo_root), old_pythonpath))
    )
    completed = subprocess.run(
        [sys.executable, "-m", "kmc", str(input_path)],
        cwd=workdir,
        env=env,
        check=True,
        capture_output=True,
        text=True,
        timeout=20,
    )
    output = workdir / "Simulation_smoke.txt"
    assert output.exists(), completed.stdout + completed.stderr
    with output.open(encoding="utf-8", newline="") as handle:
        rows = [
            row
            for row in csv.DictReader(
                line for line in handle if not line.startswith("#")
            )
            if row["Time"] != "END"
        ]
    return rows


def test_seeded_cli_run_is_small_complete_and_reproducible(tmp_path):
    first_dir = tmp_path / "first"
    second_dir = tmp_path / "second"
    first_dir.mkdir()
    second_dir.mkdir()

    first = run_simulation(first_dir)
    second = run_simulation(second_dir)

    assert len(first) == 6
    assert {row["CausaMortis"] for row in first} == {"fluor"}
    assert all(float(row["Time"]) > 0 for row in first)
    assert first == second
