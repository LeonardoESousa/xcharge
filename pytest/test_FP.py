import shutil
import subprocess
import re
from pathlib import Path

import numpy as np
import pandas as pd


BASE_DIR = Path(__file__).resolve().parent
CONFIGS = [10**power for power in (-5, -6, -7)]
THRESHOLD = 1e-2


def force_single_process(input_path: Path) -> None:
    text = input_path.read_text(encoding="utf-8")
    updated = re.sub(r"^\s*n_proc\s*=\s*\d+", "n_proc = 1", text, count=1, flags=re.MULTILINE)
    if updated != text:
        input_path.write_text(updated, encoding="utf-8")


def read_data(path):
    data = pd.read_csv(path, comment="#")
    return data[data.Time.astype(str) != "END"]


def average_generation_population(path):
    data = read_data(path)
    times = sorted(
        data[data["CausaMortis"] == "generated"]["Time"].to_numpy(dtype=float)
    )
    population = np.arange(1, len(times) + 1, dtype=float) / 1000
    return float(np.mean(population / times))


def test_frenkel_pair_generation_rate(tmp_path):
    test_data = tmp_path / "pytest"
    test_data.mkdir()
    shutil.copy2(BASE_DIR / "FP_generation_input.py", test_data)
    shutil.copy2(BASE_DIR / "cifexample.cif", test_data)
    force_single_process(test_data / "FP_generation_input.py")

    estimates = []
    for configured_rate in CONFIGS:
        subprocess.run(
            ["kmc", "pytest/FP_generation_input.py", str(configured_rate)],
            cwd=tmp_path,
            check=True,
            capture_output=True,
            text=True,
        )
        output = tmp_path / f"Simulation_generation_{configured_rate}.txt"
        estimates.append(average_generation_population(output))

    rms = float(np.sqrt(np.mean((np.asarray(CONFIGS) - estimates) ** 2)))
    assert rms < THRESHOLD
