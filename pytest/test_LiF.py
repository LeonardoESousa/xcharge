import shutil
import subprocess
import re
from pathlib import Path

import numpy as np
import pandas as pd


BASE_DIR = Path(__file__).resolve().parent
W = 2e5
DX_CM = 2.887e-8
COORDINATION = 12
DT = 1 / (COORDINATION * W)
D_ANALYTIC = DX_CM**2 / (6 * DT)
RELATIVE_TOLERANCE = 5e-2


def force_single_process(input_path: Path) -> None:
    text = input_path.read_text(encoding="utf-8")
    updated = re.sub(r"^\s*n_proc\s*=\s*\d+", "n_proc = 1", text, count=1, flags=re.MULTILINE)
    if updated != text:
        input_path.write_text(updated, encoding="utf-8")


def read_data(path):
    data = pd.read_csv(path, comment="#")
    return data[data.Time.astype(str) != "END"]


def diffusion_coefficient(data):
    selected = data[
        (data["CausaMortis"] == "alive")
        & (data["Type"] == "interstitial")
    ][["Time", "DeltaX", "DeltaY", "DeltaZ"]].to_numpy(dtype=float)
    time = selected[:, 0]
    displacement_squared = np.sum(selected[:, 1:] ** 2, axis=1)
    return float(np.mean(displacement_squared / (6 * time)) * 1e-16)


def test_lif_diffusion(tmp_path):
    test_data = tmp_path / "pytest"
    test_data.mkdir()
    shutil.copy2(BASE_DIR / "LiF_input.py", test_data)
    shutil.copy2(BASE_DIR / "LiF.cif", test_data)
    force_single_process(test_data / "LiF_input.py")

    subprocess.run(
        ["kmc", "pytest/LiF_input.py"],
        cwd=tmp_path,
        check=True,
        capture_output=True,
        text=True,
    )
    calculated = diffusion_coefficient(
        read_data(tmp_path / "Simulation_LiF.txt")
    )
    relative_error = abs(calculated - D_ANALYTIC) / D_ANALYTIC
    assert relative_error < RELATIVE_TOLERANCE
