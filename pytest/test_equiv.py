from __future__ import annotations

import os
import re
import time
import uuid
import shutil
import subprocess
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Sequence

import numpy as np
import pandas as pd
import pytest
from _pytest.tmpdir import TempPathFactory


LOG_NAME_EQUIV = "input_test1.py"
CONFIG_EQUIV = [(100, 10), (125, 4), (1000, 1)]

LOG_NAME_RF = "input_test2.py"
CONFIG_RF = [(1,), (10,), (25,), (35,)]

REF_FLUOR = 3.0
FLUOR_TOL_REL = 0.20
LD_STD_MAX = 0.4
DIFF_RMS_MAX = 5e-2
HIST_GRAN = 1.5


def repo_root_from_test_file() -> Path:
    here = Path(__file__).resolve()
    for parent in [here.parent] + list(here.parents):
        if (parent / LOG_NAME_EQUIV).exists() or (parent / LOG_NAME_RF).exists():
            return parent
    return Path.cwd()


def read_simulation_csv(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, comment="#")
    if "Time" in df.columns:
        df = df[df["Time"].astype(str) != "END"]
    return df


def get_ds(df: pd.DataFrame) -> np.ndarray:
    dx = df["DeltaX"].to_numpy(dtype=float) / 10.0
    dy = df["DeltaY"].to_numpy(dtype=float) / 10.0
    dz = df["DeltaZ"].to_numpy(dtype=float) / 10.0
    return dx * dx + dy * dy + dz * dz


def all_below(values: Iterable[float], cutoff: float) -> bool:
    return all(v <= cutoff for v in values)


def spectrum(values: np.ndarray, gran: float, vmin: float, vmax: float) -> np.ndarray:
    if vmax <= vmin:
        edges = np.array([vmin, vmin + max(gran, 1e-12)], dtype=float)
    else:
        edges = np.arange(vmin, vmax + gran, gran, dtype=float)
        if edges.size < 2:
            edges = np.array([vmin, vmax], dtype=float)
    hist, _ = np.histogram(values, bins=edges, density=True)
    return hist


@dataclass(frozen=True)
class SimResult:
    output_file: Path
    duration_s: float
    args: list[str]


def _list_outputs(dirpath: Path) -> dict[str, float]:
    # If your outputs are csv with a specific extension, tighten this:
    # e.g. dirpath.glob("Simulation_*.csv")
    candidates = {}
    for p in dirpath.glob("Simulation_*"):
        if p.is_file():
            candidates[str(p.resolve())] = p.stat().st_mtime
    return candidates


def _pick_new_output(before: dict[str, float], after: dict[str, float]) -> Path | None:
    new_paths = [Path(p) for p in after.keys() if p not in before]
    if new_paths:
        return max(new_paths, key=lambda p: p.stat().st_mtime)

    changed = [Path(p) for p, mt in after.items() if before.get(p, mt) != mt]
    if changed:
        return max(changed, key=lambda p: p.stat().st_mtime)

    return None


def _kmc_args_for(log_name: str, params: Sequence[int], unique_label: str) -> list[str]:
    """
    Build argv tail for kmc depending on which input file we run.

    - input_test1.py: seems to accept a string label before numeric params.
    - input_test2.py: expects sys.argv[2] to be int, so DO NOT insert a label there.
    """
    if Path(log_name).name == LOG_NAME_EQUIV:
        # kmc <input> <label> <p1> <p2> ...
        return [unique_label, *map(str, params)]

    if Path(log_name).name == LOG_NAME_RF:
        # kmc <input> <p1> <p2> ...
        return [*map(str, params)]

    # Default safe mode: only numeric params
    return [*map(str, params)]


def _force_single_process(input_path: Path) -> None:
    text = input_path.read_text(encoding="utf-8")
    updated = re.sub(r"^\s*n_proc\s*=\s*\d+", "n_proc = 1", text, count=1, flags=re.MULTILINE)
    if updated != text:
        input_path.write_text(updated, encoding="utf-8")


def run_kmc_once(
    *,
    kmc_exe: str,
    workdir: Path,
    input_file: Path,
    log_name: str,
    params: Sequence[int],
    unique_label: str,
) -> SimResult:
    argv_tail = _kmc_args_for(log_name, params, unique_label)
    cmd = [kmc_exe, str(input_file), *argv_tail]

    before = _list_outputs(workdir)

    t0 = time.perf_counter()
    try:
        subprocess.run(
            cmd,
            cwd=str(workdir),
            check=True,
            capture_output=True,
            text=True,
        )
    except subprocess.CalledProcessError as e:
        raise AssertionError(
            "kmc failed.\n"
            f"Command: {e.cmd}\n"
            f"Return code: {e.returncode}\n"
            f"cwd: {workdir}\n"
            f"stdout:\n{e.stdout}\n"
            f"stderr:\n{e.stderr}\n"
        ) from e
    t1 = time.perf_counter()

    after = _list_outputs(workdir)
    out = _pick_new_output(before, after)
    if out is None or not out.exists():
        raise AssertionError(
            "kmc succeeded but test could not detect output file.\n"
            f"Command: {cmd}\n"
            f"cwd: {workdir}\n"
            f"Hint: adjust _list_outputs() to match your real output pattern.\n"
        )

    return SimResult(output_file=out, duration_s=(t1 - t0), args=cmd)


def run_batch(
    *,
    kmc_exe: str,
    repo_root: Path,
    tmp_path: Path,
    log_name: str,
    configs: Sequence[Sequence[int]],
    record_property,
) -> list[Path]:
    input_path = (repo_root / log_name).resolve()
    if not input_path.exists():
        pytest.skip(f"Missing {log_name} at {input_path}")

    # Run each batch in its own unique working directory to prevent overwrites
    run_tag = "pytest" + uuid.uuid4().hex[:8]
    workdir = tmp_path / f"kmc_{Path(log_name).stem}_{run_tag}"
    workdir.mkdir(parents=True, exist_ok=True)

    # Many sims rely on relative paths next to the input file
    # So copy the input file into the workdir and run from there
    local_input = workdir / Path(log_name).name
    shutil.copy2(input_path, local_input)
    _force_single_process(local_input)

    outputs: list[Path] = []
    for params in configs:
        # label only used for input_test1.py; harmless otherwise
        unique_label = f"{run_tag}_" + "_".join(map(str, params))

        res = run_kmc_once(
            kmc_exe=kmc_exe,
            workdir=workdir,
            input_file=local_input,
            log_name=log_name,
            params=list(params),
            unique_label=unique_label,
        )

        outputs.append(res.output_file)

        record_property(
            f"kmc_time_{Path(log_name).stem}_{'_'.join(map(str, params))}",
            f"{res.duration_s:.6f}s",
        )
        print(f"\n[kmc] {Path(log_name).stem} params={params} runtime: {res.duration_s:.3f} s")

    return outputs


@pytest.fixture(scope="session")
def repo_root() -> Path:
    return repo_root_from_test_file()


@pytest.fixture(scope="session")
def kmc_exe() -> str:
    return os.environ.get("KMC_EXE", "kmc")


@pytest.fixture(scope="session")
def sim_outputs_equiv(
    tmp_path_factory: TempPathFactory,
    repo_root: Path,
    kmc_exe: str,
) -> list[Path]:
    tmp_path = tmp_path_factory.mktemp("kmc_equiv")
    return run_batch(
        kmc_exe=kmc_exe,
        repo_root=repo_root,
        tmp_path=tmp_path,
        log_name=LOG_NAME_EQUIV,
        configs=CONFIG_EQUIV,
        record_property=lambda *_args, **_kwargs: None,
    )


@pytest.fixture(scope="session")
def sim_outputs_rf(
    tmp_path_factory: TempPathFactory,
    repo_root: Path,
    kmc_exe: str,
) -> list[Path]:
    tmp_path = tmp_path_factory.mktemp("kmc_rf")
    return run_batch(
        kmc_exe=kmc_exe,
        repo_root=repo_root,
        tmp_path=tmp_path,
        log_name=LOG_NAME_RF,
        configs=CONFIG_RF,
        record_property=lambda *_args, **_kwargs: None,
    )


def test_fluor_equivalence(sim_outputs_equiv: list[Path]):
    avgs: list[float] = []
    for out in sim_outputs_equiv:
        df = read_simulation_csv(out)
        fluor = df["Time"].to_numpy(dtype=float) / 1000.0
        avgs.append(float(np.mean(fluor)))

    rel_errs = [abs((av - REF_FLUOR) / REF_FLUOR) for av in avgs]
    print("\nFluor values:", avgs)
    print("Fluor rel errors:", rel_errs)
    assert all_below(rel_errs, FLUOR_TOL_REL)


@pytest.mark.parametrize("which_fixture", ["sim_outputs_equiv", "sim_outputs_rf"])
def test_ld_stability(which_fixture: str, request):
    outputs: list[Path] = request.getfixturevalue(which_fixture)
    which = "equiv" if which_fixture.endswith("equiv") else "rf"

    lds: list[float] = []
    for out in outputs:
        df = read_simulation_csv(out)
        lds.append(float(np.sqrt(np.mean(get_ds(df)))))

    std = float(np.std(lds))
    print(f"\nLD ({which}):", lds)
    print(f"LD std ({which}):", std)
    assert std <= LD_STD_MAX


@pytest.mark.parametrize("which_fixture", ["sim_outputs_equiv", "sim_outputs_rf"])
def test_diffusion_distribution(which_fixture: str, request):
    outputs: list[Path] = request.getfixturevalue(which_fixture)
    which = "equiv" if which_fixture.endswith("equiv") else "rf"

    if len(outputs) < 2:
        pytest.skip("Need at least two simulations to compare diffusion distributions")

    ds_list = []
    for out in outputs:
        df = read_simulation_csv(out)
        ds_list.append(get_ds(df))

    vmin = float(min(np.min(ds) for ds in ds_list))
    vmax = float(max(np.max(ds) for ds in ds_list))

    hists = [spectrum(ds, HIST_GRAN, vmin, vmax) for ds in ds_list]
    ref = hists[0]

    rms_arr: list[float] = []
    for i in range(1, len(hists)):
        other = hists[i]
        mask = (ref != 0) & (other != 0)
        if not np.any(mask):
            pytest.fail("Histogram overlap mask is empty (no shared non-zero bins)")
        rms_arr.append(float(np.sqrt(np.mean((ref[mask] - other[mask]) ** 2))))

    worst = float(max(rms_arr))
    print(f"\nDiffusion RMS ({which}):", worst)
    assert worst <= DIFF_RMS_MAX
