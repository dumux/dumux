#!/usr/bin/env python3
# SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later
"""Build, run, and visualize the 1pni heat conduction benchmark."""

import glob
import subprocess
import xml.etree.ElementTree as ET
from dataclasses import dataclass
from pathlib import Path


TARGET = "test_1pni_conduction_tpfa"
INPUT_FILE = "params_conduction.input"
BASE_NAME = "1pni_1d_conduction_benchmark"
NUM_PROFILES = 5
TEMPERATURE_FIELD = "T"
EXACT_FIELD = "temperatureExact"

# the thermal front stays well within this part of the 5 m domain
X_MAX = 1.5
TEMPERATURE_LIMITS = (289.8, 300.2)

# ordinal blue ramp (light -> dark) encoding the ordered output times
PROFILE_COLORS = ("#86b6ef", "#5598e7", "#2a78d6", "#1c5cab", "#104281")


@dataclass(frozen=True)
class Result:
    time: float
    name: str
    vtu: Path


def root_dir() -> Path:
    return Path(__file__).resolve().parents[4]


def case_build_dir() -> Path:
    return root_dir() / "build-cmake/test/porousmediumflow/1p/nonisothermal"


def run(command: list[str], cwd: Path | None = None) -> None:
    print(f"+ {' '.join(command)}")
    subprocess.run(command, cwd=cwd, check=True)


def remove_old_outputs(name: str) -> None:
    for pattern in (f"{name}-*.vtu", f"{name}.pvd", f"{name}-*.pvtu"):
        for file_name in glob.glob(str(case_build_dir() / pattern)):
            Path(file_name).unlink()


def output_times(name: str) -> list[Result]:
    """Read the output times and the corresponding VTU files from the collection file."""
    collection = case_build_dir() / f"{name}.pvd"
    if not collection.exists():
        raise FileNotFoundError(f"No VTU output found for {name}")

    results = [
        Result(time=float(data_set.attrib["timestep"]),
               name=name,
               vtu=case_build_dir() / data_set.attrib["file"])
        for data_set in ET.parse(collection).getroot().iter("DataSet")
    ]

    return sorted(results, key=lambda result: result.time)


def representative_times(results: list[Result], count: int) -> list[Result]:
    """Pick a logarithmic spread of output times.

    The initial output is skipped: the exact solution is only updated inside the
    time loop and is therefore still zero in the first output file.
    """
    transient = [result for result in results if result.time > 0.0]
    if not transient:
        raise ValueError("The output contains no time step beyond the initial condition")
    if len(transient) <= count:
        return transient

    first, last = transient[0].time, transient[-1].time
    targets = [first * (last / first) ** (index / (count - 1)) for index in range(count)]
    picked = sorted({min(range(len(transient)), key=lambda i: abs(transient[i].time - target))
                     for target in targets})

    return [transient[index] for index in picked]


def build_and_run() -> list[Result]:
    build_dir = root_dir() / "build-cmake"
    run(["make", TARGET], cwd=build_dir)

    remove_old_outputs(BASE_NAME)
    run([
        str(case_build_dir() / TARGET),
        INPUT_FILE,
        "-Problem.Name", BASE_NAME,
    ], cwd=case_build_dir())

    return output_times(BASE_NAME)


def require_plot_modules():
    try:
        import matplotlib.pyplot as plt
        import pyvista as pv
    except ImportError as error:
        raise SystemExit(
            "Post-processing requires PyVista and Matplotlib. "
            "Install them with: python3 -m pip install pyvista matplotlib"
        ) from error

    return pv, plt


def sorted_profile(mesh, field: str):
    """Return the x-coordinates and values of a cell or vertex field, sorted along x.

    Vertex data of the quasi-1D grid carries several values per x-position, which
    are averaged, so that box and cell-centered schemes are handled alike.
    """
    if field in mesh.cell_data:
        positions, values = mesh.cell_centers().points[:, 0], mesh.cell_data[field]
    elif field in mesh.point_data:
        positions, values = mesh.points[:, 0], mesh.point_data[field]
    else:
        raise KeyError(f"Missing field '{field}'")

    profile: dict[float, list[float]] = {}
    for position, value in zip(positions, values):
        profile.setdefault(round(float(position), 9), []).append(float(value))

    x = sorted(profile)
    return x, [sum(profile[position]) / len(profile[position]) for position in x]


def format_time(time: float) -> str:
    for limit, factor, unit in ((60.0, 1.0, "s"), (3600.0, 60.0, "min"), (86400.0, 3600.0, "h")):
        if time < limit:
            return f"{time / factor:.3g} {unit}"
    return f"{time / 86400.0:.3g} d"


def create_line_plot(results: list[Result], image_file: Path) -> None:
    pv, plt = require_plot_modules()

    fig, ax = plt.subplots(figsize=(8.2, 4.8), constrained_layout=True)

    for index, result in enumerate(results):
        mesh = pv.read(result.vtu)
        x, temperature = sorted_profile(mesh, TEMPERATURE_FIELD)
        exact_x, exact = sorted_profile(mesh, EXACT_FIELD)
        color = PROFILE_COLORS[index % len(PROFILE_COLORS)]
        ax.plot(x, temperature, label=rf"$T$ numerical ($t$ = {format_time(result.time)})",
                linewidth=2, color=color)
        ax.plot(exact_x, exact, linewidth=1.4, color="black", linestyle="--")

    ax.plot([], [], label=r"$T$ exact", linewidth=1.4, color="black", linestyle="--")
    ax.set_xlabel("x [m]")
    ax.set_ylabel(r"Temperature $T$ [K]")
    ax.set_xlim(0, X_MAX)
    ax.set_ylim(*TEMPERATURE_LIMITS)
    ax.grid(True, alpha=0.3)
    ax.legend()
    fig.savefig(image_file, dpi=200)
    plt.close(fig)


def main() -> None:
    results = representative_times(build_and_run(), NUM_PROFILES)
    out_dir = case_build_dir()

    print("Creating line plot...")
    create_line_plot(results, out_dir / f"{BASE_NAME}_lineplot.png")


if __name__ == "__main__":
    main()
