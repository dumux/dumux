#!/usr/bin/env python3
# SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later
"""Build, run, and visualize the Buckley-Leverett benchmark."""

import glob
import subprocess
from dataclasses import dataclass
from pathlib import Path


TARGET = "test_2p_buckleyleverett_tpfa"
BASE_NAME = "buckleyleverett"
CELLS = (200, 400)
SATURATION_FIELD = "S_aq"
EXACT_FIELD = "Sw_exact"

@dataclass(frozen=True)
class Result:
    cells: int
    name: str
    vtu: Path


def root_dir() -> Path:
    return Path(__file__).resolve().parents[4]


def case_build_dir() -> Path:
    return root_dir() / "build-cmake/test/porousmediumflow/2p/buckleyleverett"


def run(command: list[str], cwd: Path | None = None) -> None:
    print(f"+ {' '.join(command)}")
    subprocess.run(command, cwd=cwd, check=True)


def remove_old_outputs(name: str) -> None:
    for pattern in (f"{name}-*.vtu", f"{name}.pvd", f"{name}-*.pvtu"):
        for file_name in glob.glob(str(case_build_dir() / pattern)):
            Path(file_name).unlink()


def latest_vtu(name: str) -> Path:
    files = sorted(case_build_dir().glob(f"{name}-*.vtu"))
    if not files:
        raise FileNotFoundError(f"No VTU output found for {name}")
    return files[-1]


def build_and_run() -> list[Result]:
    build_dir = root_dir() / "build-cmake"
    run(["make", TARGET], cwd=build_dir)

    results = []
    for cells in CELLS:
        name = f"{BASE_NAME}_{cells}x1"
        remove_old_outputs(name)
        run([
            str(case_build_dir() / TARGET),
            "params.input",
            "-Problem.Name", name,
            "-Grid.Cells", f"{cells} 1",
        ], cwd=case_build_dir())
        results.append(Result(cells=cells, name=name, vtu=latest_vtu(name)))

    return results


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


def cell_data(mesh, field: str):
    if field not in mesh.cell_data:
        raise KeyError(f"Missing cell field '{field}'")
    return mesh.cell_centers().points[:, 0], mesh.cell_data[field]


def create_line_plot(results: list[Result], image_file: Path) -> None:
    pv, plt = require_plot_modules()

    fig, ax = plt.subplots(figsize=(8.2, 4.8), constrained_layout=True)
    exact_x = exact = None

    for result in results:
        mesh = pv.read(result.vtu)
        # x, sw = sorted_profile(mesh, SATURATION_FIELD)
        # exact_x, exact = sorted_profile(mesh, EXACT_FIELD)
        x, sw = cell_data(mesh, SATURATION_FIELD)
        exact_x, exact = cell_data(mesh, EXACT_FIELD)
        ax.plot(x, sw, label=rf"$S_w$ numerical ({result.cells}x1)", linewidth=2)

    ax.plot(exact_x, exact, label=r"$S_w$ exact", linewidth=2.4, color="black", linestyle="--")
    ax.set_xlabel("x [m]")
    ax.set_ylabel(r"Wetting-phase saturation $S_w$ [-]")
    ax.set_xlim(0, 100)
    ax.set_ylim(-0.02, 1.02)
    ax.grid(True, alpha=0.3)
    ax.legend()
    fig.savefig(image_file, dpi=200)
    plt.close(fig)


def create_saturation_image(result: Result, image_file: Path) -> None:
    pv, _ = require_plot_modules()
    mesh = pv.read(result.vtu)
    if SATURATION_FIELD not in mesh.array_names:
        raise KeyError(f"Missing field '{SATURATION_FIELD}' in {result.vtu}")

    plotter = pv.Plotter(off_screen=True, window_size=(650, 400), border=False, border_color="white")
    plotter.set_background("white")
    plotter.add_mesh(
        mesh,
        scalars=SATURATION_FIELD,
        show_edges=False,
        cmap="coolwarm",
        scalar_bar_args={
            "title": "S_w",
            "vertical": True,
            "position_x": 0.86,
            "position_y": 0.20,
            "width": 0.08,
            "height": 0.60,
        },
    )
    plotter.view_xy()
    plotter.camera.zoom(1.2)
    plotter.show(screenshot=str(image_file), auto_close=True)


def main() -> None:
    results = build_and_run()
    fine = max(results, key=lambda result: result.cells)
    out_dir = case_build_dir()

    print("Creating line plot...")
    create_line_plot(results, out_dir / f"{BASE_NAME}_lineplot_comparison.png")
    print("Creating solution image...")
    create_saturation_image(fine, out_dir / f"{BASE_NAME}_sw.png")


if __name__ == "__main__":
    main()
