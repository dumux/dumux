#!/usr/bin/env python3
# SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later
"""Build, run, and visualize the heatpipe benchmark.

The semi-analytical reference solution (the steady-state ODE system of Udell
and Fitch (1985) in the form of Huang, Kolditz and Shao (2015)) is computed by
the DuMux program test_heatpipe_odesolver.cc and written to a CSV file, which
is read here for the comparison plots.
"""

import glob
import subprocess
from pathlib import Path

import numpy as np

TARGET = "test_heatpipe_box"
REFERENCE_TARGET = "test_heatpipe_odesolver"
REFERENCE_FILE = "heatpipe_reference.csv"
BASE_NAME = "heatpipe"

# Grid resolutions compared in the convergence study; the default test (and
# CTest) always runs at 120 cells (grids/heatpipe.dgf) for a fast runtime.
# See the README for why the residual mismatch near the dry-out front shrinks
# with resolution.
RESOLUTIONS = {
    120: "grids/heatpipe.dgf",
    240: "grids/heatpipe_240.dgf",
    480: "grids/heatpipe_480.dgf",
}

SATURATION_FIELD = "S_liq"
PRESSURE_FIELD = "p_gas"
TEMPERATURE_FIELD = "T"
AIR_MOLEFRACTION_FIELD = "x^Air_gas"


def root_dir() -> Path:
    return Path(__file__).resolve().parents[4]


def case_source_dir() -> Path:
    return Path(__file__).resolve().parent


def case_build_dir() -> Path:
    return root_dir() / "build-cmake/test/porousmediumflow/2p2c/heatpipe"


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


def build_and_run_all() -> dict[int, Path]:
    """Build once, then run the grid convergence study at each resolution.

    Returns a dict mapping cell count -> path of its final VTU output.
    """
    build_dir = root_dir() / "build-cmake"
    run(["cmake", "--build", ".", "--target", TARGET, REFERENCE_TARGET], cwd=build_dir)
    run([str(case_build_dir() / REFERENCE_TARGET), "params.input"], cwd=case_build_dir())

    vtus = {}
    for cells, grid_file in RESOLUTIONS.items():
        name = f"{BASE_NAME}_{cells}"
        remove_old_outputs(name)
        run(
            [
                str(case_build_dir() / TARGET),
                "params.input",
                "-Problem.Name",
                name,
                "-Grid.File",
                str(case_source_dir() / grid_file),
            ],
            cwd=case_build_dir(),
        )
        vtus[cells] = latest_vtu(name)
    return vtus


# ---------------------------------------------------------------------------
# Semi-analytical reference solution (Udell & Fitch 1985 / Huang et al. 2015)
# ---------------------------------------------------------------------------


def semianalytical_solution(x: np.ndarray) -> dict:
    """Interpolate the semi-analytical solution written by test_heatpipe_odesolver.

    Returns a dict with arrays "Sw", "p_gas", "x_a_gas", "T" of the same shape
    as x and the dry-out position "heatpipe_length". The two-phase zone and the
    dry zone are interpolated separately, since the saturation jumps at the
    dry-out front.
    """
    data = np.genfromtxt(case_build_dir() / REFERENCE_FILE, delimiter=",", names=True)
    wet_data = data[data["two_phase"] == 1]
    dry_data = data[data["two_phase"] == 0]
    x_dry = wet_data["x"][-1]

    result = {"heatpipe_length": x_dry}
    wet = x <= x_dry
    for key, column in (("Sw", "S_liq"), ("p_gas", "p_gas"), ("x_a_gas", "x_air_gas"), ("T", "T")):
        values = np.empty_like(x)
        values[wet] = np.interp(x[wet], wet_data["x"], wet_data[column])
        values[~wet] = np.interp(x[~wet], dry_data["x"], dry_data[column])
        result[key] = values
    return result


# ---------------------------------------------------------------------------
# Post-processing / plotting
# ---------------------------------------------------------------------------


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


def point_data(mesh, field: str):
    # the box method (used by this test) writes results as point data
    if field not in mesh.point_data:
        raise KeyError(f"Missing point field '{field}'")
    x = mesh.points[:, 0]
    order = np.argsort(x)
    return x[order], mesh.point_data[field][order]


def _comparison_data(vtu: Path):
    pv, _ = require_plot_modules()
    mesh = pv.read(vtu)
    x_max = mesh.points[:, 0].max()
    x_analytic = np.linspace(0, x_max, 2000)
    analytic = semianalytical_solution(x_analytic)
    print(f"semi-analytical heat-pipe length: {analytic['heatpipe_length']:.4f} m")
    return mesh, x_max, x_analytic, analytic


def create_saturation_lineplot(vtu: Path, image_file: Path) -> None:
    """Single-panel wetting-phase saturation comparison (used as the benchmark thumbnail)."""
    _, plt = require_plot_modules()
    mesh, x_max, x_analytic, analytic = _comparison_data(vtu)

    fig, ax = plt.subplots(figsize=(8.2, 4.8), constrained_layout=True)
    x, sw = point_data(mesh, SATURATION_FIELD)
    ax.plot(x, sw, label="numerical", linewidth=2)
    ax.plot(x_analytic, analytic["Sw"], "k--", linewidth=2.4, label="semi-analytical")
    ax.set_xlabel("x [m]")
    ax.set_ylabel(r"Wetting-phase saturation $S_w$ [-]")
    ax.set_xlim(0, x_max)
    ax.set_ylim(-0.02, 1.02)
    ax.grid(True, alpha=0.3)
    ax.legend()
    fig.savefig(image_file, dpi=200)
    plt.close(fig)


def create_line_plot(vtu: Path, image_file: Path) -> None:
    _, plt = require_plot_modules()
    mesh, x_max, x_analytic, analytic = _comparison_data(vtu)

    fig, axes = plt.subplots(2, 2, figsize=(11, 7.5), constrained_layout=True)
    ax_sw, ax_T, ax_pg, ax_xa = axes[0, 0], axes[0, 1], axes[1, 0], axes[1, 1]

    x, sw = point_data(mesh, SATURATION_FIELD)
    ax_sw.plot(x, sw, label="numerical", linewidth=1.8)

    _, T = point_data(mesh, TEMPERATURE_FIELD)
    ax_T.plot(x, T, label="numerical", linewidth=1.8)

    _, pg = point_data(mesh, PRESSURE_FIELD)
    ax_pg.plot(x, pg, label="numerical", linewidth=1.8)

    _, xa = point_data(mesh, AIR_MOLEFRACTION_FIELD)
    ax_xa.plot(x, np.clip(xa, 0, None), label="numerical", linewidth=1.8)

    ax_sw.plot(x_analytic, analytic["Sw"], "k--", linewidth=2, label="semi-analytical")
    ax_T.plot(x_analytic, analytic["T"], "k--", linewidth=2, label="semi-analytical")
    ax_pg.plot(x_analytic, analytic["p_gas"], "k--", linewidth=2, label="semi-analytical")
    ax_xa.plot(x_analytic, analytic["x_a_gas"], "k--", linewidth=2, label="semi-analytical")

    ax_sw.set_ylabel(r"Wetting-phase saturation $S_w$ [-]")
    ax_T.set_ylabel("Temperature $T$ [K]")
    ax_pg.set_ylabel("Gas-phase pressure $p_g$ [Pa]")
    ax_xa.set_ylabel(r"Gas-phase air mole fraction $x_g^a$ [-]")
    for ax in (ax_sw, ax_T, ax_pg, ax_xa):
        ax.set_xlabel("x [m]")
        ax.set_xlim(0, x_max)
        ax.grid(True, alpha=0.3)
        ax.legend(fontsize=8)

    fig.savefig(image_file, dpi=200)
    plt.close(fig)


def create_saturation_image(vtu: Path, image_file: Path) -> None:
    pv, _ = require_plot_modules()
    mesh = pv.read(vtu)
    if SATURATION_FIELD not in mesh.array_names:
        raise KeyError(f"Missing field '{SATURATION_FIELD}' in {vtu}")

    plotter = pv.Plotter(
        off_screen=True, window_size=(1000, 260), border=False, border_color="white"
    )
    plotter.set_background("white")
    plotter.add_mesh(
        mesh,
        scalars=SATURATION_FIELD,
        show_edges=False,
        cmap="coolwarm",
        scalar_bar_args={
            "title": "S_w",
            "vertical": True,
            "position_x": 0.90,
            "position_y": 0.15,
            "width": 0.06,
            "height": 0.70,
        },
    )
    plotter.view_xy()
    plotter.camera.zoom(1.3)
    plotter.show(screenshot=str(image_file), auto_close=True)


def create_grid_convergence_plot(vtus: dict[int, Path], image_file: Path) -> None:
    """Wetting-phase saturation and temperature at each grid resolution, zoomed
    to the dry-out front where the resolution-dependent mismatch actually shows
    up (gas-phase pressure and air mole fraction barely change with resolution,
    so they are not repeated here)."""
    _, plt = require_plot_modules()
    pv, _ = require_plot_modules()

    finest = max(vtus)
    mesh_finest = pv.read(vtus[finest])
    x_max = mesh_finest.points[:, 0].max()
    x_analytic_full = np.linspace(0, x_max, 4000)
    x_zoom_min = max(0.0, semianalytical_solution(x_analytic_full)["heatpipe_length"] - 0.6)

    # Restrict all plotted data to the zoomed window up front, rather than plotting
    # the full domain and relying on `ax.set_xlim` to crop it: matplotlib's SVG
    # backend draws the un-cropped line past the axis bounds and hides the excess
    # via an SVG clip-path, which some SVG renderers (e.g. GitLab's markdown
    # sanitizer) strip, making the "hidden" part of the line visible again.
    x_analytic = np.linspace(x_zoom_min, x_max, 2000)
    analytic = semianalytical_solution(x_analytic)

    fig, (ax_sw, ax_T) = plt.subplots(1, 2, figsize=(11, 4.8), constrained_layout=True)
    colors = plt.cm.viridis(np.linspace(0.15, 0.85, len(vtus)))
    for color, cells in zip(colors, sorted(vtus)):
        mesh = pv.read(vtus[cells])
        x, sw = point_data(mesh, SATURATION_FIELD)
        _, T = point_data(mesh, TEMPERATURE_FIELD)
        zoom = x >= x_zoom_min
        ax_sw.plot(x[zoom], sw[zoom], color=color, linewidth=1.6, label=f"{cells} cells")
        ax_T.plot(x[zoom], T[zoom], color=color, linewidth=1.6, label=f"{cells} cells")

    ax_sw.plot(x_analytic, analytic["Sw"], "k--", linewidth=2, label="semi-analytical")
    ax_T.plot(x_analytic, analytic["T"], "k--", linewidth=2, label="semi-analytical")

    ax_sw.set_xlabel("x [m]")
    ax_sw.set_ylabel(r"Wetting-phase saturation $S_w$ [-]")
    ax_T.set_xlabel("x [m]")
    ax_T.set_ylabel("Temperature $T$ [K]")
    for ax in (ax_sw, ax_T):
        ax.set_xlim(x_zoom_min, x_max)
        ax.grid(True, alpha=0.3)
        ax.legend(fontsize=9)

    fig.savefig(image_file, dpi=200)
    plt.close(fig)


def main() -> None:
    vtus = build_and_run_all()
    finest = vtus[max(vtus)]
    out_dir = case_build_dir()

    print("Creating line plot (finest resolution)...")
    create_line_plot(finest, out_dir / f"{BASE_NAME}_lineplot_comparison.svg")
    print("Creating saturation line plot (finest resolution)...")
    create_saturation_lineplot(finest, out_dir / f"{BASE_NAME}_saturation_comparison.svg")
    print("Creating saturation field image (finest resolution)...")
    create_saturation_image(finest, out_dir / f"{BASE_NAME}_sw.png")
    print("Creating grid convergence plot...")
    create_grid_convergence_plot(vtus, out_dir / f"{BASE_NAME}_grid_convergence.svg")


if __name__ == "__main__":
    main()
