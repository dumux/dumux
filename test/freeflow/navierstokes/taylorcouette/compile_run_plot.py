#!/usr/bin/env python3
# SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later
"""
compile_run_plot.py
===================

Builds the Taylor-Couette benchmark, runs it on two grids and compares both
numerical solutions with the analytical solution:

  coarse grid   params.input (the grid used by the regression test)
  refined grid  params.input with one global refinement (-Grid.Refinement 1),
                i.e. every cell is split into four

Output, written to the build directory of the test:

  analytical_comparison.png   radial velocity and pressure profiles
  l2_errors.md                table with the relative L2 errors of both runs

The analytical solution is read from the pressureExact and velocityExact
fields of the VTK output, i.e. it is evaluated by the problem class in C++.

With --update-readme, the table in README.md (between the markers
<!-- L2-ERRORS-START --> and <!-- L2-ERRORS-END -->) is replaced by the new one
and the figure is copied to images/analytical_comparison.png.

Requires PyVista and Matplotlib.

Usage:
  python3 compile_run_plot.py [--build-dir BUILD_DIR] [--skip-build] [--update-readme]
"""

from __future__ import annotations

import argparse
import json
import shutil
import subprocess
from dataclasses import dataclass
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pyvista as pv

SCRIPT_DIR = Path(__file__).resolve().parent
TARGET = "test_ff_navierstokes_taylorcouette"
TEST_SUBDIR = Path("test") / "freeflow" / "navierstokes" / "taylorcouette"

PARAMS_FILE = "params.input"

# (Problem.Name, number of global refinements) -- the first entry is the coarse grid
CASES = [
    ("test_ff_taylorcouette", 0),
    ("test_ff_taylorcouette_refined", 1),
]

R1, R2 = 1.0, 2.0  # inner and outer radius
N_BINS = 80        # radial bins used to average the numerical solution

# plot colors: analytical solution in ink, the two grids in categorical slots 1 and 2
INK = "#0b0b0b"
INK_SECONDARY = "#52514e"
SURFACE = "#fcfcfb"
# coarse grid: open ring, refined grid: small filled square on top, so both
# stay visible where the two solutions (nearly) coincide
MARKER_STYLES = [
    {"marker": "o", "ms": 7, "mfc": "none", "mec": "#2a78d6", "mew": 1.4},
    {"marker": "s", "ms": 4, "mfc": "#eb6834", "mec": SURFACE, "mew": 0.8},
]


@dataclass
class Result:
    n_radial: int
    n_angular: int
    n_cells: int
    metrics: dict
    data: dict


# =============================================================================
# Build / run
# =============================================================================


def find_build_dir(explicit: Path | None) -> Path:
    """Locate the build directory (a `build-cmake` folder in the DuMux root)."""

    if explicit is not None:
        return explicit.resolve()

    for parent in SCRIPT_DIR.parents:
        candidate = parent / "build-cmake"
        if candidate.is_dir():
            return candidate

    raise RuntimeError(
        "Could not find a build-cmake directory automatically. "
        "Pass it with --build-dir."
    )


def build(build_dir: Path) -> None:
    print(f"Building {TARGET} in {build_dir} ...")
    subprocess.run(
        ["cmake", "--build", str(build_dir), "--target", TARGET],
        check=True,
    )


def grid_resolution(params_path: Path) -> tuple[int, int]:
    """Radial (summed over all zones) and angular number of cells from a parameter file."""

    values = {}

    for line in params_path.read_text().splitlines():
        key, sep, value = line.split("#", 1)[0].partition("=")
        if sep and key.strip() in ("Cells0", "Cells1"):
            values[key.strip()] = [int(token) for token in value.split()]

    return sum(values["Cells0"]), sum(values["Cells1"])


def run_case(exe_dir: Path, params_path: Path, problem_name: str, refinement: int):
    """
    Run the benchmark with the given parameter file and number of global
    refinements. Returns the path of the
    final VTU file and the relative L2 errors written by the program.
    """

    exe = exe_dir / TARGET
    if not exe.exists():
        raise RuntimeError(f"Executable not found: {exe}")

    # every run overwrites solution_metrics.json, so read it right after the run
    metrics_file = exe_dir / "solution_metrics.json"
    metrics_file.unlink(missing_ok=True)

    print(f"Running {TARGET} with {params_path.name}, {refinement} global refinement(s) ...")
    subprocess.run(
        [
            str(exe), str(params_path),
            "-Problem.Name", problem_name,
            "-Grid.Refinement", str(refinement),
        ],
        cwd=exe_dir,
        check=True,
    )

    metrics = json.loads(metrics_file.read_text())

    return exe_dir / f"{problem_name}-00001.vtu", metrics


# =============================================================================
# Data extraction
# =============================================================================


def load(vtu_path: Path) -> dict:
    """Read the cell data of a VTU file."""

    mesh = pv.read(vtu_path)

    def field(name):
        if name not in mesh.cell_data:
            raise KeyError(
                f"Field '{name}' not found in {vtu_path.name}, "
                f"available: {list(mesh.cell_data.keys())}"
            )
        return np.asarray(mesh.cell_data[name])

    centers = mesh.cell_centers().points

    return {
        "n_cells": mesh.n_cells,
        "cx": centers[:, 0],
        "cy": centers[:, 1],
        "p": field("p"),
        "u": field("velocity_liq (m/s)"),
        "p_exact": field("pressureExact"),
        "u_exact": field("velocityExact"),
    }


# =============================================================================
# Radial profiles
# =============================================================================


def _bin_radially(r, values, bins):
    idx = np.clip(np.digitize(r, bins) - 1, 0, len(bins) - 2)
    return np.array([
        values[idx == i].mean() if (idx == i).any() else np.nan
        for i in range(len(bins) - 1)
    ])


def numerical_profile(data: dict):
    """
    Radial profile of the numerical solution: cell values averaged in radial
    bins (the flow is axisymmetric, so this only averages out the angular
    direction). Pressure is shifted to zero at the inner cylinder.
    """

    r = np.hypot(data["cx"], data["cy"])
    inside = (r >= R1) & (r <= R2)

    speed = np.linalg.norm(data["u"][inside][:, :2], axis=1)

    bins = np.linspace(R1, R2, N_BINS + 1)
    centers = 0.5 * (bins[:-1] + bins[1:])

    u_bin = _bin_radially(r[inside], speed, bins)
    p_bin = _bin_radially(r[inside], data["p"][inside], bins)
    p_bin -= p_bin[~np.isnan(p_bin)][0]

    return centers, u_bin, p_bin


def analytical_profile(data: dict):
    """
    Analytical solution from the VTK fields, sorted by radius. It is smooth, so
    it is not binned (binning would leave gaps where the mesh is coarser than
    the bin width). Pressure is shifted to zero at the inner cylinder.
    """

    r = np.hypot(data["cx"], data["cy"])
    inside = (r >= R1) & (r <= R2)

    speed = np.linalg.norm(data["u_exact"][inside][:, :2], axis=1)

    order = np.argsort(r[inside])
    p = data["p_exact"][inside][order]

    return r[inside][order], speed[order], p - p[0]


# =============================================================================
# Plot and table
# =============================================================================


def label(result: Result) -> str:
    return f"DuMux, {result.n_radial} × {result.n_angular} cells"


def generate_plot(results: list[Result], out_path: Path) -> None:

    fig, (ax_u, ax_p) = plt.subplots(1, 2, figsize=(12, 5), dpi=150, facecolor=SURFACE)

    r_ref, u_ref, p_ref = analytical_profile(results[0].data)

    for ax, ref, quantity, ylabel in (
        (ax_u, u_ref, 1, r"tangential velocity $u_\theta$ (m/s)"),
        (ax_p, p_ref, 2, "pressure $p$ (Pa)"),
    ):
        ax.set_facecolor(SURFACE)
        ax.plot(r_ref, ref, lw=1.8, color=INK, label="Analytical", zorder=1)

        for result, style in zip(results, MARKER_STYLES):
            r, u, p = numerical_profile(result.data)
            ax.plot(
                r, u if quantity == 1 else p, ls="none",
                label=label(result), zorder=2, **style,
            )

        ax.set_xlabel("radius $r$ (m)", color=INK_SECONDARY)
        ax.set_ylabel(ylabel, color=INK_SECONDARY)
        ax.tick_params(colors=INK_SECONDARY)
        ax.grid(True, color=INK, alpha=0.12, lw=0.6)
        for side in ("top", "right"):
            ax.spines[side].set_visible(False)
        for side in ("left", "bottom"):
            ax.spines[side].set_color(INK_SECONDARY)
        ax.legend(frameon=False, labelcolor=INK)

    ax_u.set_title("Velocity profile", color=INK)
    ax_p.set_title("Pressure profile", color=INK)

    fig.tight_layout()
    fig.savefig(out_path, bbox_inches="tight", facecolor=SURFACE)
    plt.close(fig)

    print(f"Saved: {out_path}")


README_START = "<!-- L2-ERRORS-START -->"
README_END = "<!-- L2-ERRORS-END -->"


def update_readme(table: str, figure: Path) -> None:
    """Insert the error table into README.md and copy the figure to images/."""

    readme = SCRIPT_DIR / "README.md"
    text = readme.read_text()

    if README_START not in text or README_END not in text:
        raise RuntimeError(f"Markers {README_START} / {README_END} not found in {readme}")

    head, rest = text.split(README_START, 1)
    _, tail = rest.split(README_END, 1)
    readme.write_text(f"{head}{README_START}\n{table}{README_END}{tail}")

    images = SCRIPT_DIR / "images"
    images.mkdir(exist_ok=True)
    shutil.copy(figure, images / figure.name)

    print(f"Updated: {readme} and {images / figure.name}")


def write_error_table(results: list[Result], out_path: Path) -> str:

    lines = [
        "| Grid (radial × angular cells) | Total cells | Rel. L2 error pressure | Rel. L2 error velocity |",
        "|:--|--:|--:|--:|",
    ]

    for result in results:
        lines.append(
            f"| {result.n_radial} × {result.n_angular} | {result.n_cells} "
            f"| {result.metrics['l2_error_pressure_rel']:.3e} "
            f"| {result.metrics['l2_error_velocity_rel']:.3e} |"
        )

    table = "\n".join(lines) + "\n"
    out_path.write_text(table)

    print(f"\n{table}")
    print(f"Saved: {out_path}")

    return table


# =============================================================================
# Main
# =============================================================================


def main():

    parser = argparse.ArgumentParser(
        description="Build, run (coarse and refined grid) and plot the Taylor-Couette benchmark."
    )
    parser.add_argument(
        "--build-dir",
        type=Path,
        default=None,
        help="DuMux build directory (default: build-cmake in the DuMux root)",
    )
    parser.add_argument(
        "--skip-build",
        action="store_true",
        help="do not build the executable before running it",
    )
    parser.add_argument(
        "--update-readme",
        action="store_true",
        help="write the error table into README.md and copy the figure to images/",
    )
    args = parser.parse_args()

    build_dir = find_build_dir(args.build_dir)
    if not args.skip_build:
        build(build_dir)

    exe_dir = build_dir / TEST_SUBDIR

    results = []
    params_path = SCRIPT_DIR / PARAMS_FILE
    base_radial, base_angular = grid_resolution(params_path)

    for problem_name, refinement in CASES:
        vtu_path, metrics = run_case(exe_dir, params_path, problem_name, refinement)
        data = load(vtu_path)
        # every global refinement halves the cell size in both directions
        n_radial = base_radial * 2**refinement
        n_angular = base_angular * 2**refinement
        results.append(Result(n_radial, n_angular, data["n_cells"], metrics, data))

    figure = exe_dir / "analytical_comparison.png"
    generate_plot(results, figure)
    table = write_error_table(results, exe_dir / "l2_errors.md")

    if args.update_readme:
        update_readme(table, figure)


if __name__ == "__main__":
    main()
