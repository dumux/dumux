#!/usr/bin/env python3
# SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

"""
Horton-Rogers-Lapwood benchmark: growth rate of a convection cell under grid refinement
and as a function of the Rayleigh number, compared with linear stability theory.

Usage: python3 run_benchmark.py ./test_1p2c_hortonrogerslapwood_box [output directory of the figures]
"""

import math
import re
import subprocess
import sys
import xml.etree.ElementTree as ET

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

RESULT = re.compile(
    r"Ra = (\S+), k = (\S+), growth rate (\S+), exact (\S+), error (\S+), "
    r"max deviation from linear profile (\S+)"
)


def run(executable, name, *args):
    """Run the benchmark and return Ra, k, computed and exact growth rate, error and max deviation"""
    out = subprocess.run(
        [executable, "params.input", "-Problem.Name", name, *args],
        capture_output=True, text=True, check=True,
    ).stdout
    return [float(v) for v in RESULT.search(out).groups()]


def read_vtu(filename):
    """Return the points and the point data of an ASCII VTU file"""
    piece = ET.parse(filename).getroot().find("UnstructuredGrid/Piece")
    to_array = lambda array: np.array(array.text.split(), dtype=float).reshape(
        -1, int(array.get("NumberOfComponents", "1")))
    points = to_array(piece.find("Points/DataArray"))
    point_data = {array.get("Name"): to_array(array) for array in piece.find("PointData")}
    return points, point_data


def structured(points, values, cells):
    """Arrange vertex values of a structured grid with the given number of cells per direction"""
    lower, upper = points.min(axis=0), points.max(axis=0)
    i = np.rint((points[:, 0] - lower[0]) / (upper[0] - lower[0]) * cells[0]).astype(int)
    j = np.rint((points[:, 1] - lower[1]) / (upper[1] - lower[1]) * cells[1]).astype(int)
    field = np.empty((cells[1] + 1, cells[0] + 1))
    field[j, i] = values
    return field


def plot_cells(cases, filename):
    """Perturbation of the linear profile and streamlines of the Darcy velocity"""
    fig, axes = plt.subplots(1, len(cases), figsize=(10, 3.2), layout="constrained",
                             gridspec_kw={"width_ratios": [c["length"] for c in cases]})
    for ax, case in zip(axes, cases):
        points, data = read_vtu(case["vtu"])
        cells = case["cells"]
        x = np.linspace(0.0, case["length"], cells[0] + 1)
        z = np.linspace(0.0, 1.0, cells[1] + 1)
        perturbation = structured(points, data["X^Solute_liq"][:, 0] - points[:, 1], cells)
        perturbation /= np.abs(perturbation).max()
        u = structured(points, data["velocity_liq (m/s)"][:, 0], cells)
        w = structured(points, data["velocity_liq (m/s)"][:, 1], cells)
        image = ax.contourf(x, z, perturbation, levels=np.linspace(-1.0, 1.0, 21), cmap="RdBu_r")
        ax.streamplot(x, z, u, w, color="k", linewidth=0.6, density=0.8, arrowsize=0.8)
        ax.set_xlim(0.0, case["length"])
        ax.set_ylim(0.0, 1.0)
        ax.set_aspect("equal")
        ax.set_xlabel(r"$x / H$")
        ax.set_title(case["title"], fontsize=10)
    axes[0].set_ylabel(r"$z / H$")
    fig.colorbar(image, ax=axes, label=r"$(X - z/H) \,/\, \max |X - z/H|$", ticks=[-1, -0.5, 0, 0.5, 1])
    fig.savefig(filename, dpi=150)


def exact_growth_rate(ra, k):
    return ra * k**2 / (k**2 + math.pi**2) - (k**2 + math.pi**2)


def main():
    if len(sys.argv) < 2:
        sys.exit(__doc__)
    executable = sys.argv[1]
    image_dir = sys.argv[2] if len(sys.argv) > 2 else "."

    # the permeability is proportional to the Rayleigh number
    ra_input, *_ = run(executable, "hrl_input", "-TimeLoop.NumTimeSteps", "1")
    permeability_input = float(
        re.search(r"Permeability\s*=\s*(\S+)", open("params.input").read()).group(1)
    )

    def permeability(ra):
        return f"{permeability_input * ra / ra_input:.17g}"

    ra_critical = 4 * math.pi**2

    print("Grid convergence, Ra = 8 pi^2, k = pi")
    print("| Grid | growth rate | error | rate |")
    print("|------|------|------|------|")
    cells = [8, 16, 32, 64, 128]
    errors = []
    for n in cells:
        _, _, sigma, _, error, _ = run(
            executable, f"hrl_grid_{n}", "-Grid.Cells", f"{n} {n}",
            "-SpatialParams.Permeability", permeability(2 * ra_critical),
        )
        rate = f"{math.log(errors[-1] / error, 2):.2f}" if errors else "-"
        errors.append(error)
        print(f"| ${n} \\times {n}$ | {sigma:.5f} | ${error:.2e}$ | {rate} |")

    print("\nRayleigh number sweep, 64 cells per unit length")
    print("| Ra / Ra_c | k | growth rate | exact | error |")
    print("|------|------|------|------|------|")
    sweep = {math.pi: ([1, 1], [64, 64]), math.pi / 2: ([2, 1], [128, 64])}
    ratios = [0.5, 0.75, 1.0, 1.5, 2.0, 3.0, 4.0]
    computed = {}
    for k, (upper_right, grid) in sweep.items():
        computed[k] = []
        for ratio in ratios:
            ra, k_run, sigma, sigma_exact, error, _ = run(
                executable, "hrl_sweep", "-Grid.UpperRight", " ".join(map(str, upper_right)),
                "-Grid.Cells", " ".join(map(str, grid)),
                "-SpatialParams.Permeability", permeability(ratio * ra_critical),
            )
            computed[k].append(sigma)
            print(f"| {ratio:.2f} | {k_run:.4f} | {sigma:.4f} | {sigma_exact:.4f} | ${error:.1e}$ |")

    _, _, _, _, _, deviation = run(
        executable, "hrl_rest", "-Problem.Amplitude", "0",
        "-SpatialParams.Permeability", permeability(2 * ra_critical),
    )
    print(f"\nMaximum deviation from the linear profile at rest after 20 time steps: {deviation:.1e}")

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 4))
    ra_fine = [ra_critical * (0.4 + 3.7 * i / 200) for i in range(201)]
    for k, label, color in [(math.pi, r"$k = \pi$", "C0"), (math.pi / 2, r"$k = \pi/2$", "C1")]:
        ax1.plot([ra / ra_critical for ra in ra_fine], [exact_growth_rate(ra, k) for ra in ra_fine],
                 color=color, label=f"linear stability, {label}")
        ax1.plot(ratios, computed[k], "o", color=color, label=f"DuMux, {label}")
    ax1.axhline(0.0, color="0.5", lw=0.8)
    ax1.set_xlabel(r"$Ra / Ra_c$, $Ra_c = 4\pi^2$")
    ax1.set_ylabel(r"growth rate $\sigma$ $[D_\mathrm{eff} / (\phi H^2)]$")
    ax1.legend(fontsize=8)

    ax2.loglog(cells, errors, "o-", label=r"$|\sigma_h - \sigma| / (k^2 + \pi^2)$")
    ax2.loglog(cells, [errors[0] * (cells[0] / n) ** 2 for n in cells], "k--", lw=0.8, label="slope 2")
    ax2.set_xlabel("cells per unit length")
    ax2.set_ylabel("error of the growth rate")
    ax2.legend(fontsize=8)

    fig.tight_layout()
    fig.savefig(f"{image_dir}/hortonrogerslapwood_growthrate.png", dpi=150)

    cases = [
        {"name": "hrl_cell_k_pi", "length": 1, "cells": [64, 64], "title": r"$k = \pi$, $Ra = 2 Ra_c$"},
        {"name": "hrl_cell_k_pi_2", "length": 2, "cells": [128, 64], "title": r"$k = \pi/2$, $Ra = 2 Ra_c$"},
    ]
    for case in cases:
        run(executable, case["name"], "-Grid.UpperRight", f"{case['length']} 1",
            "-Grid.Cells", " ".join(map(str, case["cells"])), "-Vtk.Precision", "Float64",
            "-SpatialParams.Permeability", permeability(2 * ra_critical))
        case["vtu"] = f"{case['name']}-00001.vtu"
    plot_cells(cases, f"{image_dir}/hortonrogerslapwood_cells.png")


if __name__ == "__main__":
    main()
