#!/usr/bin/env python3
# SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

"""
Growth benchmarks with exact solutions: runs the convergence studies of the strip bent by growth
and of the growing two-layer cylinder and sphere, prints the errors and rates of the benchmark
page and, if an output directory is given, writes the figures (errors against the mesh size and
the deformed shapes against the exact solutions). Runs in the build directory of the growth
tests, about three minutes with the studies in parallel.

usage: run_benchmark.py [output directory of the figures]
"""

import math
import os
import subprocess
import sys
import xml.etree.ElementTree as ET
from concurrent.futures import ThreadPoolExecutor

from convergencetest import FIELDS, run

STUDIES = [
    # (group, label, directory, executable, parameter file, cells of the coarsest mesh)
    ("strip", "box, 2D", "bending_strip", "test_hyperelastic_growth_bending_strip_box_2d",
     "params.input", [40, 4]),
    ("strip", "Q2 hybrid, 2D", "bending_strip",
     "test_hyperelastic_growth_bending_strip_pq2hybrid_2d", "params.input", [40, 4]),
    ("strip", "Q2 FE, 2D", "bending_strip", "test_hyperelastic_growth_bending_strip_pq2fe_2d",
     "params.input", [40, 4]),
    ("strip", "box, 3D", "bending_strip", "test_hyperelastic_growth_bending_strip_box_3d",
     "params_3d.input", [20, 2, 2]),
    ("strip", "MINI", "bending_strip", "test_hyperelastic_growth_bending_strip_mixed_p1b_2d",
     "params_incompressible.input", [40, 4]),
    ("strip", "Q2–Q1", "bending_strip", "test_hyperelastic_growth_bending_strip_mixed_q2_2d",
     "params_incompressible.input", [40, 4]),
    ("shell", "box, 2D", "layered_shell", "test_hyperelastic_growth_layered_shell_box_2d",
     "params.input", [6, 4, 16]),
    ("shell", "Q2 hybrid, 2D", "layered_shell",
     "test_hyperelastic_growth_layered_shell_pq2hybrid_2d", "params.input", [6, 4, 16]),
    ("shell", "box, 3D", "layered_shell", "test_hyperelastic_growth_layered_shell_box_3d",
     "params_3d.input", [3, 2, 4]),
    ("shell", "MINI", "layered_shell", "test_hyperelastic_growth_layered_shell_mixed_p1b_2d",
     "params_incompressible.input", [6, 4, 16]),
    ("shell", "Q2–Q1", "layered_shell", "test_hyperelastic_growth_layered_shell_mixed_q2_2d",
     "params_incompressible.input", [6, 4, 16]),
]
CURVATURE, LENGTH, THICKNESS = math.pi, 1.0, 0.1
CYLINDER_RADII = (0.601108, 0.891269, 1.138978)


def study(entry):
    _, _, directory, executable, params, cells = entry
    return [run(f"{directory}/{executable}", f"{directory}/{params}", [c * factor for c in cells])
            for factor in [1, 2, 4]]


def read_vtu(filename):
    piece = ET.parse(filename).getroot().find("UnstructuredGrid/Piece")

    def to_array(array):
        values = [float(v) for v in array.text.split()]
        n = int(array.get("NumberOfComponents", "1"))
        return [values[i:i + n] for i in range(0, len(values), n)]
    points = to_array(piece.find("Points/DataArray"))
    data = {array.get("Name"): to_array(array) for array in piece.find("PointData")}
    cells = piece.find("Cells")
    connectivity = [int(v) for v in cells.find("DataArray[@Name='connectivity']").text.split()]
    offsets = [int(v) for v in cells.find("DataArray[@Name='offsets']").text.split()]
    elements = [connectivity[b:e] for b, e in zip([0] + offsets[:-1], offsets)]
    return points, data, elements


def shapes(directory, executable, name, cells):
    subprocess.run([f"./{executable}", "params.input", "-Grid.Cells", " ".join(map(str, cells)),
                    "-Problem.EnableVtkOutput", "true", "-Problem.Name", name],
                   cwd=directory, check=True, stdout=subprocess.DEVNULL)
    vtu = sorted(f for f in os.listdir(directory) if f.startswith(name) and f.endswith(".vtu"))[-1]
    return read_vtu(os.path.join(directory, vtu))


def slope_triangle(ax, x0, y0, order, width=0.5):
    """a slope triangle with legs "1" and the order, below the point (x0, y0) in log-log axes"""
    x1 = x0 * 2**width
    y1 = y0 * 2**(width * order)
    ax.plot([x0, x1, x1, x0], [y0, y0, y1, y0], "k-", lw=1.0)
    ax.annotate("1", ((x0 * x1)**0.5, y0), textcoords="offset points", xytext=(0, -14),
                ha="center", fontsize=11)
    ax.annotate(f"{order}", (x1, (y0 * y1)**0.5), textcoords="offset points", xytext=(4, -4),
                fontsize=11)


def plot(results, outputDir):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import matplotlib.ticker
    import numpy as np

    plt.rcParams.update({"font.size": 13, "xtick.labelsize": 12, "ytick.labelsize": 12,
                         "legend.fontsize": 12})
    styles = {"box, 2D": ("o", "-", "#1f77b4"), "Q2 hybrid, 2D": ("s", "--", "#d62728"),
              "Q2 FE, 2D": ("^", ":", "#2ca02c"), "box, 3D": ("D", "-.", "#9467bd"),
              "MINI": ("v", (0, (5, 2, 1, 2, 1, 2)), "#8c564b"),
              "Q2–Q1": ("P", (0, (1, 1)), "#e377c2")}
    triangles = {("strip", 0): [(0.045, 2e-8, 4), (0.045, 1e-6, 2)],
                 ("strip", 1): [(0.045, 4e-6, 2), (0.045, 8e-5, 1)],
                 ("shell", 0): [(0.15, 4e-6, 2), (0.15, 3e-5, 1)],
                 ("shell", 1): [(0.15, 1e-4, 1), (0.15, 5e-4, 0.5)]}
    ticks = {"strip": [0.01, 0.02, 0.05, 0.1], "shell": [0.03, 0.1, 0.3]}
    fig, axes = plt.subplots(2, 2, figsize=(12, 9), layout="constrained")
    for row, group in enumerate(["strip", "shell"]):
        for col, field in enumerate(FIELDS):
            ax = axes[row, col]
            for entry, res in zip(STUDIES, results):
                if entry[0] != group:
                    continue
                marker, line, color = styles[entry[1]]
                ax.loglog([r["h"] for r in res], [r[field] for r in res], marker=marker, ls=line,
                          color=color, mfc="none", ms=7, label=entry[1])
            for x0, y0, order in triangles[(group, col)]:
                slope_triangle(ax, x0, y0, order)
            ax.set_ylim(bottom=ax.get_ylim()[0] / 2.5)
            ax.set_xticks(ticks[group])
            ax.set_xticklabels([f"{t:g}" for t in ticks[group]])
            ax.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
            ax.set_xlabel(r"$h_\mathrm{max}$")
            ax.set_ylabel(f"{'$L^2$' if col == 0 else '$H^1$'} error of the displacement")
            ax.set_title("strip" if group == "strip" else "cylinder and sphere")
    handles, labels = axes[0, 0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="outside lower center", ncol=6)
    fig.savefig(os.path.join(outputDir, "growth_convergence.png"), dpi=110)

    fig, axes = plt.subplots(1, 2, figsize=(12, 5.2), layout="constrained",
                             gridspec_kw={"width_ratios": [1.45, 1.0]})
    points, data, elements = shapes("bending_strip",
                                    "test_hyperelastic_growth_bending_strip_box_2d",
                                    "strip_shape", [80, 8])
    x = np.array(points)[:, :2] + np.array(data["d"])[:, :2]
    X = np.array(points)[:, :2]
    for element in elements:
        ax = axes[0]
        ax.fill(X[element, 0], X[element, 1], facecolor="none", edgecolor="0.65", lw=0.4)
        ax.fill(x[element, 0], x[element, 1], facecolor="none", edgecolor="k", lw=0.4)
    s = np.linspace(0.0, LENGTH, 200)
    for t in [-0.5 * THICKNESS, 0.5 * THICKNESS]:
        r = 1.0 / CURVATURE + t
        axes[0].plot(r * np.sin(CURVATURE * s), 0.5 * THICKNESS + r * np.cos(CURVATURE * s)
                     - 1.0 / CURVATURE, color="#d62728", lw=2.2, ls="--")
    axes[0].set_aspect("equal")
    axes[0].set_xlabel("$x$")
    axes[0].set_ylabel("$y$")
    axes[0].set_title("strip, $80 \\times 8$ cells")

    points, data, elements = shapes("layered_shell",
                                    "test_hyperelastic_growth_layered_shell_box_2d",
                                    "shell_shape", [3, 2, 8])
    x = np.array(points)[:, :2] + np.array(data["d"])[:, :2]
    X = np.array(points)[:, :2]
    for element in elements:
        axes[1].fill(X[element, 0], X[element, 1], facecolor="none", edgecolor="0.65", lw=0.6)
        axes[1].fill(x[element, 0], x[element, 1], facecolor="none", edgecolor="k", lw=0.6)
    phi = np.linspace(0.0, 0.5 * np.pi, 200)
    for radius in CYLINDER_RADII:
        axes[1].plot(radius * np.cos(phi), radius * np.sin(phi), color="#d62728", lw=2.2, ls="--")
    axes[1].set_aspect("equal")
    axes[1].set_xlabel("$x$")
    axes[1].set_title("cylinder, $3 \\times 2 \\times 8$ cells")
    from matplotlib.lines import Line2D
    fig.legend([Line2D([], [], color="0.65", lw=1.0), Line2D([], [], color="k", lw=1.0),
                Line2D([], [], color="#d62728", lw=2.2, ls="--")],
               ["reference mesh", "deformed mesh", "exact solution: boundaries, layer interface"],
               loc="outside lower center", ncol=3)
    fig.savefig(os.path.join(outputDir, "growth_shapes.png"), dpi=110)


def main():
    with ThreadPoolExecutor(max_workers=len(STUDIES)) as pool:
        results = list(pool.map(study, STUDIES))

    failed = False
    for entry, res in zip(STUDIES, results):
        print(f"\n{entry[0]}, {entry[1]} ({entry[3]})")
        print("| $h_\\mathrm{max}$ | $L^2$ error | rate | $H^1$ error | rate |")
        print("| ---: | ---: | ---: | ---: | ---: |")
        for i, r in enumerate(res):
            row = f"| {r['h']:.4f} |"
            for field in FIELDS:
                failed = failed or not math.isfinite(r[field])
                rate = ""
                if i > 0:
                    order = (math.log(res[i - 1][field] / r[field])
                             / math.log(res[i - 1]["h"] / r["h"]))
                    rate = f"{order:.2f}"
                row += f" {r[field]:.2e} | {rate} |"
            print(row)

    if len(sys.argv) > 1:
        plot(results, sys.argv[1])
    if failed:
        sys.exit("non-finite errors")


if __name__ == "__main__":
    main()
