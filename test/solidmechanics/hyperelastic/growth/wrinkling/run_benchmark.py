#!/usr/bin/env python3
# SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

"""
Onset of growth-induced wrinkling: critical growth stretch and wrinkled branch for the box and the
Q2 hybrid scheme on three nested meshes. Prints the tables of the benchmark page and, if an output
directory is given, writes the figure (wrinkled film, branch, convergence of the onset). Runs in
the build directory of the wrinkling test, about two minutes with the runs in parallel.

usage: run_benchmark.py [output directory of the figure]
"""

import math
import os
import subprocess
import sys
import xml.etree.ElementTree as ET
from concurrent.futures import ThreadPoolExecutor

from wrinklingtest import gridParameters

# critical growth stretch of the linearized bilayer problem (reference.py)
REFERENCE_STRETCH = 1.1430689262
SCHEMES = [("box", "test_hyperelastic_growth_wrinkling_box"),
           ("Q2 hybrid", "test_hyperelastic_growth_wrinkling_pq2hybrid")]
FACTORS = [1, 2, 4]
AMPLITUDES = [0.01, 0.02, 0.04, 0.08, 0.16]
SHAPE = ("box", 2)


def run(scheme, executable, factor):
    args = [f"./{executable}", "params.input",
            "-Problem.BranchAmplitudes", " ".join(map(str, AMPLITUDES)),
            "-Problem.Name", f"wrinkling_{executable[-3:]}_{factor}"]
    for key, value in gridParameters("params.input").items():
        if key.startswith("Cells"):
            args += [f"-Grid.{key}", " ".join(str(int(v) * factor) for v in value)]
        elif key.startswith("Grading"):
            args += [f"-Grid.{key}",
                     " ".join(str(math.copysign(abs(v) ** (1.0 / factor), v)) for v in value)]
    if (scheme, factor) == SHAPE:
        args += ["-Problem.EnableVtkOutput", "true"]
    output = subprocess.run(args, check=True, capture_output=True, text=True).stdout
    result = {"branch": []}
    for line in output.splitlines():
        if line.startswith("Branch amplitude "):
            result["branch"].append((float(line.split()[2].rstrip(":")), float(line.split()[-1])))
        elif ": " in line:
            key, value = line.split(": ", 1)
            try:
                result[key] = float(value.split()[0])
            except ValueError:
                pass
    return result


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
    return points, data, [connectivity[b:e] for b, e in zip([0] + offsets[:-1], offsets)]


def slope_triangle(ax, x0, y0, order, width=0.5):
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
    fig = plt.figure(figsize=(12, 8.2), layout="constrained")
    grid = fig.add_gridspec(2, 2, height_ratios=[1.0, 1.15])
    top = fig.add_subplot(grid[0, :])
    axBranch, axOnset = fig.add_subplot(grid[1, 0]), fig.add_subplot(grid[1, 1])

    scheme, factor = SHAPE
    vtu = sorted(f for f in os.listdir(".")
                 if f.startswith(f"wrinkling_box_{factor}") and f.endswith(".vtu"))[-1]
    points, data, elements = read_vtu(vtu)
    X = np.array(points)[:, :2]
    x = X + np.array(data["d"])[:, :2]
    filmBottom = X[:, 1].max() - 1.0
    for element in elements:
        if X[element, 1].min() < filmBottom - 3.0:
            continue
        film = X[element, 1].min() >= filmBottom - 1e-8
        top.fill(x[element, 0], x[element, 1], facecolor="0.8" if film else "none",
                 edgecolor="k", lw=0.4)
    top.set_aspect("equal")
    top.set_ylim(filmBottom - 3.0, filmBottom + 1.6)
    top.set_xlabel("$x$")
    top.set_ylabel("$y$")
    top.set_title(f"film (shaded) and substrate, $a = {AMPLITUDES[-1]}$, box, "
                  f"{factor * 18} cells along $x$")

    styles = [("o", "-", "#1f77b4"), ("s", "--", "#d62728"), ("^", ":", "#2ca02c")]
    for (marker, line, color), factor in zip(styles, FACTORS):
        for (name, _), filled in zip(SCHEMES, [False, True]):
            r = results[(name, factor)]
            a = np.array([p[0] for p in r["branch"]])
            stretch = np.array([p[1] for p in r["branch"]]) - r["Critical growth stretch"]
            label = f"{name}, $h_\\mathrm{{max}} = {r['Max element diameter']:.2f}$"
            axBranch.loglog(a, stretch, marker=marker, ls=line if not filled else "none",
                            color=color, mfc=color if filled else "none", ms=7, label=label)
    slope_triangle(axBranch, 0.03, 2e-5, 2)
    axBranch.set_xlabel("amplitude $a$")
    axBranch.set_ylabel(r"$\vartheta - \vartheta_{c,h}$")
    axBranch.set_ylim(bottom=axBranch.get_ylim()[0] / 2.5)

    for (name, _), (marker, line, color) in zip(SCHEMES, styles[:2]):
        h = [results[(name, f)]["Max element diameter"] for f in FACTORS]
        error = [abs(results[(name, f)]["Critical growth stretch"] - REFERENCE_STRETCH)
                 for f in FACTORS]
        axOnset.loglog(h, error, marker=marker, ls=line, color=color, mfc="none", ms=7,
                       label=name)
    slope_triangle(axOnset, 2.6, 1.5e-4, 2)
    meshSizes = [results[(SCHEMES[0][0], f)]["Max element diameter"] for f in FACTORS]
    axOnset.set_xticks(meshSizes)
    axOnset.set_xticklabels([f"{h:.2f}" for h in meshSizes])
    axOnset.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
    axOnset.set_xlabel(r"$h_\mathrm{max}$")
    axOnset.set_ylabel(r"$\vert \vartheta_{c,h} - \vartheta_c \vert$")
    axOnset.set_ylim(bottom=axOnset.get_ylim()[0] / 2.5)
    axOnset.legend(loc="upper left")
    axBranch.legend(loc="upper left", fontsize=10)
    fig.savefig(os.path.join(outputDir, "growth_wrinkling.png"), dpi=110)


def main():
    jobs = [(name, executable, factor) for name, executable in SCHEMES for factor in FACTORS]
    with ThreadPoolExecutor(max_workers=len(jobs)) as pool:
        values = list(pool.map(lambda job: run(*job), jobs))
    results = {(name, factor): value for (name, _, factor), value in zip(jobs, values)}

    failed = False
    for name, _ in SCHEMES:
        print(f"\nonset, {name}")
        print("| $h_\\mathrm{max}$ | critical stretch | error | rate | half wavelengths |")
        print("| ---: | ---: | ---: | ---: | ---: |")
        for i, factor in enumerate(FACTORS):
            r = results[(name, factor)]
            error = abs(r["Critical growth stretch"] - REFERENCE_STRETCH)
            failed = failed or not math.isfinite(error)
            rate = ""
            if i > 0:
                p = results[(name, FACTORS[i - 1])]
                order = (math.log(abs(p["Critical growth stretch"] - REFERENCE_STRETCH) / error)
                         / math.log(p["Max element diameter"] / r["Max element diameter"]))
                rate = f"{order:.2f}"
            print(f"| {r['Max element diameter']:.2f} | {r['Critical growth stretch']:.7f} | "
                  f"{error:.2e} | {rate} | {int(r['Wrinkle half-wavelengths'])} |")

    print("\nbranch: coefficient c2 and deviation of the intercept from the critical stretch")
    print("| scheme | " + " | ".join(f"mesh {f}" for f in FACTORS) + " |")
    print("| :--- |" + " ---: |" * len(FACTORS))
    for name, _ in SCHEMES:
        cells = [f"{results[(name, f)]['Branch coefficient c2']:.5f} "
                 f"({results[(name, f)]['Branch intercept deviation']:.1e})" for f in FACTORS]
        print(f"| {name} | " + " | ".join(cells) + " |")

    if len(sys.argv) > 1:
        plot(results, sys.argv[1])
    if failed:
        sys.exit("non-finite results")


if __name__ == "__main__":
    main()
