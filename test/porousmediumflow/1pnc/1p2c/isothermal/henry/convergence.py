#!/usr/bin/env python3
# SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later
"""
Convergence of the Henry problem: RMSE of the horizontal offset x_computed(z) - x_table(z) of the 10, 50 and 90 %
isochlors to the reference tables computed with henry_legendre.py, on a sequence of grids.

Run in the build directory:  python3 convergence.py
Runs test_1p2c_henry_diffusion_box with both input files on each grid (output in henry_convergence_runs) and
writes henry_convergence.png. Requires numpy, matplotlib and pyvista.

The test itself cannot be used for the error, since it needs grid lines at all heights z of the tables and fails
above its tolerance on the coarse grids. Its exit code is ignored here and only the output of the finished
simulation is used. Heights between grid lines are handled by linear interpolation in z.
"""

import csv
import os
import re
import shutil
import subprocess

import numpy as np
import pyvista as pv
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter

EXECUTABLE = os.path.abspath("test_1p2c_henry_diffusion_box")
WORK_DIR = "henry_convergence_runs"
DAY = 86400.0
# grids as cells in x and z, h = 3 m / cells in x; the end times are long enough for the isochlors to be steady at
# the scale of the error on the finest grid
TEST_CASES = [
    dict(name="diffusion", label="test case 1, molecular diffusion", paramFile="params_diffusion.input",
         table="reference_diffusion.csv", endTime=DAY, grids=[(60, 20), (120, 40), (240, 80), (480, 160)],
         style=dict(color="#4C72B0", marker="o", linestyle="-")),
    dict(name="dispersion", label="test case 2, velocity-dependent dispersion", paramFile="params_dispersion.input",
         table="reference_dispersion.csv", endTime=2*DAY, grids=[(60, 20), (120, 40), (240, 80), (360, 120)],
         style=dict(color="#9E3D22", marker="s", linestyle="--")),
]
LEVELS = [0.1, 0.5, 0.9]


def readTable(fileName):
    with open(fileName, newline="", encoding="utf-8") as f:
        return np.array([[float(v) for v in row] for row in csv.reader(line for line in f if line[:1].isdigit())])


def run(name, paramFile, table, cells, endTime):
    """Run the test on a grid and return the mesh of the last time step"""
    os.makedirs(WORK_DIR, exist_ok=True)
    subprocess.run([EXECUTABLE, os.path.abspath(paramFile), "-Grid.Cells", f"{cells[0]} {cells[1]}",
                    "-Problem.Name", name, "-Problem.ReferenceTable", os.path.abspath(table),
                    "-TimeLoop.TEnd", str(endTime)],
                   cwd=WORK_DIR, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=False)
    with open(os.path.join(WORK_DIR, name + ".pvd"), encoding="utf-8") as f:
        return pv.read(os.path.join(WORK_DIR, re.findall(r'file="([^"]+)"', f.read())[-1]))


def offsets(mesh, table):
    """x_computed(z) - x_table(z) for all heights z of the table and all isochlors"""
    points = np.round(mesh.points, 9)
    xs, zs = np.unique(points[:, 0]), np.unique(points[:, 1])
    c = np.zeros((len(zs), len(xs)))
    c[np.searchsorted(zs, points[:, 1]), np.searchsorted(xs, points[:, 0])] = np.ravel(mesh.point_data["X^Salt_liq"])
    result = []
    for z, *reference in table:
        line = np.array([np.interp(z, zs, c[:, i]) for i in range(len(xs))])
        row = []
        for level, x in zip(LEVELS, reference):  # first crossing from the freshwater side
            i = next(i for i in range(1, len(xs)) if line[i] >= level)
            row.append(xs[i - 1] + (level - line[i - 1])*(xs[i] - xs[i - 1])/(line[i] - line[i - 1]) - x)
        result.append(row)
    return np.array(result)


def slopeTriangle(ax, n0, n1, y0, order):
    """Right triangle below the corner (n0, y0) whose hypotenuse has the slope -order in the log-log plot"""
    y1 = y0*(n0/n1)**order
    ax.plot([n0, n1, n0, n0], [y0, y1, y1, y0], color="#4D4D4D", lw=1)
    ax.text(n0/1.06, np.sqrt(y0*y1), str(order), ha="right", va="center")


def plot(rmse):
    """Plot the RMSE in mm per test case, given as {name: list over the grids of the test case}"""
    if shutil.which("latex") and shutil.which("dvipng"):
        plt.rcParams.update({"text.usetex": True, "font.family": "serif", "font.serif": ["Computer Modern Roman"]})
    else:
        plt.rcParams.update({"font.family": "serif", "mathtext.fontset": "cm"})
    plt.rcParams.update({"font.size": 13, "legend.fontsize": 12, "xtick.labelsize": 12, "ytick.labelsize": 12})

    fig, ax = plt.subplots(figsize=(8.0, 4.2))
    for case in TEST_CASES:
        n = np.array([cells[0] for cells in case["grids"]])
        ax.loglog(n, rmse[case["name"]], label=case["label"], markersize=7, **case["style"])

    # below the last segment of test case 1
    (nA, _), (nB, _) = TEST_CASES[0]["grids"][-2:]
    n0, n1 = 1.2*nA, nB/1.2
    slopeTriangle(ax, n0, n1, 0.5*rmse[TEST_CASES[0]["name"]][-2]*(nA/n0)**2, 2)

    numCells = sorted({cells[0] for case in TEST_CASES for cells in case["grids"]})
    ax.set_xlabel(r"number of cells in $x$ direction $n_x$ ($h = 3\,\mathrm{m}/n_x$)")
    ax.set_ylabel("RMSE of the horizontal\nisochlor offset [mm]")
    ax.set_xticks(numCells)
    ax.set_xticklabels([str(n) for n in numCells])
    ax.minorticks_off()
    ax.yaxis.set_major_formatter(FuncFormatter(lambda value, _: f"{value:g}"))
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    ax.legend(frameon=False, loc="upper right")
    fig.tight_layout()
    fig.savefig("henry_convergence.png", dpi=200)
    print("wrote henry_convergence.png")


def main():
    rmse = {}
    for case in TEST_CASES:
        table = readTable(case["table"])
        rmse[case["name"]] = []
        for cells in case["grids"]:
            mesh = run(f"{case['name']}_{cells[0]}x{cells[1]}", case["paramFile"], case["table"], cells,
                       case["endTime"])
            error = offsets(mesh, table)
            rmse[case["name"]].append(1e3*np.sqrt(np.mean(error**2)))
            print(f"{case['label']}: {cells[0]}x{cells[1]} cells, RMSE = {rmse[case['name']][-1]:.3f} mm, "
                  f"maximum = {1e3*np.max(np.abs(error)):.2f} mm")
        n, e = np.array([cells[0] for cells in case["grids"]]), np.array(rmse[case["name"]])
        print(f"{case['label']}: observed orders", np.round(np.log(e[:-1]/e[1:])/np.log(n[1:]/n[:-1]), 2))
    plot(rmse)


if __name__ == "__main__":
    main()
