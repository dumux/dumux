#!/usr/bin/env python3
# SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later
"""
Convergence of the Henry problem: RMSE of the horizontal offset x_computed(z) - x_table(z) of the 10, 50 and 90 %
isochlors to the tables of Fahs et al. (2016), for grids with 9x3 to 240x80 cells.

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

EXECUTABLE = os.path.abspath("test_1p2c_henry_diffusion_box")
WORK_DIR = "henry_convergence_runs"
GRIDS = [(9, 3), (15, 5), (30, 10), (60, 20), (120, 40), (240, 80)]  # cells in x and z, h = 3 m / cells in x
# (name, title, input file, reference table, color, marker)
TEST_CASES = [
    ("diffusion", "Test Case 1 (purely diffusive)", "params_diffusion.input", "fahs2016_table_d1.csv", "#4C72B0", "o"),
    ("dispersion", r"Test Case 2 ($\alpha_L=0.1$ m, $\alpha_T=0.01$ m)", "params_dispersion.input", "fahs2016_table_d2.csv",
     "#9E3D22", "s"),
]
LEVELS = [0.1, 0.5, 0.9]


def readTable(fileName):
    with open(fileName, newline="", encoding="utf-8") as f:
        return np.array([[float(v) for v in row] for row in csv.reader(line for line in f if line[:1].isdigit())])


def run(name, paramFile, table, cells):
    """Run the test on a grid and return the mesh of the last time step"""
    os.makedirs(WORK_DIR, exist_ok=True)
    subprocess.run([EXECUTABLE, os.path.abspath(paramFile), "-Grid.Cells", f"{cells[0]} {cells[1]}",
                    "-Problem.Name", name, "-Problem.ReferenceTable", os.path.abspath(table)],
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


def main():
    if shutil.which("latex") and shutil.which("dvipng"):
        plt.rcParams.update({"text.usetex": True, "font.family": "serif",
                             "font.serif": ["Computer Modern Roman"], "font.size": 12})
    else:
        plt.rcParams.update({"font.family": "serif", "mathtext.fontset": "cm", "font.size": 12})

    fig, ax = plt.subplots(figsize=(8.0, 4.2))
    numCells = [cells[0] for cells in GRIDS]
    coarsest = []
    for name, title, paramFile, table, color, marker in TEST_CASES:
        rmse = []
        for cells in GRIDS:
            error = offsets(run(f"{name}_{cells[0]}x{cells[1]}", paramFile, table, cells), readTable(table))
            rmse.append(1e3*np.sqrt(np.mean(error**2)))
            print(f"{title}: {cells[0]}x{cells[1]} cells, RMSE = {rmse[-1]:.2f} mm")
        ax.loglog(numCells, rmse, marker + "-", color=color, label=title)
        coarsest.append(rmse[0])

    # the lines of the orders of convergence start at the coarsest grid: the first order at the larger of the two
    # points, the second order at the point of the second test case
    n = np.array([numCells[0], numCells[-1]], dtype=float)
    ax.loglog(n, max(coarsest)*(n/n[0])**-1, ":", color="#4D4D4D", lw=1, label=r"$\mathcal{O}(h)$")
    ax.loglog(n, coarsest[1]*(n/n[0])**-2, "--", color="#4D4D4D", lw=1, label=r"$\mathcal{O}(h^2)$")
    ax.set_xlabel(r"number of cells in $x$ direction $n_x$ ($h = 3\,\mathrm{m}/n_x$)")
    ax.set_ylabel("RMSE of the horizontal\nisochlor offset [mm]")
    ax.set_xticks(numCells)
    ax.set_xticklabels([str(n) for n in numCells])
    ax.minorticks_off()
    ax.legend(frameon=False, fontsize=9, loc="lower left")
    fig.tight_layout()
    fig.savefig("henry_convergence.png", dpi=200)
    print("wrote henry_convergence.png")


if __name__ == "__main__":
    main()
