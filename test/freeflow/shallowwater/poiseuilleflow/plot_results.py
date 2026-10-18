#!/usr/bin/env python3
# SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later
"""Poiseuille flow in a channel: velocity profile and errors under grid refinement.

Runs the channel with lateral viscous friction (executable of this folder) and the film flow with
viscous bottom friction (executable of the folder vertical).

Usage:
  python3 plot_results.py <executable> <executable of the folder vertical> [output directory of the figure]
"""

import importlib.util
import os
import re
import sys
import xml.etree.ElementTree as ET

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker
import numpy as np

from convergencetest import rates, run

ALONG = [(10, 320), (20, 320), (40, 320), (80, 320)]
ACROSS = [(320, 10), (320, 20), (320, 40), (320, 80)]
PROFILES = [(320, 10), (320, 20), (320, 40)]
FILM_CELLS = [20, 40, 80, 160]


def cell_ticks(ax, cells):
    """Label the logarithmic axis with the numbers of cells"""
    ax.set_xticks(cells, [str(c) for c in cells])
    ax.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())


def read_cells(filename):
    """Cell centers and cell data of an ASCII VTU file"""
    piece = ET.parse(filename).getroot().find("UnstructuredGrid/Piece")
    values = lambda array: np.array(array.text.split(), dtype=float)
    points = values(piece.find("Points/DataArray")).reshape(-1, 3)
    cells = {array.get("Name"): array for array in piece.find("Cells")}
    connectivity = values(cells["connectivity"]).astype(int)
    offsets = values(cells["offsets"]).astype(int)
    centers = np.array([points[connectivity[b:e]].mean(axis=0)
                        for b, e in zip(np.concatenate([[0], offsets[:-1]]), offsets)])
    return centers, {array.get("Name"): values(array) for array in piece.find("CellData")}


def parameters(filename="params.input"):
    """Gravity, bed slope, turbulent viscosity and channel width"""
    values = dict(re.findall(r"^\s*(\w+)\s*=\s*([^#\n]+)", open(filename).read(), re.MULTILINE))
    width = float(values["UpperRight"].split()[1]) - float(values["LowerLeft"].split()[1])
    return float(values["Gravity"]), float(values["BedSlope"]), float(values["TurbulentViscosity"]), width


def table(title, grids, errors):
    print(f"\n{title}\n| grid | L2 error u (m/s) | rate |\n|--:|--:|--:|")
    r = rates(errors)
    for i, (cells, error) in enumerate(zip(grids, errors)):
        print(f"| {cells[0]} x {cells[1]} | {error:.3e} | {r[i - 1]:.2f} |" if i > 0
              else f"| {cells[0]} x {cells[1]} | {error:.3e} | - |")


def main():
    if len(sys.argv) < 3:
        sys.exit(__doc__)
    executable, filmExecutable = (os.path.abspath(e) for e in sys.argv[1:3])
    imageDir = sys.argv[3] if len(sys.argv) > 3 else "."
    vtk = ["-Vtk.Precision", "Float64"]

    along = [run(executable, cells, vtk)["L2(u)"] for cells in ALONG]
    across = [run(executable, cells, vtk)["L2(u)"] for cells in ACROSS]
    table("Refinement along the channel", ALONG, along)
    table("Refinement across the channel", ACROSS, across)

    # the film flow reads its own input file in the folder vertical
    spec = importlib.util.spec_from_file_location("filmconvergencetest", "vertical/convergencetest.py")
    film = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(film)
    cwd = os.getcwd()
    os.chdir("vertical")
    filmErrors = [film.run(filmExecutable, cells)["L2(u)"] for cells in FILM_CELLS]
    table("Film flow, refinement along the channel", [(n, 2) for n in FILM_CELLS], filmErrors)

    # density and dynamic viscosity of the fluid system of the shallow water model
    density, viscosity = 1000.0, 1e-3
    values = dict(re.findall(r"^\s*(\w+)\s*=\s*([^#\n]+)", open("params.input").read(), re.MULTILINE))
    g, length = float(values["Gravity"]), float(values["UpperRight"].split()[0]) - float(values["LowerLeft"].split()[0])
    print(f"\nFilm flow, {FILM_CELLS[0]} cells along the channel\n| S | h_0 (m) | u (m/s) | relative error u | S dx/(2 h_0) |\n|--:|--:|--:|--:|--:|")
    for slope, depth in [(float(values["BedSlope"]), float(values["WaterDepth"])),
                         (10*float(values["BedSlope"]), float(values["WaterDepth"])),
                         (float(values["BedSlope"]), 2*float(values["WaterDepth"]))]:
        error = film.run(filmExecutable, FILM_CELLS[0], ["-Problem.BedSlope", str(slope), "-Problem.WaterDepth", str(depth)])["L2(u)"]
        u = density*g*slope*depth**2/(3.0*viscosity)
        print(f"| {slope:g} | {depth:g} | {u:.4e} | {error/u:.3e} | {slope*length/FILM_CELLS[0]/(2.0*depth):.3e} |")
    os.chdir(cwd)

    g, slope, viscosity, width = parameters()
    exact = lambda y: g*slope/(8.0*viscosity)*(width**2 - 4.0*y**2)
    print(f"\nMaximum velocity {exact(0.0):.4f} m/s")

    fig, (axP, axE, axF) = plt.subplots(1, 3, figsize=(12, 3.6), layout="constrained")
    for cells in PROFILES:
        centers, data = read_cells(f"poiseuilleflow_{cells[0]}x{cells[1]}-00001.vtu")
        column = np.abs(centers[:, 0] - centers[np.argmin(np.abs(centers[:, 0] - 200.0)), 0]) < 1e-9
        order = np.argsort(centers[column, 1])
        y = centers[column, 1][order]
        axP.plot(y, data["velocityX"][column][order] - exact(y), "o-", ms=3, lw=0.8, label=f"{cells[1]} cells across")
    axP.axhline(0.0, color="k", lw=0.8)
    axP.set_xlabel("$y$ (m)")
    axP.set_ylabel(r"$u - u_\mathrm{exact}$ at $x = L/2$ (m/s)")
    axP.legend(fontsize=8)

    axE.loglog([c[0] for c in ALONG], along, "o-", label=f"along, {ALONG[0][1]} cells across")
    axE.loglog([c[1] for c in ACROSS], across, "s-", label=f"across, {ACROSS[0][0]} cells along")
    n = np.array([10.0, 80.0])
    axE.loglog(n, along[0]*(n[0]/n), "k--", lw=0.8, label="slope 1")
    axE.loglog(n, across[0]*(n[0]/n)**2, "k:", lw=0.8, label="slope 2")
    cell_ticks(axE, [c[0] for c in ALONG])
    axE.set_xlabel("number of cells in the refined direction")
    axE.set_ylabel(r"$\|u - u_\mathrm{exact}\|_2$ (m/s)")
    axE.legend(fontsize=8)

    axF.loglog(FILM_CELLS, filmErrors, "o-", label="film flow")
    axF.loglog(FILM_CELLS, [filmErrors[0]*FILM_CELLS[0]/c for c in FILM_CELLS], "k--", lw=0.8, label="slope 1")
    cell_ticks(axF, FILM_CELLS)
    axF.set_xlabel("cells along the channel")
    axF.set_ylabel(r"$\|u - u_\mathrm{exact}\|_2$ (m/s)")
    axF.legend(fontsize=8)
    fig.savefig(f"{imageDir}/shallowwater_poiseuille.png", dpi=130)


if __name__ == "__main__":
    main()
