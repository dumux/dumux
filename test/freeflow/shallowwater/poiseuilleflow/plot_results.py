#!/usr/bin/env python3
# SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later
"""Poiseuille flow in a channel: velocity profile and relative errors under grid refinement.

Runs the channel with lateral viscous friction (executable of this folder) and the film flow with
viscous bottom friction (executable of the folder vertical).

Usage:
  python3 plot_results.py <executable> <executable of the folder vertical> [output directory of the figure]
"""

import math
import os
import re
import shlex
import sys
import xml.etree.ElementTree as ET

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
from plotstyle import cell_ticks, plt, series, slope_triangle
sys.path.insert(0, os.path.join(os.path.dirname(os.path.realpath(__file__)), "../../../../bin/testing"))
from dumux_convergencetest import convergenceRates, runLevels

ALONG = [(10, 320), (20, 320), (40, 320), (80, 320)]
ACROSS = [(320, 10), (320, 20), (320, 40), (320, 80)]
PROFILES = [(320, 10), (320, 20), (320, 40)]
FILM_CELLS = [20, 40, 80, 160]
STEADY_STATE = "-TimeLoop.TEnd 1e6 -TimeLoop.MaxTimeStepSize 1e5 -TimeLoop.PrintoutTimes 1e6"


def channel_errors(executable, grids, direction):
    """Mesh sizes in the given direction and relative L2 errors of the velocity of the steady flow"""
    rows = runLevels(f"{shlex.quote(executable)} params.input {STEADY_STATE} -Vtk.Precision Float64 -Problem.MeshSizeDirection {direction}",
                     [f"{nx} {ny}" for nx, ny in grids], "-Problem.Name poiseuilleflow_{0}x{1} -Grid.Cells '{0} {1}'",
                     "poiseuilleflow_{0}x{1}-errors.csv")
    return [r["h"] for r in rows], [r["L2Rel(u)"] for r in rows]


def film_errors(executable, cells, args=""):
    """Mesh sizes and relative L2 errors of the velocity of the film flow"""
    rows = runLevels(f"{shlex.quote(executable)} params.input {args}", [str(n) for n in cells],
                     "-Problem.Name filmflow_{0} -Grid.Cells '{0} 2'", "filmflow_{0}-errors.csv")
    return [r["h"] for r in rows], [r["L2Rel(u)"] for r in rows]


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
    return dict(re.findall(r"^\s*(\w+)\s*=\s*([^#\n]+)", open(filename).read(), re.MULTILINE))


def middle_column(cells):
    """Cell centers and data of the column of cells closest to the middle of the channel"""
    centers, data = read_cells(f"poiseuilleflow_{cells[0]}x{cells[1]}-00001.vtu")
    column = np.abs(centers[:, 0] - centers[np.argmin(np.abs(centers[:, 0] - 200.0)), 0]) < 1e-9
    order = np.argsort(centers[column, 1])
    return centers[column][order], {name: values[column][order] for name, values in data.items()}


def table(title, grids, meshSizes, errors):
    print(f"\n{title}\n| grid | rel. L2 error u | rate |\n|--:|--:|--:|")
    r = convergenceRates(meshSizes, errors)
    for i, (cells, error) in enumerate(zip(grids, errors)):
        print(f"| {cells[0]} x {cells[1]} | {error:.3e} | {r[i - 1]:.2f} |" if i > 0
              else f"| {cells[0]} x {cells[1]} | {error:.3e} | - |")


def main():
    if len(sys.argv) < 3:
        sys.exit(__doc__)
    executable, filmExecutable = (os.path.abspath(e) for e in sys.argv[1:3])
    imageDir = sys.argv[3] if len(sys.argv) > 3 else "."

    values = parameters()
    g, slopeS, viscosity = float(values["Gravity"]), float(values["BedSlope"]), float(values["TurbulentViscosity"])
    depth = float(values["HBoundary"])
    lower, upper = (np.array(values[k].split(), dtype=float) for k in ("LowerLeft", "UpperRight"))
    length, width = upper - lower
    exact = lambda y: g*slopeS/(8.0*viscosity)*(width**2 - 4.0*y**2)
    uMax = exact(0.0)
    print(f"Maximum velocity {uMax:.4f} m/s")

    hAlong, along = channel_errors(executable, ALONG, 0)
    hAcross, across = channel_errors(executable, ACROSS, 1)
    table("Refinement along the channel", ALONG, hAlong, along)
    table("Refinement across the channel", ACROSS, hAcross, across)

    # the Godunov mass flux between the reconstructed depths h and h - S dx exceeds h u by (c - u) S dx / 2
    print("\nMean velocity deficit in the middle of the channel\n"
          "| grid | measured (m/s) | (c - u) S dx / (2 h0) averaged over the width (m/s) |\n|--:|--:|--:|")
    for cells in ALONG:
        centers, data = middle_column(cells)
        predicted = np.mean((math.sqrt(g*depth) - exact(centers[:, 1]))*slopeS*(length/cells[0]))/(2.0*depth)
        print(f"| {cells[0]} x {cells[1]} | {np.mean(data['velocityX'] - exact(centers[:, 1])):+.3e} | {-predicted:+.3e} |")

    # the film flow reads its own input file in the folder vertical
    cwd = os.getcwd()
    os.chdir("vertical")
    hFilm, filmErrors = film_errors(filmExecutable, FILM_CELLS)
    table("Film flow, refinement along the channel", [(n, 2) for n in FILM_CELLS], hFilm, filmErrors)

    filmValues = parameters()
    filmLength = float(filmValues["UpperRight"].split()[0]) - float(filmValues["LowerLeft"].split()[0])
    baseSlope, baseDepth = float(filmValues["BedSlope"]), float(filmValues["WaterDepth"])
    print(f"\nFilm flow, {FILM_CELLS[0]} cells along the channel\n| S | h_0 (m) | rel. L2 error u | S dx/(2 h_0) |\n|--:|--:|--:|--:|")
    for s, h0 in [(baseSlope, baseDepth), (10*baseSlope, baseDepth), (baseSlope, 2*baseDepth)]:
        _, (error,) = film_errors(filmExecutable, FILM_CELLS[:1], f"-Problem.BedSlope {s} -Problem.WaterDepth {h0}")
        print(f"| {s:g} | {h0:g} | {error:.3e} | {s*filmLength/FILM_CELLS[0]/(2.0*h0):.3e} |")
    os.chdir(cwd)

    fig, (axP, axE, axF) = plt.subplots(1, 3, figsize=(12, 4.0), layout="constrained")
    for i, cells in enumerate(PROFILES):
        centers, data = middle_column(cells)
        y = centers[:, 1]
        axP.plot(y, (data["velocityX"] - exact(y))/uMax, label=f"{cells[1]} cells across", **series(i, len(y)))
    axP.axhline(0.0, color="black", lw=0.8)
    axP.set_xlabel("$y$ (m)")
    axP.set_ylabel(r"$(u - u_\mathrm{exact}) / u_\mathrm{max}$ at $x = L/2$")
    axP.legend()

    axE.loglog([c[0] for c in ALONG], along, label=f"along, {ALONG[0][1]} cells across", **series(0))
    axE.loglog([c[1] for c in ACROSS], across, label=f"across, {ACROSS[0][0]} cells along", **series(1))
    cells = [c[0] for c in ALONG]
    slope_triangle(axE, cells[-2], cells[-1], 0.6*along[-2], 1)
    slope_triangle(axE, cells[-2], cells[-1], 1.8*across[-2], 2, above=True)
    cell_ticks(axE, cells)
    axE.set_xlabel("number of cells in the refined direction")
    axE.set_ylabel(r"relative $L^2$ error of $u$")
    axE.legend(loc="lower left")

    axF.loglog(FILM_CELLS, filmErrors, label="film flow", **series(0))
    slope_triangle(axF, FILM_CELLS[-2], FILM_CELLS[-1], 0.6*filmErrors[-2], 1)
    cell_ticks(axF, FILM_CELLS)
    axF.set_xlabel("cells along the channel")
    axF.set_ylabel(r"relative $L^2$ error of $u$")
    axF.legend()
    fig.savefig(f"{imageDir}/shallowwater_poiseuille.png")


if __name__ == "__main__":
    main()
