#!/usr/bin/env python3
# SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later
"""Wet dam break: profiles at the final time and L1 errors under grid refinement.

The exact solution (Stoker) is computed here independently of the DuMux implementation.

Usage:
  python3 plot_results.py <executable> [output directory of the figure]
"""

import math
import sys
import xml.etree.ElementTree as ET

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker
import numpy as np
from scipy.optimize import brentq

from convergencetest import TIME_STEP_SIZE_100, rates, run

CELLS = [100, 200, 400, 800, 1600]
PROFILE_CELLS = [100, 400, 1600]
DEPTH_LEFT, DEPTH_RIGHT, GATE, GRAVITY, TIME = 4.0, 1.0, 10.0, 9.81, 1.0


def stoker(x, t):
    """Depth and velocity of the dam break on a wet bed, left rarefaction and right shock"""
    cl = math.sqrt(GRAVITY*DEPTH_LEFT)

    def mismatch(hm):
        shock = (hm - DEPTH_RIGHT)*math.sqrt(0.5*GRAVITY*(hm + DEPTH_RIGHT)/(hm*DEPTH_RIGHT))
        return 2.0*(cl - math.sqrt(GRAVITY*hm)) - shock

    hm = brentq(mismatch, DEPTH_RIGHT, DEPTH_LEFT, xtol=1e-14)
    um = 2.0*(cl - math.sqrt(GRAVITY*hm))
    shockSpeed = hm*um/(hm - DEPTH_RIGHT)
    xi = (np.asarray(x) - GATE)/t
    fan = (xi > -cl) & (xi <= um - math.sqrt(GRAVITY*hm))
    h = np.where(xi <= -cl, DEPTH_LEFT, np.where(fan, (2.0*cl - xi)**2/(9.0*GRAVITY),
                 np.where(xi <= shockSpeed, hm, DEPTH_RIGHT)))
    u = np.where(fan, 2.0/3.0*(cl + xi), np.where((xi > -cl) & (xi <= shockSpeed), um, 0.0))
    return h, u, hm, um, shockSpeed


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


def main():
    if len(sys.argv) < 2:
        sys.exit(__doc__)
    executable = sys.argv[1] if "/" in sys.argv[1] else "./" + sys.argv[1]
    imageDir = sys.argv[2] if len(sys.argv) > 2 else "."

    results = [run(executable, cells, ["-Vtk.Precision", "Float64"]) for cells in CELLS]
    _, _, hm, um, shockSpeed = stoker(0.0, TIME)
    print(f"Stoker: h_m = {hm:.6f} m, u_m = {um:.6f} m/s, shock speed {shockSpeed:.6f} m/s\n")

    rh = rates([r["L1(h)"] for r in results])
    ru = rates([r["L1(u)"] for r in results])
    print("| cells | time step size (s) | L1 error h (m) | rate | L1 error u (m/s) | rate |")
    print("|--:|--:|--:|--:|--:|--:|")
    for i, (cells, r) in enumerate(zip(CELLS, results)):
        rateH = f"{rh[i - 1]:.2f}" if i > 0 else "-"
        rateU = f"{ru[i - 1]:.2f}" if i > 0 else "-"
        print(f"| {cells} | {TIME_STEP_SIZE_100*100/cells:.3e} | {r['L1(h)']:.3e} | {rateH} "
              f"| {r['L1(u)']:.3e} | {rateU} |")

    # contributions to the L1 error of the water depth from the vicinity of the waves
    cl = math.sqrt(GRAVITY*DEPTH_LEFT)
    waves = {"bore": GATE + shockSpeed*TIME, "head of the rarefaction": GATE - cl*TIME,
             "tail of the rarefaction": GATE + (um - math.sqrt(GRAVITY*hm))*TIME}
    print("\n| cells | max deviation of exactWaterDepth (m) | "
          + " | ".join(f"within 1 m of the {w} (m)" for w in waves) + " | elsewhere (m) |")
    print("|--:|--:|" + "--:|"*(len(waves) + 1))
    for cells in CELLS:
        centers, data = read_cells(f"dambreak_{cells}-00001.vtu")
        h, *_ = stoker(centers[:, 0], TIME)
        contribution = np.abs(data["waterDepth"] - h)/len(h)
        near = {w: np.abs(centers[:, 0] - position) < 1.0 for w, position in waves.items()}
        elsewhere = ~np.any(list(near.values()), axis=0)
        print(f"| {cells} | {np.abs(data['exactWaterDepth'] - h).max():.1e} | "
              + " | ".join(f"{contribution[m].sum():.3e}" for m in near.values())
              + f" | {contribution[elsewhere].sum():.3e} |")

    fig, (axH, axU, axE) = plt.subplots(1, 3, figsize=(12, 3.6), layout="constrained")
    xFine = np.linspace(0.0, 20.0, 4001)
    hExact, uExact, *_ = stoker(xFine, TIME)
    axH.plot(xFine, hExact, "k", lw=1.0, label="exact")
    axU.plot(xFine, uExact, "k", lw=1.0, label="exact")
    for cells in PROFILE_CELLS:
        centers, data = read_cells(f"dambreak_{cells}-00001.vtu")
        order = np.argsort(centers[:, 0])
        axH.plot(centers[order, 0], data["waterDepth"][order], lw=1.0, label=f"{cells} cells")
        axU.plot(centers[order, 0], data["velocityX"][order], lw=1.0, label=f"{cells} cells")
    axH.set_xlabel("$x$ (m)")
    axH.set_ylabel("water depth $h$ (m)")
    axU.set_xlabel("$x$ (m)")
    axU.set_ylabel("velocity $u$ (m/s)")
    axH.legend(fontsize=8)

    axE.loglog(CELLS, [r["L1(h)"] for r in results], "o-", label=r"$\|h - h_\mathrm{exact}\|_1$ (m)")
    axE.loglog(CELLS, [r["L1(u)"] for r in results], "s-", label=r"$\|u - u_\mathrm{exact}\|_1$ (m/s)")
    axE.loglog(CELLS, [results[0]["L1(u)"]*CELLS[0]/n for n in CELLS], "k--", lw=0.8, label="slope 1")
    cell_ticks(axE, CELLS)
    axE.set_xlabel("cells along the channel")
    axE.set_ylabel("$L^1$ error")
    axE.legend(fontsize=8)
    fig.savefig(f"{imageDir}/shallowwater_dambreak.png", dpi=130)


if __name__ == "__main__":
    main()
