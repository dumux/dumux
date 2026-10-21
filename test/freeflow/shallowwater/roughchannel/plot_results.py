#!/usr/bin/env python3
# SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later
"""Flow in a rough channel: deviation from the normal flow and relative errors under grid refinement.

With --second-order, the grids are also run with the executable with the second-order flux.

Usage:
  python3 plot_results.py <executable> [--second-order <executable>] [output directory of the figure]
"""

import math
import os
import re
import sys
import xml.etree.ElementTree as ET

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
from plotstyle import cell_ticks, plt, series, slope_triangle
from convergencetest import rates, run

CELLS = [250, 500, 1000, 2000]


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


def executable_path(name):
    return name if "/" in name else "./" + name


def main():
    args = sys.argv[1:]
    secondOrder = None
    if "--second-order" in args:
        i = args.index("--second-order")
        secondOrder = executable_path(args[i + 1])
        del args[i:i + 2]
    if not args:
        sys.exit(__doc__)
    executable = executable_path(args[0])
    imageDir = args[1] if len(args) > 1 else "."

    values = dict(re.findall(r"^\s*(\w+)\s*=\s*([^#\n]+)", open("params.input").read(), re.MULTILINE))
    n, q, slopeS, g = float(values["ManningN"]), abs(float(values["Discharge"])), float(values["BedSlope"]), float(values["Gravity"])
    length = float(values["Positions0"].split()[-1]) - float(values["Positions0"].split()[0])
    normalDepth = (n*q/slopeS**0.5)**0.6
    normalVelocity = q/normalDepth
    celerity = math.sqrt(g*normalDepth)
    print(f"Normal flow: h = {normalDepth:.6f} m, u = {normalVelocity:.6f} m/s\n")

    results = [run(executable, cells, ["-Vtk.Precision", "Float64"]) for cells in CELLS]
    rh, ru = rates([r["h"] for r in results]), rates([r["u"] for r in results])
    print("| cells | rel. L2 error h | rate | rel. L2 error u | rate |")
    print("|--:|--:|--:|--:|--:|")
    for i, (cells, r) in enumerate(zip(CELLS, results)):
        rate = (f"{rh[i - 1]:.2f}", f"{ru[i - 1]:.2f}") if i > 0 else ("-", "-")
        print(f"| {cells} | {r['h']:.3e} | {rate[0]} | {r['u']:.3e} | {rate[1]} |")

    # the Godunov mass flux between the reconstructed depths h and h - S dx exceeds h u by (c - u) S dx / 2
    print("\n| cells | min of (h - h_n)/h_n | (h - h_n)/h_n in the outflow cell "
          "| cell-centre discharge - q, 20 m < x < 120 m (m^2/s) | -(c - u_n) S dx / 2 (m^2/s) |")
    print("|--:|--:|--:|--:|--:|")
    fig, (axP, axE) = plt.subplots(1, 2, figsize=(10, 4.0), layout="constrained")
    for i, cells in enumerate(CELLS):
        centers, data = read_cells(f"roughchannel_{cells}-00001.vtu")
        row = np.abs(centers[:, 1] - centers[:, 1].min()) < 1e-9
        order = np.argsort(centers[row, 0])
        x = centers[row, 0][order]
        deviation = (data["waterDepth"][row][order] - normalDepth)/normalDepth
        discharge = (data["waterDepth"]*data["velocityX"])[row][order]
        upstream = (x > 20.0) & (x < 120.0)
        print(f"| {cells} | {deviation.min():.2e} | {deviation[-1]:.1e} | {discharge[upstream].mean() - q:+.4e} "
              f"| {-(celerity - normalVelocity)*slopeS*(length/cells)/2.0:+.4e} |")
        axP.plot(x, deviation, label=f"{cells} cells", **series(i, cells))
    axP.axhline(0.0, color="black", lw=0.8)
    axP.set_xlabel("$x$ (m)")
    axP.set_ylabel(r"$(h - h_n) / h_n$")
    axP.legend()

    axE.loglog(CELLS, [r["h"] for r in results], label=r"$h$", **series(0))
    axE.loglog(CELLS, [r["u"] for r in results], label=r"$u$", **series(1))
    slope_triangle(axE, CELLS[-2], CELLS[-1], 0.6*results[-2]["h"], 1)
    cell_ticks(axE, CELLS)
    axE.set_xlabel("cells along the channel")
    axE.set_ylabel(r"relative $L^2$ error at $t = 900$ s")
    axE.legend()
    fig.savefig(f"{imageDir}/shallowwater_roughchannel.png")

    if secondOrder:
        print("\nSecond-order flux:\n")
        print("| cells | rel. L2 error h | rel. L2 error u |")
        print("|--:|--:|--:|")
        for cells in CELLS:
            r = run(secondOrder, cells, ["-Vtk.Precision", "Float64"], name="roughchannel_secondorder")
            print(f"| {cells} | {r['h']:.1e} | {r['u']:.1e} |")


if __name__ == "__main__":
    main()
