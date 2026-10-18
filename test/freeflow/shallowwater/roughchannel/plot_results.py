#!/usr/bin/env python3
# SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later
"""Flow in a rough channel: deviation from the normal flow and errors under grid refinement.

Usage:
  python3 plot_results.py <executable> [output directory of the figure]
"""

import re
import sys
import xml.etree.ElementTree as ET

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker
import numpy as np

from convergencetest import rates, run

CELLS = [250, 500, 1000, 2000]


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

    values = dict(re.findall(r"^\s*(\w+)\s*=\s*([^#\n]+)", open("params.input").read(), re.MULTILINE))
    n, q, slope = float(values["ManningN"]), abs(float(values["Discharge"])), float(values["BedSlope"])
    normalDepth = (n*q/slope**0.5)**0.6
    print(f"Normal flow: h = {normalDepth:.6f} m, u = {q/normalDepth:.6f} m/s\n")

    results = [run(executable, cells, ["-Vtk.Precision", "Float64"]) for cells in CELLS]
    rh, ru = rates([r["L1(h)"] for r in results]), rates([r["L1(u)"] for r in results])
    print("| cells | L1 error h (m) | rate | L1 error u (m/s) | rate |")
    print("|--:|--:|--:|--:|--:|")
    for i, (cells, r) in enumerate(zip(CELLS, results)):
        rates_ = (f"{rh[i - 1]:.2f}", f"{ru[i - 1]:.2f}") if i > 0 else ("-", "-")
        print(f"| {cells} | {r['L1(h)']:.3e} | {rates_[0]} | {r['L1(u)']:.3e} | {rates_[1]} |")

    print("\n| cells | min of h - h_n (m) | h - h_n in the outflow cell (m) |\n|--:|--:|--:|")
    fig, (axP, axE) = plt.subplots(1, 2, figsize=(10, 3.6), layout="constrained")
    for cells in CELLS:
        centers, data = read_cells(f"roughchannel_{cells}-00001.vtu")
        row = np.abs(centers[:, 1] - centers[:, 1].min()) < 1e-9
        order = np.argsort(centers[row, 0])
        deviation = data["waterDepth"][row][order] - normalDepth
        print(f"| {cells} | {deviation.min():.2e} | {deviation[-1]:.1e} |")
        axP.plot(centers[row, 0][order], deviation, lw=1.0, label=f"{cells} cells")
    axP.axhline(0.0, color="k", lw=0.8)
    axP.set_xlabel("$x$ (m)")
    axP.set_ylabel(r"$h - h_n$ (m)")
    axP.legend(fontsize=8)

    axE.loglog(CELLS, [r["L1(h)"] for r in results], "o-", label=r"$\|h - h_n\|_1$ (m)")
    axE.loglog(CELLS, [r["L1(u)"] for r in results], "s-", label=r"$\|u - u_n\|_1$ (m/s)")
    axE.loglog(CELLS, [results[0]["L1(u)"]*CELLS[0]/c for c in CELLS], "k--", lw=0.8, label="slope 1")
    cell_ticks(axE, CELLS)
    axE.set_xlabel("cells along the channel")
    axE.set_ylabel("$L^1$ error")
    axE.legend(fontsize=8)
    fig.savefig(f"{imageDir}/shallowwater_roughchannel.png", dpi=130)


if __name__ == "__main__":
    main()
