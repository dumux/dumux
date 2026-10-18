#!/usr/bin/env python3
# SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later
"""Lake at rest over a bed with a smooth hump and a discontinuous sill, submerged and partially dry.

Usage:
  python3 plot_results.py <executable> [output directory of the figure]
"""

import re
import subprocess
import sys
import xml.etree.ElementTree as ET

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

CASES = {
    "submerged": {"WaterLevel": 1.0, "MinWaterDepth": 0.0, "ErrorTolerance": 1e-10},
    "emerged": {"WaterLevel": 0.5, "MinWaterDepth": 1e-5, "ErrorTolerance": 1e-5},
}

DEVIATION = re.compile(r"free surface elevation: (\S+) m \(max\).*\n.*specific discharge:\s+(\S+) m\^2/s \(max\)")


def bed(x, y, hump, sill):
    return hump*np.exp(-5.0*(x - 0.9)**2 - 50.0*(y - 0.5)**2) + np.where((x > 1.4) & (x < 1.6), sill, 0.0)


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
    hump, sill = float(values["HumpHeight"]), float(values["SillHeight"])

    print("| case | water level H (m) | film thickness (m) | max deviation free surface (m) | max deviation discharge (m^2/s) |")
    print("|:--|--:|--:|--:|--:|")
    fig, (axS, axD) = plt.subplots(1, 2, figsize=(10, 3.6), layout="constrained")
    for case, params in CASES.items():
        name = f"lakeatrest_{case}"
        args = [a for key, value in params.items() for a in (f"-Problem.{key}", str(value))]
        output = subprocess.run(
            [executable, "params.input", "-Problem.Name", name, "-Vtk.Precision", "Float64", *args],
            check=True, capture_output=True, text=True,
        ).stdout
        surface, discharge = (float(v) for v in DEVIATION.search(output).groups())
        print(f"| {case} | {params['WaterLevel']} | {params['MinWaterDepth']:g} | {surface:.1e} | {discharge:.1e} |")

        centers, data = read_cells(f"{name}-00001.vtu")
        yRow = centers[np.argmin(np.abs(centers[:, 1] - 0.5)), 1]
        row = np.abs(centers[:, 1] - yRow) < 1e-9
        order = np.argsort(centers[row, 0])
        x = centers[row, 0][order]
        depth = data["waterDepth"][row][order]
        z = bed(x, centers[row, 1][order], hump, sill)
        exactDepth = np.maximum(params["WaterLevel"] - z, params["MinWaterDepth"])
        axS.plot(x, depth + z, lw=1.0, label=f"free surface, {case}")
        axD.semilogy(x, np.maximum(np.abs(depth - exactDepth), 1e-17), lw=1.0, label=case)

    xFine = np.linspace(0.0, 2.0, 2001)
    axS.fill_between(xFine, bed(xFine, yRow, hump, sill), color="0.8", label=f"bed at $y = {yRow:.2f}$ m")
    axS.set_xlabel("$x$ (m)")
    axS.set_ylabel("elevation (m)")
    axS.set_ylim(0.0, 1.2)
    axS.legend(fontsize=8, loc="upper left")
    axD.set_xlabel("$x$ (m)")
    axD.set_ylabel(r"$|h - h_\mathrm{exact}|$ (m)")
    axD.legend(fontsize=8)
    fig.savefig(f"{imageDir}/shallowwater_lakeatrest.png", dpi=130)


if __name__ == "__main__":
    main()
