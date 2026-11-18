#!/usr/bin/env python3
# SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

"""
DFG flow around a cylinder benchmark (Schäfer & Turek, 1996): plots of the grid convergence study

Reads the benchmark indicators of the runs of benchmark.py (runs/2d1/<scheme>_r<level>/*.json)
and plots the drag coefficient, the lift coefficient and the pressure difference over the total
number of unknowns, together with the FeatFlow reference values and the bounds of the paper.
The runs on the default meshes of the regression tests (<scheme>_default) are shown as unfilled
markers.
The figure is written to results/convergence_2d1.png (the errors and the orders of convergence
are in results/convergence_2d1.md, written by benchmark.py).

Usage (in the build directory test/freeflow/navierstokes/flow_around_cylinder):
    python3 visualize.py
"""

import argparse
import glob
import json
import os

import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from matplotlib.ticker import NullFormatter

from benchmark import REFERENCE, RESULTS_DIR, RUNS_DIR, SCHEMES

plt.rcParams.update({"mathtext.fontset": "cm", "font.family": "serif", "font.size": 11})

LABELS = {
    "diamond_tpfa": "Diamond / CCTpfa",
    "pq1bubble_box": "PQ1Bubble / Box",
    "pq1bubble_diamond": "PQ1Bubble / Diamond (quadrilaterals)",
}
MARKERS = {"diamond_tpfa": "^", "pq1bubble_box": "v", "pq1bubble_diamond": "s"}
SYMBOLS = {"cDrag": r"$C_D$", "cLift": r"$C_L$", "pDiff": r"$\Delta p$ [Pa]"}

# bounds of the paper for case 2D-1 (Schäfer & Turek, 1996, lower and upper bound of Table 3)
PAPER_BOUNDS = {"cDrag": (5.57, 5.59), "cLift": (0.0104, 0.0110), "pDiff": (0.1172, 0.1176)}
REFERENCE_STYLE = {"color": "k", "linestyle": "-", "linewidth": 1.0}
BOUNDS_STYLE = {"color": "0.6", "alpha": 0.4, "linewidth": 0.0}


def readRuns(scheme):
    """The results of the refined coarse meshes of a scheme, sorted by the number of unknowns"""
    runs = []
    for fileName in glob.glob(os.path.join(RUNS_DIR, f"{scheme}_r*", f"{scheme}_r*.json")):
        with open(fileName, encoding="utf-8") as jsonFile:
            runs.append(json.load(jsonFile))
    return sorted(runs, key=lambda run: run["numDofs"]["total"])


def plotReferenceValues(axes):
    """FeatFlow reference values (lines) and bounds of the paper (shaded) of the indicators"""
    for axis, key in zip(axes, SYMBOLS):
        axis.axhspan(*PAPER_BOUNDS[key], **BOUNDS_STYLE)
        axis.axhline(REFERENCE[key], **REFERENCE_STYLE)


def legendEntries(axis):
    """Legend entries of the schemes and the reference values"""
    handles, labels = axis.get_legend_handles_labels()
    defaultMarker = plt.Line2D(
        [], [], linestyle="none", marker="o", markersize=10, markerfacecolor="none", color="k"
    )
    handles += [defaultMarker, plt.Line2D([], [], **REFERENCE_STYLE), Patch(**BOUNDS_STYLE)]
    labels += [
        "default mesh (regression tests)",
        "FeatFlow reference value",
        "Schäfer & Turek (1996), bounds",
    ]
    return handles, labels


def readDefaultRun(scheme):
    """The result of a scheme on the default mesh of the regression tests, None if not available"""
    fileName = os.path.join(RUNS_DIR, f"{scheme}_default", f"{scheme}_default.json")
    if not os.path.exists(fileName):
        return None
    with open(fileName, encoding="utf-8") as jsonFile:
        return json.load(jsonFile)


def plotSchemes(axes):
    """
    Indicators of all schemes over the number of unknowns: lines for the refined coarse meshes,
    unfilled markers (same color) for the default meshes of the regression tests
    """
    for i, scheme in enumerate(SCHEMES):
        color = f"C{i}"
        runs = readRuns(scheme)
        default = readDefaultRun(scheme)
        for axis, key in zip(axes, SYMBOLS):
            if runs:
                dofs = [run["numDofs"]["total"] for run in runs]
                values = [run[key] for run in runs]
                axis.semilogx(
                    dofs, values, color=color, marker=MARKERS[scheme], label=LABELS[scheme]
                )
            if default is not None:
                axis.semilogx(
                    default["numDofs"]["total"],
                    default[key],
                    linestyle="none",
                    marker=MARKERS[scheme],
                    markersize=10,
                    markerfacecolor="none",
                    markeredgewidth=1.5,
                    color=color,
                )


def main():
    """Plot the benchmark indicators of case 2D-1 over the number of unknowns"""
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("--results-dir", default=RESULTS_DIR, help="output folder of the figure")
    args = parser.parse_args()

    figure, axes = plt.subplots(1, len(SYMBOLS), figsize=(13, 4.6))
    plotReferenceValues(axes)
    plotSchemes(axes)
    for axis, symbol in zip(axes, SYMBOLS.values()):
        axis.set_xlabel(r"number of unknowns $N$")
        axis.set_ylabel(symbol)
        axis.grid(True, which="both", linewidth=0.3)
        axis.xaxis.set_minor_formatter(NullFormatter())
    handles, labels = legendEntries(axes[0])
    figure.legend(handles, labels, loc="lower center", ncol=3, fontsize=9)
    figure.suptitle(r"Case 2D-1 ($\mathrm{Re} = 20$): grid convergence")
    figure.tight_layout(rect=(0.0, 0.12, 1.0, 1.0))

    os.makedirs(args.results_dir, exist_ok=True)
    fileName = os.path.join(args.results_dir, "convergence_2d1.png")
    figure.savefig(fileName, dpi=150)
    print(f"Wrote {fileName}")


if __name__ == "__main__":
    main()
