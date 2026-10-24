#!/usr/bin/env python3
# SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later
"""Plot the steady profiles and the grid convergence of the heat pipe test,
with full upwinding and with central weighting of the mobilities.

Usage (run from the test directory in the build tree):
  python3 plot_results.py <executable> [<executable> ...]
"""

import os
import sys
import subprocess
import numpy as np
import matplotlib.pyplot as plt

plt.rcParams.update({"font.size": 13, "axes.labelsize": 14, "xtick.labelsize": 12,
                     "ytick.labelsize": 12, "legend.fontsize": 12})

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from convergencetest import FIELDS, collect_data, compute_rates, print_table

if len(sys.argv) < 2:
    sys.stderr.write(__doc__)
    sys.exit(1)

executables = sys.argv[1:]
weightings = [("upwind", "1.0", "o"), ("central", "0.5", "s")]
schemeColors = {"box": "#eb6834", "tpfa": "#1baf7a"}
referenceColor = "#2a78d6"
inkColor = "#52514e"
columns = ["saturation", "temperature", "gas_pressure", "liquid_pressure"]
profileScale = {"saturation": (1.0, "-"), "temperature": (1.0, "K"),
                "gas_pressure": (1e-5, "$10^5$ Pa"), "liquid_pressure": (1e-5, "$10^5$ Pa")}
labels = {
    "saturation": ("liquid saturation", "-"),
    "temperature": ("temperature", "K"),
    "gas_pressure": ("gas pressure", "Pa"),
    "liquid_pressure": ("liquid pressure", "Pa"),
}


def scheme(executable):
    return os.path.basename(executable).rsplit("_", 1)[-1]


def slope_triangle(ax, h, error, order):
    """Draw a slope triangle of the given order with its right corner below (h, error) on log-log axes."""
    hLeft = h/4
    errorLeft = error/4**order
    ax.plot([hLeft, h, h, hLeft], [errorLeft, errorLeft, error, errorLeft], color=inkColor, linewidth=1.2, zorder=1)
    ax.text(np.sqrt(hLeft*h), errorLeft/1.5, "1", ha="center", va="top", color=inkColor)
    ax.text(h*1.12, np.sqrt(errorLeft*error), str(order), ha="left", va="center", color=inkColor)


def style(ax):
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    ax.grid(True, which="major", color="#e4e3df", linewidth=0.8)
    ax.tick_params(colors=inkColor)


configurations = [(k, executable, weighting, upwindWeight, schemeColors[scheme(executable)], marker)
                  for weighting, upwindWeight, marker in weightings
                  for k, executable in enumerate(executables)]

# grid convergence of the discrete L2 errors, the finest central run also provides the reference profile
results = []
for k, executable, weighting, upwindWeight, color, marker in configurations:
    name = f"{os.path.basename(executable)}-{weighting}"
    hs, errors = collect_data(executable, refinements=range(7), name=name,
                              extra_args=["-Flux.UpwindWeight", upwindWeight])
    rates = {field: compute_rates(hs, errors[field]) for field in FIELDS}
    print(f"\nConvergence rates for {executable}, {weighting} mobility weighting:")
    print_table(hs, errors, rates)
    results.append((hs, errors))
reference = np.genfromtxt(name + "-profile.csv", delimiter=",", names=True)
xRef, uniqueRef = np.unique(reference["x"], return_index=True)

fig, axes = plt.subplots(2, 2, figsize=(10, 8), sharex=True)
for (k, executable, weighting, _, color, marker), (hs, errors) in zip(configurations, results):
    for ax, field in zip(axes.flat, ["saturation", "temperature", "gas pressure", "liquid pressure"]):
        ax.loglog(hs, errors[field], marker=marker, markersize=8 - 3*k, linewidth=2 - 0.5*k, color=color,
                  markeredgecolor="white", markeredgewidth=1, label=f"{scheme(executable)}, {weighting}", zorder=2 + k)
for ax, field in zip(axes.flat, ["saturation", "temperature", "gas pressure", "liquid pressure"]):
    for (hs, errors), order, i in ((results[0], 1, 3), (results[-1], 2, 2)):
        slope_triangle(ax, hs[i], 0.4*errors[field][i], order)

for ax, column in zip(axes.flat, columns):
    name, unit = labels[column]
    ax.set_ylabel(f"$L^2$ error {name} ({unit})" if unit != "-" else f"$L^2$ error {name}", color=inkColor)
    style(ax)
for ax in axes[1]:
    ax.set_xlabel("mesh size h (m)", color=inkColor)
handles, legendLabels = axes[0, 0].get_legend_handles_labels()
fig.legend(handles, legendLabels, loc="upper center", ncol=len(legendLabels), frameon=False)
fig.tight_layout(rect=(0, 0, 1, 0.95))
fig.savefig("convergence.png", dpi=150)
print("Saved convergence.png")

# steady profiles on the coarsest grid
fig, axes = plt.subplots(2, 2, figsize=(10, 8), sharex=True)
for ax, column in zip(axes.flat, columns):
    ax.plot(xRef, profileScale[column][0]*reference[column + "_exact"][uniqueRef], color=referenceColor, linewidth=1.5, label="Udell (1985)", zorder=3)
for k, executable, weighting, upwindWeight, color, marker in configurations:
    exe = executable if os.path.isabs(executable) else "./" + executable
    name = f"{os.path.basename(executable)}-{weighting}-coarse"
    subprocess.run([exe, "params.input", "-Problem.Name", name, "-Grid.Refinement", "0",
                    "-Flux.UpwindWeight", upwindWeight], check=True, stdout=subprocess.DEVNULL)
    data = np.genfromtxt(name + "-profile.csv", delimiter=",", names=True)
    x, unique = np.unique(data["x"], return_index=True)
    for ax, column in zip(axes.flat, columns):
        ax.plot(x, profileScale[column][0]*data[column][unique], linestyle="none", marker=marker, markersize=6, color=color,
                markeredgecolor="white", markeredgewidth=1, label=f"{scheme(executable)}, {weighting}, {len(x)} dofs",
                zorder=2)

for ax, column in zip(axes.flat, columns):
    ax.set_ylabel(f"{labels[column][0]} ({profileScale[column][1]})", color=inkColor)
    style(ax)
for ax in axes[1]:
    ax.set_xlabel("distance from the cold end x (m)", color=inkColor)
handles, legendLabels = axes[0, 0].get_legend_handles_labels()
fig.legend(handles, legendLabels, loc="upper center", ncol=3, frameon=False)
fig.tight_layout(rect=(0, 0, 1, 0.91))
fig.savefig("profile.png", dpi=150)
print("Saved profile.png")
