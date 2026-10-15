#!/usr/bin/env python3
# SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later
"""Plot the steady profiles and the grid convergence of the heat pipe test.

Usage (run from the test directory in the build tree):
  python3 plot_convergence.py <executable> [<executable> ...]
"""

import os
import sys
import subprocess
import numpy as np
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from convergencetest import FIELDS, collect_data, compute_rates, print_table

if len(sys.argv) < 2:
    sys.stderr.write(__doc__)
    sys.exit(1)

executables = sys.argv[1:]
colors = ["#eb6834", "#1baf7a"]
markers = ["o", "s"]
referenceColor = "#2a78d6"
inkColor = "#52514e"
columns = ["saturation", "temperature", "gas_pressure", "liquid_pressure"]
labels = {
    "saturation": ("liquid saturation", "-"),
    "temperature": ("temperature", "K"),
    "gas_pressure": ("gas pressure", "Pa"),
    "liquid_pressure": ("liquid pressure", "Pa"),
}


def scheme(executable):
    return os.path.basename(executable).rsplit("_", 1)[-1]


def style(ax):
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    ax.grid(True, which="major", color="#e4e3df", linewidth=0.8)
    ax.tick_params(colors=inkColor, labelsize=9)


# grid convergence of the discrete L2 errors, the finest run also provides the reference profile
results = []
for executable in executables:
    hs, errors = collect_data(executable, refinements=range(7))
    rates = {field: compute_rates(hs, errors[field]) for field in FIELDS}
    print(f"\nConvergence rates for {executable}:")
    print_table(hs, errors, rates)
    results.append((hs, errors))
reference = np.genfromtxt(os.path.basename(executables[0]) + "-profile.csv", delimiter=",", names=True)
xRef, uniqueRef = np.unique(reference["x"], return_index=True)

fig, axes = plt.subplots(2, 2, figsize=(9, 6), sharex=True)
for k, (executable, (hs, errors)) in enumerate(zip(executables, results)):
    for ax, field in zip(axes.flat, ["saturation", "temperature", "gas pressure", "liquid pressure"]):
        ax.loglog(hs, errors[field], marker=markers[k], markersize=8 - 3*k, linewidth=2 - 0.5*k, color=colors[k],
                  markeredgecolor="white", markeredgewidth=1, label=scheme(executable), zorder=2 + k)
for ax, field in zip(axes.flat, ["saturation", "temperature", "gas pressure", "liquid pressure"]):
    hs, errors = results[0]
    hRef = np.array([hs[0], hs[-1]])
    ax.loglog(hRef, 0.6*errors[field][0]*hRef/hs[0], "--", color=inkColor, linewidth=1, label=r"$\mathcal{O}(h)$", zorder=1)

for ax, column in zip(axes.flat, columns):
    name, unit = labels[column]
    ax.set_ylabel(f"$L^2$ error {name} [{unit}]" if unit != "-" else f"$L^2$ error {name}", color=inkColor)
    style(ax)
for ax in axes[1]:
    ax.set_xlabel("mesh size h [m]", color=inkColor)
axes[0, 0].legend(frameon=False, fontsize=9)
fig.tight_layout()
fig.savefig("convergence.png", dpi=150)
print("Saved convergence.png")

# steady profiles on the coarsest grid
fig, axes = plt.subplots(2, 2, figsize=(9, 6), sharex=True)
for ax, column in zip(axes.flat, columns):
    ax.plot(xRef, reference[column + "_exact"][uniqueRef], color=referenceColor, linewidth=2, label="Udell (1985)", zorder=1)
for k, executable in enumerate(executables):
    exe = executable if os.path.isabs(executable) else "./" + executable
    name = os.path.basename(executable) + "-coarse"
    subprocess.run([exe, "params.input", "-Problem.Name", name, "-Grid.Refinement", "0"],
                   check=True, stdout=subprocess.DEVNULL)
    data = np.genfromtxt(name + "-profile.csv", delimiter=",", names=True)
    x, unique = np.unique(data["x"], return_index=True)
    for ax, column in zip(axes.flat, columns):
        ax.plot(x, data[column][unique], linestyle="none", marker=markers[k], markersize=6, color=colors[k],
                markeredgecolor="white", markeredgewidth=1, label=f"{scheme(executable)}, {len(x)} dofs", zorder=2)

for ax, column in zip(axes.flat, columns):
    name, unit = labels[column]
    ax.set_ylabel(f"{name} [{unit}]", color=inkColor)
    style(ax)
for ax in axes[1]:
    ax.set_xlabel("distance from the cold end x [m]", color=inkColor)
axes[0, 0].legend(frameon=False, fontsize=9)
fig.tight_layout()
fig.savefig("profile.png", dpi=150)
print("Saved profile.png")
