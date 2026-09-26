#!/usr/bin/env python3
# SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later
"""Plot the VE solution of the radial injection benchmark against the similarity solution.

Usage (run from the build directory of the test, which contains params.input):
  python3 plot_convergence.py <executable>

Produces
- profile.png: interface height over the similarity variable at the end of the case of params.input and of
  two cases with the same injected volume at ten and a hundred times the injection rate,
- plume.png: fine-level saturation of the injected fluid at the end of the case of params.input,
- convergence.png: interface error over the radial cell size, refined together with the time step size.
"""

import configparser
import glob
import os
import subprocess
import sys
import xml.etree.ElementTree as ET

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import NullFormatter

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from convergencetest import (BASE_END_TIME, BASE_INITIAL_TIME_STEP, BASE_MAX_TIME_STEP, BASE_RATE,
                             compute_rates, grid_study, print_table, remove_outputs)


def read_parameters(file_name):
    parser = configparser.ConfigParser(inline_comment_prefixes=("#",))
    parser.optionxform = str
    parser.read(file_name)
    return parser


def read_vtu(file_name):
    """Return the cell centers and the cell data of an ASCII VTU file with quadrilateral cells."""
    piece = ET.parse(file_name).getroot().find(".//Piece")
    points = np.array(piece.find("Points/DataArray").text.split(), dtype=float).reshape(-1, 3)
    connectivity = np.array(piece.find("Cells/DataArray[@Name='connectivity']").text.split(), dtype=int)
    centers = points[connectivity.reshape(-1, 4)].mean(axis=1)
    cell_data = {array.get("Name"): np.array(array.text.split(), dtype=float) for array in piece.find("CellData")}
    return centers, cell_data


def read_pvd(file_name):
    """Return the times and file names of a PVD collection."""
    datasets = ET.parse(file_name).getroot().findall(".//DataSet")
    return [(float(d.get("timestep")), d.get("file")) for d in datasets]


def similarity_interface_height(chi, mobility_ratio):
    """Height of the interface above the bottom relative to the aquifer height, eq. (14) of Nordbotten and Celia (2006)."""
    thickness = (np.sqrt(2.0*mobility_ratio/chi) - 1.0)/(mobility_ratio - 1.0)
    return 1.0 - np.clip(thickness, 0.0, 1.0)


if len(sys.argv) < 2:
    sys.stderr.write("Usage: python3 plot_convergence.py <executable>\n")
    sys.exit(1)

testname = str(sys.argv[1])
exe = testname if os.path.isabs(testname) else "./" + testname
params = read_parameters("params.input")
name = params["Problem"]["Name"]

height = float(params["Grid"]["UpperRight"].split()[1]) - float(params["Grid"]["LowerLeft"].split()[1])
porosity = float(params["SpatialParams"]["Porosity"])
residual_saturation = float(params["SpatialParams"]["Swr"])
injection_rate = float(params["Problem"]["InjectionRate"])
mobility_ratio = float(params["1.Component"]["LiquidDynamicViscosity"])/float(params["2.Component"]["GasDynamicViscosity"])

# the same injected volume at increasing injection rates, which decreases the gravity number by the same factor
profile_rate_factors = [1, 10, 100]
profile_colors = ["#86b6ef", "#2a78d6", "#104281"]
fig, ax = plt.subplots(figsize=(6, 4))
chi_sqrt = np.linspace(1e-3, 5.0, 1000)
ax.plot(chi_sqrt, similarity_interface_height(chi_sqrt**2, mobility_ratio), color="black", lw=1.5,
        label=r"similarity solution, $\Gamma \to 0$")
for factor, color in zip(profile_rate_factors, profile_colors):
    case_name = name if factor == 1 else f"{name}_rate{factor}"
    output = subprocess.check_output([exe, "params.input",
                                      "-Problem.Name", case_name,
                                      "-Problem.InjectionRate", str(BASE_RATE*factor),
                                      "-TimeLoop.TEnd", str(BASE_END_TIME/factor),
                                      "-TimeLoop.MaxTimeStepSize", str(BASE_MAX_TIME_STEP/factor),
                                      "-TimeLoop.DtInitial", str(BASE_INITIAL_TIME_STEP/factor)], text=True)
    gravity_number = next(float(line.split()[-1]) for line in output.split("\n") if line.startswith("Mobility ratio: "))
    time, file_name = read_pvd(case_name + ".pvd")[-1]
    centers, cell_data = read_vtu(file_name)
    chi = 2.0*np.pi*height*porosity*(1.0 - residual_saturation)*centers[:, 0]**2/(injection_rate*factor*time)
    ax.plot(np.sqrt(chi), cell_data["interfaceHeight"]/height, color=color, lw=1.5,
            label=rf"VE model, $\Gamma = {gravity_number:.2g}$")
    if factor != 1:
        remove_outputs(case_name)
ax.set_xlim(0.0, 5.0)
ax.set_ylim(0.0, 1.05)
ax.set_xlabel(r"dimensionless distance $\chi^{1/2}$")
ax.set_ylabel("dimensionless interface height")
ax.legend(loc="lower right")
ax.grid(True, alpha=0.3)
fig.tight_layout()
fig.savefig("profile.png", dpi=150)
print("Saved profile.png")

# fine-level saturation of the injected fluid at the end of the simulation
time, file_name = read_pvd("fine_" + name + ".pvd")[-1]
centers, cell_data = read_vtu(file_name)
radii, heights = np.unique(centers[:, 0]), np.unique(centers[:, 1])
order = np.lexsort((centers[:, 0], centers[:, 1]))
saturation = cell_data["S_gas"][order].reshape(len(heights), len(radii))
fig, ax = plt.subplots(figsize=(8, 2.5))
mesh = ax.pcolormesh(radii, heights, saturation, cmap="viridis", vmin=0.0, vmax=1.0, shading="nearest")
ax.set_xlim(0.0, 1.3*np.sqrt(2.0*mobility_ratio*injection_rate*time/(2.0*np.pi*height*porosity*(1.0 - residual_saturation))))
ax.set_xlabel("radial distance $r$ [m]")
ax.set_ylabel("height $z$ [m]")
fig.colorbar(mesh, ax=ax, label="saturation of the injected fluid")
fig.tight_layout()
fig.savefig("plume.png", dpi=150)
print("Saved plume.png")

for f in glob.glob(name + "*.vtu") + glob.glob("fine_" + name + "*.vtu") + [name + ".pvd", "fine_" + name + ".pvd"]:
    if os.path.exists(f):
        os.remove(f)

# interface error for decreasing radial cell size and time step size at a negligible gravity number and capillary fringe
cell_sizes, errors = grid_study(testname)
print_table("dr [m]", cell_sizes, errors, compute_rates(cell_sizes, errors))

fig, ax = plt.subplots(figsize=(6, 4))
ax.loglog(cell_sizes, errors, "o-", color="#2a78d6", lw=1.5, label="VE model")
reference = errors[0]*np.array(cell_sizes)/cell_sizes[0]
ax.loglog(cell_sizes, reference, "k--", lw=1.0, label=r"$\mathcal{O}(\Delta r)$")
ax.set_xlabel(r"radial cell size $\Delta r$ [m]")
ax.set_ylabel("relative interface error")
ax.set_title(r"$\Gamma = 1.4 \cdot 10^{-4}$, $p_e = 1$ Pa, $\Delta t_{max} \propto \Delta r$")
ax.legend()
ax.grid(True, which="both", alpha=0.3)
ax.xaxis.set_minor_formatter(NullFormatter())
ax.set_xticks(cell_sizes)
ax.set_xticklabels([f"{size:.3g}" for size in cell_sizes])
fig.tight_layout()
fig.savefig("convergence.png", dpi=150)
print("Saved convergence.png")
