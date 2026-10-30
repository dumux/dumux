#!/usr/bin/env python3
# SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later
"""
Build and run the wellbore heat transport benchmark on a sequence of grids and plot the results
together with the results of OpenGeoSys. Each run writes its temperatures and those of the reference
solution (analyticalsolution.hh), the errors are measured against the reference solution (b).
"""

import subprocess
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from ogs_solution import OGS_DEPTH, OGS_OUTLET, OGS_PROFILE, OGS_TIME


METHODS = ["surface", "kernel"]
REFINEMENTS = [0, 1]#, 2, 3]
KERNEL_CELLS = 68  # lateral cells of the kernel coupling grid at refinement 0, of params.input

ROOT_DIR = Path(__file__).resolve().parents[5]
BUILD_DIR = ROOT_DIR / "build-cmake/test/multidomain/embedded/1d3d/nonisothermal"


def build_and_run():
    targets = [f"test_wellbore_heat_transport_{method}" for method in METHODS]
    subprocess.run(["cmake", "--build", ROOT_DIR / "build-cmake", "--target", *targets], check=True)

    for method in METHODS:
        for refinement in REFINEMENTS:
            # the kernel coupling does not support refined rock grids, so its cells are doubled instead
            cells = KERNEL_CELLS*2**refinement
            rock_grid = (["-Rock.Grid.Cells", f"{cells} {cells} {2**refinement}"] if method == "kernel"
                         else ["-Rock.Grid.Refinement", str(refinement)])
            subprocess.run([f"./test_wellbore_heat_transport_{method}", "params.input", *rock_grid,
                            "-Well.Grid.Refinement", str(refinement),
                            "-Vtk.OutputName", f"{method}_refinement_{refinement}",
                            "-Vtk.EnableBulkOutput", "false"],
                           cwd=BUILD_DIR, check=True)


def read(method, refinement, kind):
    """The columns of an output file of one run."""
    return np.genfromtxt(BUILD_DIR / f"{method}_refinement_{refinement}_{kind}.csv",
                         delimiter=",", skip_header=1).T


def plot_solution(kind, reference, ogs_x, ogs_T, xlabel, ylabel, image_file):
    """The temperatures on the finest grid and below their errors, which span several orders of magnitude."""
    fig, (ax, error_ax) = plt.subplots(2, 1, sharex=True, figsize=(7, 5.5),
                                       gridspec_kw={"height_ratios": (3, 2)})

    x, reference_a, reference_b = reference
    ax.plot(x, reference_a, "k:", label="Reference (a) benchmark")
    ax.plot(x, reference_b, "k-.", label="Reference (b) embedded")
    error_ax.plot(x, reference_a - reference_b, "k:")

    ax.plot(ogs_x, ogs_T, "C2", label="OGS")
    error_ax.plot(ogs_x, ogs_T - np.interp(ogs_x, x, reference_b), "C2")

    # the kernel coupling is dashed, so both stay visible where the solutions coincide
    for method, style in zip(METHODS, ("C0-", "C1--")):
        x, T, _, reference_b = read(method, REFINEMENTS[-1], kind)
        ax.plot(x, T, style, label=f"{method} refinement {REFINEMENTS[-1]}")
        error_ax.plot(x, T - reference_b, style)

    error_ax.set_yscale("symlog", linthresh=0.01)
    error_ax.set_xlabel(xlabel)
    error_ax.set_ylabel("Error vs reference (b) [K]")
    ax.set_ylabel(ylabel)
    ax.grid(True, alpha=0.3)
    error_ax.grid(True, alpha=0.3)
    ax.legend()
    fig.tight_layout()
    fig.savefig(BUILD_DIR / image_file, dpi=300)


def plot_convergence(image_file):
    """The maximum errors over the refinement. The first day is left out for one of the outlet
    errors, as its much larger error would hide the late times."""
    fig, ax = plt.subplots(figsize=(7, 5.5))
    largest_error = 0.0
    for method, color in zip(METHODS, ("C0", "C1")):
        outlet, late_outlet, profile = [], [], []
        for refinement in REFINEMENTS:
            t, T, _, reference_b = read(method, refinement, "1d")
            outlet.append(np.abs(T - reference_b).max())
            late_outlet.append(np.abs(T - reference_b)[t >= 1.0].max())
            _, T, _, reference_b = read(method, refinement, "profile")
            profile.append(np.abs(T - reference_b).max())
        ax.plot(REFINEMENTS, outlet, "o-", color=color, label=f"{method}, outlet")
        ax.plot(REFINEMENTS, late_outlet, "s--", color=color, label=f"{method}, outlet from 1 d on")
        ax.plot(REFINEMENTS, profile, "D:", color=color, label=f"{method}, profile")
        largest_error = max(largest_error, outlet[0])

    # each refinement halves the cell size, so second order convergence is a factor of 4
    ax.plot(REFINEMENTS, largest_error*4.0**-np.array(REFINEMENTS), "k:", label="2nd order")
    ax.set_yscale("log")
    ax.set_xticks(REFINEMENTS)
    ax.set_xlabel("Refinement")
    ax.set_ylabel("Max. error vs reference (b) [K]")
    ax.grid(True, alpha=0.3)
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.12), ncol=3)
    fig.tight_layout()
    fig.savefig(BUILD_DIR / image_file, dpi=300)


def main():
    build_and_run()

    # t = 0 is left out, where the pipe still holds the cold initial water
    plot_solution("1d", read("surface", REFINEMENTS[-1], "outlet_reference"),
                  OGS_TIME[1:]/86400, OGS_OUTLET[1:],
                  "Time [d]", "Outlet temperature [°C]", "outlet_temperature.png")
    plot_solution("profile", read("surface", REFINEMENTS[-1], "profile_reference"),
                  OGS_DEPTH, OGS_PROFILE,
                  "Distance along borehole [m]", "Fluid temperature at the final time [°C]",
                  "temperature_profile.png")
    plot_convergence("refinement_convergence.png")


if __name__ == "__main__":
    main()
