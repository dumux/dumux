#!/usr/bin/env python3
# SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later
"""
Grid convergence test for the wellbore heat transport benchmark, comparing the
coupling modes of the 1D wellbore and the 3D rock domain.

The rock grids at refinement 0 are read from params.input: a graded tensor-product
grid for the surface coupling, an equidistant grid for the kernel coupling.
The kernel coupling manager rejects refined bulk grids, so instead of using
Rock.Grid.Refinement the script doubles the cells of the kernel grid per refinement.

Needs the Python packages numpy, matplotlib and vtk.
"""
import configparser
import glob
import math
import os
import re
import shutil
import subprocess
import sys

import matplotlib.pyplot as plt
import numpy as np
import vtk
from vtk.util.numpy_support import vtk_to_numpy

from ogs_solution import OGS_DEPTH, OGS_OUTLET, OGS_PROFILE, OGS_TIME
from ramey import RameySolution

# ── Simulation variants (different couplings of the 1D and 3D domain) ─────────
VARIANTS = {
    "surface": "test_wellbore_heat_transport_surface",
    "kernel": "test_wellbore_heat_transport_kernel",
}
REFINEMENTS = [0, 1, 2, 3]

# ── Directories ───────────────────────────────────────────────────────────────
# this directory in the source tree
SOURCE_DIR = os.path.dirname(os.path.realpath(__file__))
# input file of the benchmark, in the source directory
INPUT_FILE = os.path.join(SOURCE_DIR, "params.input")
# directory for the runs and the figures, relative to the build directory
OUTPUT_DIR = "convergence"
# build directory in the module root as created by dunecontrol
BUILD_DIR_NAME = "build-cmake"

# ── Plot style ────────────────────────────────────────────────────────────────
# one hue per variant
VARIANT_COLORS = {"surface": "#1C4E8C", "kernel": "#E8590C"}
REFERENCE_COLOR = "#000000"
GUIDE_COLOR = "#9a9a9a"
OGS_COLOR = "#00A087"
LINE_WIDTH = 1.5


def find_build_dir():
    """The directory matching the source directory in the build directory of the module."""
    root = SOURCE_DIR
    while not os.path.exists(os.path.join(root, "dune.module")):
        parent = os.path.dirname(root)
        if parent == root:
            sys.exit("error: no dune.module found above the source directory")
        root = parent
    return os.path.join(root, BUILD_DIR_NAME, os.path.relpath(SOURCE_DIR, root))


def build_executables(build_dir, targets):
    """Build the executables with CMake in the build directory of the module."""
    cmake_dir = build_dir
    while not os.path.exists(os.path.join(cmake_dir, "CMakeCache.txt")):
        parent = os.path.dirname(cmake_dir)
        if parent == cmake_dir:
            sys.exit(f"error: '{build_dir}' is not inside a CMake build directory")
        cmake_dir = parent
    command = ["cmake", "--build", cmake_dir, "--target", *targets]
    print(f"[build] {' '.join(command)}", flush=True)
    if subprocess.run(command).returncode != 0:
        sys.exit("error: building the executables failed")


# ── Input file handling ───────────────────────────────────────────────────────
def read_input(path):
    """Read a DuMux input file into a nested dict {group: {key: value}}."""
    parser = configparser.ConfigParser(inline_comment_prefixes=("#",), strict=False)
    parser.optionxform = str  # keep the case of the keys
    with open(path) as f:
        parser.read_string("[__root__]\n" + f.read())
    return {section: dict(parser[section]) for section in parser.sections()}


def read_dgf(path):
    """Return (borehole length, inner radius, outer radius)."""
    vertices, elements = [], []
    block = None
    with open(path) as f:
        for raw in f:
            line = raw.split("%")[0].strip()
            if not line:
                continue
            upper = line.upper()
            if upper in ("DGF", "#"):
                block = None if upper == "#" else block
                continue
            if upper.split()[0] in ("VERTEX", "SIMPLEX", "CUBE", "BOUNDARYDOMAIN",
                                    "BOUNDARYSEGMENTS", "INTERVAL"):
                block = upper.split()[0]
                continue
            if line.lower().startswith("parameters"):
                continue
            if block == "VERTEX":
                vertices.append([float(v) for v in line.split()])
            elif block == "SIMPLEX":
                elements.append(line.split())

    lower = [min(coords) for coords in zip(*vertices)]
    upper = [max(coords) for coords in zip(*vertices)]
    length = math.dist(lower, upper)
    first = elements[0]
    r_inner, r_outer = float(first[2]), float(first[3])
    return length, r_inner, r_outer


def collect_setup(input_file):
    """Everything the script needs to know about the simulation setup."""
    params = read_input(input_file)
    input_dir = os.path.dirname(os.path.abspath(input_file))

    grid_file = params["Well.Grid"]["File"]
    if not os.path.isabs(grid_file):
        grid_file = os.path.normpath(os.path.join(input_dir, grid_file))
    length, r_inner, r_outer = read_dgf(grid_file)

    setup = {
        "grid_file": grid_file,
        "length": length,
        "r_inner": r_inner,
        "r_outer": r_outer,
        "t_end": float(params["TimeLoop"]["TEnd"]),
        "injection_rate": float(params["Well.Problem"]["InjectionRate"]),
        "injection_temperature": float(params["Well.Problem"]["InjectionTemperature"]) - 273.15,
        "formation_temperature": float(params["Rock.Problem"]["InitialTemperature"]) - 273.15,
        "lambda_pg": float(params["MixedDimension"]["PipeWallThermalConductivity"]),
        "rho_f": float(params["Component"]["LiquidDensity"]),
        "c_p_f": float(params["Component"]["LiquidHeatCapacity"]),
        "mu_f": float(params["Component"]["LiquidDynamicViscosity"]),
        "lambda_f": float(params["Component"]["LiquidThermalConductivity"]),
        "lambda_s": float(params["Component"]["SolidThermalConductivity"]),
        "rho_s": float(params["Component"]["SolidDensity"]),
        "c_p_s": float(params["Component"]["SolidHeatCapacity"]),
        # equidistant rock grid of the kernel coupling at refinement 0
        "kernel_cells": [int(count) for count in params["Rock.Grid"]["Cells"].split()],
    }
    return setup


# ── Rock grid ─────────────────────────────────────────────────────────────────
def rock_grid_args(variant, refinement, setup):
    """Command line arguments for the rock grid of one variant at one refinement."""
    if variant == "kernel":
        # the kernel coupling manager does not support refined bulk grids, so the cells
        # are doubled in every direction like Rock.Grid.Refinement would do it
        cells = (count * 2**refinement for count in setup["kernel_cells"])
        return ["-Rock.Grid.Cells", " ".join(str(count) for count in cells)]
    return ["-Rock.Grid.Refinement", str(refinement)]


# ── Running the simulation ────────────────────────────────────────────────────
def run_simulation(variant, executable, input_file, setup, refinement, run_dir):
    """Run one variant at one refinement in its own directory."""
    if os.path.exists(run_dir):
        shutil.rmtree(run_dir)
    os.makedirs(run_dir)

    command = [executable, os.path.abspath(input_file)]
    command += rock_grid_args(variant, refinement, setup)
    command += ["-Well.Grid.Refinement", str(refinement)]
    # the grid file is given relative to the input file, make it absolute
    command += ["-Well.Grid.File", setup["grid_file"]]
    # the 3D bulk output is not evaluated here and is very large on the fine grids
    command += ["-Vtk.EnableBulkOutput", "false"]

    print(f"\033[1;36m[{variant} refinement {refinement}]\033[0;36m {' '.join(command)}\033[0m",
          flush=True)
    if subprocess.run(command, cwd=run_dir).returncode != 0:
        sys.exit(f"error: the simulation of {variant} refinement {refinement} failed")


# ── Reading the simulation results ────────────────────────────────────────────
def load_outlet_temperature(run_dir):
    """Outlet temperature over time from the csv written by the 1D problem."""
    files = sorted(glob.glob(os.path.join(run_dir, "*.csv")))
    if not files:
        sys.exit(f"error: no csv output found in {run_dir}")
    raw = np.genfromtxt(files[0], delimiter=",", skip_header=1)
    if raw.ndim == 1:
        raw = raw[np.newaxis, :]
    return raw[:, 0] * 86400.0, raw[:, 1]  # time in s, temperature in °C


def load_profile(run_dir):
    """Temperature profile along the borehole at the final time from the 1D vtp."""
    files = [f for f in glob.glob(os.path.join(run_dir, "*.vtp"))
             if re.search(r"-(\d+)\.vtp$", f)]
    if not files:
        sys.exit(f"error: no vtp output found in {run_dir}")
    vtp = max(files, key=lambda f: int(re.search(r"-(\d+)\.vtp$", f).group(1)))

    reader = vtk.vtkXMLPolyDataReader()
    reader.SetFileName(vtp)
    reader.Update()
    poly = reader.GetOutput()
    # the wellbore is vertical with the inlet at z = 0, so the distance along
    # the borehole is -z
    x = -vtk_to_numpy(poly.GetPoints().GetData())[:, 2]
    T = vtk_to_numpy(poly.GetPointData().GetArray("T")) - 273.15  # K -> °C
    order = np.argsort(x)
    return x[order], T[order]


# ── Results of one run ────────────────────────────────────────────────────────
def load_results(run_dir):
    """Outlet temperature over time and temperature profile of one run."""
    return {"outlet_data": load_outlet_temperature(run_dir),
            "profile_data": load_profile(run_dir)}


# ── Output ────────────────────────────────────────────────────────────────────
FIGURE_SIZE = (7, 5.5)
SOLUTION_FIGURE_SIZE = (7, 5.5)  # temperature panel above an error panel
FIGURE_DPI = 300


def save_figure(fig, path):
    fig.tight_layout()
    fig.savefig(path, dpi=FIGURE_DPI)
    print(f"Saved: {path}")


def solution_figure():
    """Figure with the temperatures on top and their errors vs Ramey below."""
    fig, (ax, error_ax) = plt.subplots(2, 1, sharex=True, figsize=SOLUTION_FIGURE_SIZE,
                                       gridspec_kw={"height_ratios": (3, 2)})
    return fig, ax, error_ax


def finish_solution_figure(fig, ax, error_ax, path, **legend_args):
    """Common styling of the figures with solutions and errors, then save."""
    ax.grid(color=GUIDE_COLOR, linewidth=0.8, alpha=0.7)
    error_ax.grid(color=GUIDE_COLOR, linewidth=0.8, alpha=0.7)
    error_ax.set_ylabel("Error vs Ramey [°C]")
    ax.legend(frameon=False, fontsize=9, **legend_args)
    save_figure(fig, path)


def plot_outlet(variants, refinement, results, setup, ramey, path):
    """Outlet temperature over time on the finest grid."""
    fig, ax, error_ax = solution_figure()

    # The error over time leaves out the first output: there the fluid that was in the pipe
    # initially is still flowing out, which Ramey neglects, and its error of about 1 °C would
    # set the scale of the error axis. The outlet temperature itself shows all outputs.
    def plot_error(t, T, color, linestyle="-"):
        error_ax.plot(t[1:] / 86400, T[1:] - ramey.outlet(t[1:]), color=color,
                      linewidth=LINE_WIDTH, linestyle=linestyle)

    times = np.linspace(0.0, setup["t_end"], 200)[1:]
    ax.plot(times / 86400, ramey.outlet(times), color=REFERENCE_COLOR, linewidth=LINE_WIDTH,
            linestyle=":", label="Ramey (1962)")
    # t = 0: the pipe still holds the cold initial water, which Ramey neglects
    mask = OGS_TIME > 0.0
    ax.plot(OGS_TIME[mask] / 86400, OGS_OUTLET[mask], color=OGS_COLOR,
            linewidth=LINE_WIDTH, label="OGS")
    plot_error(OGS_TIME[mask], OGS_OUTLET[mask], OGS_COLOR)
    for variant in variants:
        t, T = results[variant, refinement]["outlet_data"]
        color = VARIANT_COLORS[variant]
        # the kernel coupling is drawn dashed on top of the surface coupling, so both stay
        # visible where the solutions coincide
        linestyle = "--" if variant == "kernel" else "-"
        ax.plot(t / 86400, T, color=color, linewidth=LINE_WIDTH, linestyle=linestyle,
                label=f"{variant} refinement {refinement}")
        plot_error(t, T, color, linestyle)
    error_ax.set_xlabel("Time [d]")
    ax.set_ylabel("Outlet temperature [°C]")
    ax.set_title("Outlet temperature over time")

    finish_solution_figure(fig, ax, error_ax, path)


def plot_profile(variants, refinement, results, setup, ramey, path):
    """Temperature profile along the borehole at the final time on the finest grid."""
    fig, ax, error_ax = solution_figure()

    depths = np.linspace(0.0, setup["length"], 200)
    ax.plot(depths, ramey(depths, setup["t_end"]), color=REFERENCE_COLOR, linewidth=LINE_WIDTH,
            linestyle=":", label="Ramey (1962)")
    ax.plot(OGS_DEPTH, OGS_PROFILE, color=OGS_COLOR, linewidth=LINE_WIDTH, label="OGS")
    error_ax.plot(OGS_DEPTH, OGS_PROFILE - ramey(OGS_DEPTH, setup["t_end"]),
                  color=OGS_COLOR, linewidth=LINE_WIDTH)
    for variant in variants:
        x, T = results[variant, refinement]["profile_data"]
        color = VARIANT_COLORS[variant]
        # the kernel coupling is drawn dashed on top of the surface coupling, so both stay
        # visible where the solutions coincide
        linestyle = "--" if variant == "kernel" else "-"
        ax.plot(x, T, color=color, linewidth=LINE_WIDTH, linestyle=linestyle,
                label=f"{variant} refinement {refinement}")
        error_ax.plot(x, T - ramey(x, setup["t_end"]), color=color, linewidth=LINE_WIDTH,
                      linestyle=linestyle)
    error_ax.set_xlabel("Distance along borehole [m]")
    ax.set_ylabel("Fluid temperature [°C]")
    ax.set_xlim(0, setup["length"])
    ax.set_title(f"Temperature profile at t = {setup['t_end'] / 86400:.3g} d")

    # the curves rise to the upper right
    finish_solution_figure(fig, ax, error_ax, path, loc="upper left")


def refinement_changes(variant, refinements, results, setup):
    """Maximum change of the outlet temperature (over time) and of the profile (along the
    borehole) between successive refinements, as {quantity: (finer refinements, changes)}."""
    changes = {"outlet": ([], []), "profile": ([], [])}
    depths = np.linspace(0.0, setup["length"], 1001)
    for coarse, fine in zip(refinements[:-1], refinements[1:]):
        coarse_results, fine_results = results[variant, coarse], results[variant, fine]
        t, T = fine_results["outlet_data"]
        change = np.abs(T - np.interp(t, *coarse_results["outlet_data"])).max()
        changes["outlet"][0].append(fine)
        changes["outlet"][1].append(change)
        change = np.abs(np.interp(depths, *fine_results["profile_data"])
                        - np.interp(depths, *coarse_results["profile_data"])).max()
        changes["profile"][0].append(fine)
        changes["profile"][1].append(change)
    return changes


def plot_convergence(variants, refinements, results, setup, path):
    """Change of the solution with regard to the previous refinement. There is no exact solution
    (Ramey's is an approximation), so the convergence is measured between successive grids."""
    fig, ax = plt.subplots(figsize=FIGURE_SIZE)
    first_changes = []
    for variant in variants:
        changes = refinement_changes(variant, refinements, results, setup)
        for quantity, linestyle, marker in (("outlet", "-", "o"), ("profile", "--", "s")):
            fine, change = changes[quantity]
            ax.plot(fine, change, color=VARIANT_COLORS[variant], linewidth=LINE_WIDTH,
                    linestyle=linestyle, marker=marker, label=f"{variant}, {quantity}")
            first_changes.append(change[0])

    # each refinement halves the cell size, so second order convergence is a factor of 4
    steps = np.array(refinements[1:])
    ax.plot(steps, max(first_changes) * 4.0**-(steps - steps[0]), color=GUIDE_COLOR,
            linewidth=LINE_WIDTH, linestyle=":", label="2nd order")

    ax.set_yscale("log")
    ax.set_xticks(steps)
    ax.set_xlabel("Refinement")
    ax.set_ylabel("Max. change of the temperature [°C]")
    ax.set_title("Change with regard to the previous refinement")
    ax.grid(color=GUIDE_COLOR, linewidth=0.8, alpha=0.7)
    ax.legend(frameon=False, fontsize=9)
    save_figure(fig, path)


# ── Main ──────────────────────────────────────────────────────────────────────
def main():
    build_dir = find_build_dir()
    if not os.path.isdir(build_dir):
        sys.exit(f"error: the build directory '{build_dir}' does not exist")
    input_file = INPUT_FILE

    setup = collect_setup(input_file)
    ramey = RameySolution(
        injection_temperature=setup["injection_temperature"],
        formation_temperature=setup["formation_temperature"],
        injection_rate=setup["injection_rate"],
        length=setup["length"],
        r_inner=setup["r_inner"],
        r_outer=setup["r_outer"],
        lambda_pg=setup["lambda_pg"],
        rho_f=setup["rho_f"],
        c_p_f=setup["c_p_f"],
        mu_f=setup["mu_f"],
        lambda_f=setup["lambda_f"],
        lambda_s=setup["lambda_s"],
        rho_s=setup["rho_s"],
        c_p_s=setup["c_p_s"],
    )

    # every run gets its own directory below this one, next to the figures
    output_dir = os.path.join(build_dir, OUTPUT_DIR)
    os.makedirs(output_dir, exist_ok=True)

    refinements = REFINEMENTS
    variants = list(VARIANTS)
    run_dirs = {(variant, refinement):
                    os.path.join(output_dir, f"{variant}_refinement_{refinement}")
                for variant in variants for refinement in refinements}

    build_executables(build_dir, [VARIANTS[variant] for variant in variants])
    for (variant, refinement), run_dir in run_dirs.items():
        executable = os.path.join(build_dir, VARIANTS[variant])
        run_simulation(variant, executable, input_file, setup, refinement, run_dir)

    results = {key: load_results(run_dir) for key, run_dir in run_dirs.items()}

    # the solutions are shown on the finest grid, the coarser ones enter the convergence plot
    finest = refinements[-1]
    plot_outlet(variants, finest, results, setup, ramey,
                os.path.join(output_dir, "outlet_temperature.png"))
    plot_profile(variants, finest, results, setup, ramey,
                 os.path.join(output_dir, "temperature_profile.png"))
    plot_convergence(variants, refinements, results, setup,
                     os.path.join(output_dir, "refinement_convergence.png"))


if __name__ == "__main__":
    main()
