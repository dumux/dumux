#!/usr/bin/env python3
# SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later
"""
Grid convergence test for the wellbore heat transport benchmark, comparing the
coupling modes of the 1D wellbore and the 3D rock domain.

The rock grids at refinement 0 come from this script: a graded tensor-product grid
for the surface coupling, an equidistant grid for the kernel coupling.
The kernel coupling manager rejects refined bulk grids, so instead of using
Rock.Grid.Refinement the script doubles the cells of the kernel grid per refinement.

Options:
--variant: run only this coupling mode (default: both)
--refinements: which refinements to run (default: 0 1 2 3)
--build-dir: build directory of this test (default: the current directory, or the
             matching directory in build-cmake when started from the source directory)
--no-build: do not build the executables before running them
--no-run: only evaluate and plot existing results
--hide-ogs: do not plot the OGS results

Arguments that are not recognized are forwarded to the simulation, e.g.

    python3 convergence.py --refinements 0 1 -TimeLoop.TEnd 86400

Building and running the simulations needs the Python standard library only. Evaluating
and plotting the results needs numpy and matplotlib, the temperature profile also vtk.
Without numpy or matplotlib the simulations are still run, but not evaluated.
"""
import argparse
import configparser
import glob
import math
import os
import re
import shutil
import subprocess
import sys
import time

# only the evaluation needs these, building and running the simulations does not
try:
    import numpy as np
    import matplotlib.pyplot as plt

    from ogs_solution import OGS_DEPTH, OGS_OUTLET, OGS_PROFILE, OGS_TIME
    from ramey import RameySolution
except ImportError as error:
    MISSING_MODULE = error.name
else:
    MISSING_MODULE = None

# ── Rock grid at refinement 0 ─────────────────────────────────────────────────
# graded tensor-product grid of the surface coupling, per direction
BASE_CELLS = ("34 34", "34 34", "1")
GRADING = ("-1.2 1.2", "-1.2 1.2", "1")
# equidistant grid of the kernel coupling
KERNEL_BASE_CELLS = (68, 68, 1)

# ── Simulation variants (different couplings of the 1D and 3D domain) ─────────
VARIANTS = {
    "surface": "test_wellbore_heat_transport_surface",
    "kernel": "test_wellbore_heat_transport_kernel",
}

# ── Directories ───────────────────────────────────────────────────────────────
# this directory in the source tree (the build directory only holds a symlink to the script)
SOURCE_DIR = os.path.dirname(os.path.realpath(__file__))
# input file of the benchmark, in the source directory
INPUT_FILE = os.path.join(SOURCE_DIR, "params.input")
# directory for the runs and the figures, relative to the build directory
OUTPUT_DIR = "convergence"
# build directory in the module root as created by dunecontrol
BUILD_DIR_NAME = "build-cmake"

# ── Plot style ────────────────────────────────────────────────────────────────
# one hue per variant
VARIANT_COLORS = {"surface": "#0072B2", "kernel": "#D55E00"}
REFERENCE_COLOR = "#333333"
GUIDE_COLOR = "#9a9a9a"
OGS_COLOR = "#6A3D9A"
SOLUTION_WIDTH = 1.5  # line width of the temperatures
ERROR_WIDTH = 1.0     # line width of the errors vs Ramey (dashed)
GRID_STYLE = {"color": GUIDE_COLOR, "linewidth": 0.5, "alpha": 0.4}


def find_build_dir():
    """The current directory, or, when started from the source directory, the matching
    directory in the build directory of the module."""
    cwd = os.getcwd()
    if os.path.realpath(cwd) != SOURCE_DIR:
        return cwd
    root = SOURCE_DIR
    while not os.path.exists(os.path.join(root, "dune.module")):
        parent = os.path.dirname(root)
        if parent == root:
            sys.exit("error: no dune.module found above the source directory, "
                     "use --build-dir to point to the build directory")
        root = parent
    return os.path.join(root, BUILD_DIR_NAME, os.path.relpath(SOURCE_DIR, root))


def build_executables(build_dir, targets):
    """Build the executables with CMake in the build directory of the module."""
    cmake_dir = build_dir
    while not os.path.exists(os.path.join(cmake_dir, "CMakeCache.txt")):
        parent = os.path.dirname(cmake_dir)
        if parent == cmake_dir:
            sys.exit(f"error: '{build_dir}' is not inside a CMake build directory, "
                     "use --build-dir to point to it or --no-build to skip building")
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


def get(params, group, key, cast=str):
    """Read individual values from the input file."""
    try:
        value = params[group][key]
    except KeyError:
        sys.exit(f"error: [{group}] {key} is not set in the input file")
    try:
        return cast(value)
    except ValueError:
        sys.exit(f"error: [{group}] {key} = '{value}' is not a {cast.__name__}")


def read_dgf(path):
    """Return (borehole length, inner radius, outer radius, number of elements)."""
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
    return length, r_inner, r_outer, len(elements)


def collect_setup(input_file):
    """Everything the script needs to know about the simulation setup."""
    params = read_input(input_file)
    input_dir = os.path.dirname(os.path.abspath(input_file))

    grid_file = get(params, "Well.Grid", "File")
    if not os.path.isabs(grid_file):
        grid_file = os.path.normpath(os.path.join(input_dir, grid_file))
    length, r_inner, r_outer, num_segments = read_dgf(grid_file)

    setup = {
        "params": params,
        "grid_file": grid_file,
        "length": length,
        "r_inner": r_inner,
        "r_outer": r_outer,
        "num_segments": num_segments,
        "t_end": get(params, "TimeLoop", "TEnd", float),
        "injection_rate": get(params, "Well.Problem", "InjectionRate", float),
        "injection_temperature": get(params, "Well.Problem", "InjectionTemperature", float) - 273.15,
        "formation_temperature": get(params, "Rock.Problem", "InitialTemperature", float) - 273.15,
        "lambda_pg": get(params, "MixedDimension", "PipeWallThermalConductivity", float),
        "rho_f": get(params, "Component", "LiquidDensity", float),
        "c_p_f": get(params, "Component", "LiquidHeatCapacity", float),
        "mu_f": get(params, "Component", "LiquidDynamicViscosity", float),
        "lambda_f": get(params, "Component", "LiquidThermalConductivity", float),
        "lambda_s": get(params, "Component", "SolidThermalConductivity", float),
        "rho_s": get(params, "Component", "SolidDensity", float),
        "c_p_s": get(params, "Component", "SolidHeatCapacity", float),
    }
    return setup


# ── Rock grid ─────────────────────────────────────────────────────────────────
def rock_grid_args(variant, refinement):
    """Command line arguments for the rock grid of one variant at one refinement."""
    if variant == "kernel":
        # the kernel coupling manager does not support refined bulk grids, so the cells
        # are doubled in every direction like Rock.Grid.Refinement would do it
        cells = (count * 2**refinement for count in KERNEL_BASE_CELLS)
        return ["-Rock.Grid.Cells", " ".join(str(count) for count in cells)]

    args = []
    for direction, (cells, grading) in enumerate(zip(BASE_CELLS, GRADING)):
        args += [f"-Rock.Grid.Cells{direction}", cells]
        args += [f"-Rock.Grid.Grading{direction}", grading]
    return args + ["-Rock.Grid.Refinement", str(refinement)]


# ── Running the simulation ────────────────────────────────────────────────────
def run_simulation(variant, executable, input_file, setup, refinement, run_dir,
                   extra_args, tag):
    """Run one variant at one refinement in its own directory."""
    if os.path.exists(run_dir):
        shutil.rmtree(run_dir)
    os.makedirs(run_dir)

    command = [executable, os.path.abspath(input_file)]
    # the rock grid at refinement 0 comes from this script, not from params.input
    command += rock_grid_args(variant, refinement)
    command += ["-Well.Grid.Refinement", str(refinement)]
    # the grid file is given relative to the input file, make it absolute
    command += ["-Well.Grid.File", setup["grid_file"]]
    # the 3D bulk output is not evaluated here and is very large on the fine grids
    command += ["-Vtk.EnableBulkOutput", "false"]
    command += extra_args

    print(f"[{tag}] {' '.join(command)}", flush=True)
    log_file = os.path.join(run_dir, "run.log")
    start = time.time()
    # the DuMux output goes to the terminal as it appears and into run.log
    with subprocess.Popen(command, cwd=run_dir, stdout=subprocess.PIPE,
                          stderr=subprocess.STDOUT, text=True, bufsize=1) as process, \
            open(log_file, "w") as log:
        for line in process.stdout:
            sys.stdout.write(line)
            sys.stdout.flush()
            log.write(line)
    print(f"[{tag}] finished after {time.time() - start:.1f} s", flush=True)

    if process.returncode != 0:
        sys.exit(f"error: the simulation of {tag} failed, see {log_file}")


# ── Reading the simulation results ────────────────────────────────────────────
def load_outlet_temperature(run_dir):
    """Outlet temperature over time from the csv written by the 1D problem."""
    files = sorted(glob.glob(os.path.join(run_dir, "*.csv")))
    if not files:
        return None, None
    path = files[0]
    raw = np.genfromtxt(path, delimiter=",", skip_header=1)
    if raw.ndim == 1:
        raw = raw[np.newaxis, :]
    return raw[:, 0] * 86400.0, raw[:, 1]  # time in s, temperature in °C


def load_profile(run_dir):
    """Temperature profile along the borehole at the final time from the 1D vtp."""
    try:
        import vtk
        from vtk.util.numpy_support import vtk_to_numpy
    except ImportError:
        return None, None

    files = [f for f in glob.glob(os.path.join(run_dir, "*.vtp"))
             if re.search(r"-(\d+)\.vtp$", f)]
    if not files:
        return None, None
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
    results = {}

    times, temperatures = load_outlet_temperature(run_dir)
    if times is not None:
        results["outlet_data"] = (times, temperatures)

    depth, profile = load_profile(run_dir)
    if depth is not None:
        results["profile_data"] = (depth, profile)

    return results


# ── Output ────────────────────────────────────────────────────────────────────
FIGURE_SIZE = (7, 5.5)
FIGURE_DPI = 300


def save_figure(fig, path):
    fig.tight_layout()
    fig.savefig(path, dpi=FIGURE_DPI)
    print(f"Saved: {path}")


def finish_solution_figure(fig, ax, error_ax, path, **legend_args):
    """Common styling of the figures with solutions and errors, then save."""
    ax.grid(**GRID_STYLE)
    error_ax.set_ylabel("Absolute error vs Ramey (°C, dashed)")
    # the legend holds the solution curves only; it is put on the twin
    # because that one is drawn on top
    error_ax.legend(*ax.get_legend_handles_labels(), frameon=False, fontsize=9, **legend_args)
    save_figure(fig, path)


def plot_outlet(variants, refinement, results, setup, ramey, show_ogs, path):
    """Outlet temperature over time on the finest grid."""
    fig, ax = plt.subplots(figsize=FIGURE_SIZE)
    error_ax = ax.twinx()

    # The error over time leaves out the first output: there the fluid that was in the pipe
    # initially is still flowing out, which Ramey neglects, and its error of about 1 °C would
    # set the scale of the error axis. The outlet temperature itself shows all outputs.
    def plot_error(t, T, color):
        error_ax.plot(t[1:] / 86400, T[1:] - ramey.outlet(t[1:]), color=color,
                      linewidth=ERROR_WIDTH, linestyle="--")

    times = np.linspace(0.0, setup["t_end"], 200)[1:]
    ax.plot(times / 86400, ramey.outlet(times), color=REFERENCE_COLOR, linewidth=SOLUTION_WIDTH,
            linestyle="--", label="Ramey (1962)")
    if show_ogs:
        # t = 0: the pipe still holds the cold initial water, which Ramey neglects
        mask = OGS_TIME > 0.0
        ax.plot(OGS_TIME[mask] / 86400, OGS_OUTLET[mask], color=OGS_COLOR,
                linewidth=SOLUTION_WIDTH, label="OGS")
        plot_error(OGS_TIME[mask], OGS_OUTLET[mask], OGS_COLOR)
    for variant in variants:
        if "outlet_data" not in results[variant, refinement]:
            continue
        t, T = results[variant, refinement]["outlet_data"]
        color = VARIANT_COLORS[variant]
        ax.plot(t / 86400, T, color=color, linewidth=SOLUTION_WIDTH,
                label=f"{variant} ref {refinement}")
        plot_error(t, T, color)
    ax.set_xlabel("Time (d)")
    ax.set_ylabel("Outlet temperature (°C)")
    ax.set_title("Outlet temperature over time")
    # Both the temperatures and the errors end at the bottom of their axes, so all curves
    # would meet in the lower right corner. Lift the errors a bit above the temperatures.
    low, high = error_ax.get_ylim()
    error_ax.set_ylim(low - 0.3 * (high - low), high)

    finish_solution_figure(fig, ax, error_ax, path)
    return fig


def plot_profile(variants, refinement, results, setup, ramey, show_ogs, path):
    """Temperature profile along the borehole at the final time on the finest grid."""
    fig, ax = plt.subplots(figsize=FIGURE_SIZE)
    error_ax = ax.twinx()

    depths = np.linspace(0.0, setup["length"], 200)
    ax.plot(depths, ramey(depths, setup["t_end"]), color=REFERENCE_COLOR, linewidth=SOLUTION_WIDTH,
            linestyle="--", label="Ramey (1962)")
    if show_ogs:
        ax.plot(OGS_DEPTH, OGS_PROFILE, color=OGS_COLOR, linewidth=SOLUTION_WIDTH, label="OGS")
        error_ax.plot(OGS_DEPTH, OGS_PROFILE - ramey(OGS_DEPTH, setup["t_end"]),
                      color=OGS_COLOR, linewidth=ERROR_WIDTH, linestyle="--")
    for variant in variants:
        if "profile_data" not in results[variant, refinement]:
            continue
        x, T = results[variant, refinement]["profile_data"]
        color = VARIANT_COLORS[variant]
        ax.plot(x, T, color=color, linewidth=SOLUTION_WIDTH, label=f"{variant} ref {refinement}")
        error_ax.plot(x, T - ramey(x, setup["t_end"]), color=color, linewidth=ERROR_WIDTH,
                      linestyle="--")
    ax.set_xlabel("Distance along borehole (m)")
    ax.set_ylabel("Fluid temperature (°C)")
    ax.set_xlim(0, setup["length"])
    ax.set_title(f"Temperature profile at t = {setup['t_end'] / 86400:.3g} d")

    # the curves rise to the upper right and the errors start in the upper left
    finish_solution_figure(fig, ax, error_ax, path, loc="upper center", ncol=2)
    return fig


def refinement_changes(variant, refinements, results, setup):
    """Maximum change of the outlet temperature (over time) and of the profile (along the
    borehole) between successive refinements, as {quantity: (finer refinements, changes)}."""
    changes = {"outlet": ([], []), "profile": ([], [])}
    depths = np.linspace(0.0, setup["length"], 1001)
    for coarse, fine in zip(refinements[:-1], refinements[1:]):
        coarse_results, fine_results = results[variant, coarse], results[variant, fine]
        if "outlet_data" in coarse_results and "outlet_data" in fine_results:
            t, T = fine_results["outlet_data"]
            change = np.abs(T - np.interp(t, *coarse_results["outlet_data"])).max()
            changes["outlet"][0].append(fine)
            changes["outlet"][1].append(change)
        if "profile_data" in coarse_results and "profile_data" in fine_results:
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
            if not fine:
                continue
            ax.plot(fine, change, color=VARIANT_COLORS[variant], linewidth=SOLUTION_WIDTH,
                    linestyle=linestyle, marker=marker, label=f"{variant}, {quantity}")
            first_changes.append(change[0])

    # each refinement halves the cell size, so second order convergence is a factor of 4
    steps = np.array(refinements[1:])
    if first_changes:
        ax.plot(steps, max(first_changes) * 4.0**-(steps - steps[0]), color=GUIDE_COLOR,
                linewidth=ERROR_WIDTH, linestyle=":", label="2nd order")

    ax.set_yscale("log")
    ax.set_xticks(steps)
    ax.set_xlabel("Refinement")
    ax.set_ylabel("Max. change of the temperature (°C)")
    ax.set_title("Change with regard to the previous refinement")
    ax.grid(**GRID_STYLE)
    ax.legend(frameon=False, fontsize=9)
    save_figure(fig, path)
    return fig


# ── Main ──────────────────────────────────────────────────────────────────────
def main():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--refinements", type=int, nargs="+", default=[0, 1, 2, 3],
                        help="refinements to run (default: 0 1 2 3)")
    parser.add_argument("--variant", choices=list(VARIANTS),
                        help="run only this coupling mode (default: both)")
    parser.add_argument("--build-dir",
                        help="build directory of this test (default: the current directory, or "
                             f"the matching directory in {BUILD_DIR_NAME} when started from the "
                             "source directory)")
    parser.add_argument("--no-build", action="store_true",
                        help="do not build the executables before running them")
    parser.add_argument("--no-run", action="store_true",
                        help="only evaluate and plot existing results")
    parser.add_argument("--hide-ogs", action="store_true",
                        help="do not plot the OGS results")
    args, extra_args = parser.parse_known_args()
    # the simulation only takes -Group.Key value pairs, anything else is a mistyped option
    keys = extra_args[::2]
    if len(extra_args) % 2 or any(not k.startswith("-") or k.startswith("--") for k in keys):
        parser.error(f"unrecognized arguments: {' '.join(extra_args)}")

    build_dir = os.path.abspath(args.build_dir) if args.build_dir else find_build_dir()
    if not os.path.isdir(build_dir):
        sys.exit(f"error: the build directory '{build_dir}' does not exist, "
                 "use --build-dir to point to it")
    input_file = INPUT_FILE

    setup = collect_setup(input_file)
    if MISSING_MODULE is not None:
        if args.no_run:
            sys.exit(f"error: evaluating the results needs the Python module '{MISSING_MODULE}'")
        print(f"note: the Python module '{MISSING_MODULE}' is not available, the simulations "
              "are run but not evaluated")
    else:
        try:
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
        except ValueError as e:
            sys.exit(f"error: {e}")
    print(f"Input file: {input_file}")
    print(f"Build directory: {build_dir}")
    print(f"Borehole: L = {setup['length']} m, r_inner = {setup['r_inner']} m, "
          f"r_outer = {setup['r_outer']} m, {setup['num_segments']} elements at refinement 0")
    print(f"Rock grid at refinement 0: Cells = {' | '.join(BASE_CELLS)}, "
          f"Grading = {' | '.join(GRADING)} (surface), "
          f"Cells = {' '.join(str(count) for count in KERNEL_BASE_CELLS)} (kernel)")
    if MISSING_MODULE is None:
        print(f"Ramey: Re = {ramey.Re:.1f}, Pr = {ramey.Pr:.3f}, Nu = {ramey.Nu:.4f}, "
              f"h = {ramey.h:.4f} W/m²K, U = {ramey.U:.4f} W/m²K")

    # every run gets its own directory below this one, next to the figures
    output_dir = os.path.join(build_dir, OUTPUT_DIR)
    os.makedirs(output_dir, exist_ok=True)

    refinements = sorted(set(args.refinements))
    variants = [args.variant] if args.variant else list(VARIANTS)

    executables = {variant: os.path.join(build_dir, VARIANTS[variant]) for variant in variants}
    if not args.no_run:
        if not args.no_build:
            build_executables(build_dir, [VARIANTS[variant] for variant in variants])
        for path in executables.values():
            if not os.access(path, os.X_OK):
                sys.exit(f"error: '{path}' is not an executable, "
                         "use --build-dir to point to the build directory")

    run_dirs = {(variant, refinement):
                    os.path.join(output_dir, f"{variant}_refinement_{refinement}")
                for variant in variants for refinement in refinements}
    if not args.no_run:
        for (variant, refinement), run_dir in run_dirs.items():
            run_simulation(variant, executables[variant], input_file, setup, refinement,
                           run_dir, extra_args, f"{variant} refinement {refinement}")

    if MISSING_MODULE is not None:
        print(f"\nnote: the Python module '{MISSING_MODULE}' is not available, the results "
              f"in {output_dir} are not evaluated")
        return

    results = {}
    for variant in variants:
        for refinement in refinements:
            run_dir = run_dirs[variant, refinement]
            run_results = load_results(run_dir)
            if not run_results:
                sys.exit(f"error: no results found for {variant} refinement {refinement} in {run_dir}")
            results[variant, refinement] = run_results

    if not any("profile_data" in r for r in results.values()):
        print("\nnote: python3-vtk not available or no vtp output — "
              "the temperature profile along the borehole is skipped")

    # the solutions are shown on the finest grid, the coarser ones enter the convergence plot
    finest = refinements[-1]
    plot_outlet(variants, finest, results, setup, ramey, not args.hide_ogs,
                os.path.join(output_dir, "outlet_temperature.png"))
    if any("profile_data" in r for r in results.values()):
        plot_profile(variants, finest, results, setup, ramey, not args.hide_ogs,
                     os.path.join(output_dir, "temperature_profile.png"))
    if len(refinements) > 1:
        plot_convergence(variants, refinements, results, setup,
                         os.path.join(output_dir, "refinement_convergence.png"))

    plt.show()


if __name__ == "__main__":
    main()
