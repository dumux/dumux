#!/usr/bin/env python3
# SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later
"""
Run the 1pni heat conduction benchmark and plot the results against the analytical solution.

Usage (run from the build directory of the test, which contains params.input):
  python3 plot_results.py <executable> [<executable> ...] [--output-dir <directory>]

The executables are test_1pni_conduction_box, test_1pni_conduction_tpfa and/or test_1pni_conduction_mpfa.
The script writes to the output directory (default: current directory):
  - 1pni_1d_conduction_benchmark_lineplot.png: numerical and exact temperature profiles
    at five output times (from the first cell-centered executable, else the first one)
  - 1pni_1d_conduction_convergence_box.png: space and time convergence of the box scheme
    (continuous L2 norm of the error)
  - 1pni_1d_conduction_convergence_cc.png: space and time convergence of the cell-centered
    schemes (discrete L2 norm of the error)
and prints the convergence tables. Box and cell-centered errors use different norms
and are therefore shown in separate figures.
"""

import argparse
import glob
import os
import subprocess
import xml.etree.ElementTree as ET
from pathlib import Path

from convergencetest import END_TIME, compute_rates, print_table, run_study

INPUT_FILE = "params.input"
BASE_NAME = "1pni_1d_conduction_benchmark"
NUM_PROFILES = 5
TEMPERATURE_FIELD = "T"
EXACT_FIELD = "temperatureExact"

# the thermal front stays well within this part of the 5 m domain
X_MAX = 1.5
TEMPERATURE_LIMITS = (289.8, 300.2)

# ordinal blue ramp (light -> dark) encoding the ordered output times
PROFILE_COLORS = ("#86b6ef", "#5598e7", "#2a78d6", "#1c5cab", "#104281")
# categorical colors in fixed order: first scheme blue, second scheme orange
SCHEME_COLORS = ("#2a78d6", "#eb6834")
SCHEME_MARKERS = ("o", "s")
REFERENCE_COLOR = "#7a7a75"


def require_plot_modules():
    try:
        import matplotlib.pyplot as plt
        import pyvista as pv
    except ImportError as error:
        raise SystemExit(
            "Post-processing requires PyVista and Matplotlib. "
            "Install them with: python3 -m pip install pyvista matplotlib"
        ) from error
    return pv, plt


def scheme_of(executable):
    return os.path.basename(executable).rsplit("_", 1)[-1]


# ---------------------------------------------------------------- line plot

def run_profiles(executable):
    """Run the benchmark setup with frequent output and return (time, vtu file) of all outputs."""
    for pattern in (f"{BASE_NAME}-*.vtu", f"{BASE_NAME}.pvd"):
        for file_name in glob.glob(pattern):
            os.remove(file_name)

    exe = executable if os.path.isabs(executable) else "./" + executable
    subprocess.run([exe, INPUT_FILE, "-Problem.Name", BASE_NAME, "-Problem.OutputInterval", "5"],
                   check=True, stdout=subprocess.DEVNULL)

    collection = ET.parse(f"{BASE_NAME}.pvd").getroot()
    results = [(float(d.attrib["timestep"]), d.attrib["file"]) for d in collection.iter("DataSet")]
    return sorted(results)


def representative_times(results, count):
    """Pick a logarithmic spread of output times.

    The initial output is skipped: the exact solution is only updated inside the
    time loop and is therefore still zero in the first output file.
    """
    transient = [result for result in results if result[0] > 0.0]
    if len(transient) <= count:
        return transient

    first, last = transient[0][0], transient[-1][0]
    targets = [first*(last/first)**(index/(count - 1)) for index in range(count)]
    picked = sorted({min(range(len(transient)), key=lambda i: abs(transient[i][0] - target))
                     for target in targets})
    return [transient[index] for index in picked]


def sorted_profile(mesh, field):
    """Return the x-coordinates and values of a cell or vertex field, sorted along x.

    Vertex data of the quasi-1D grid carries several values per x-position, which
    are averaged, so that box and cell-centered schemes are handled alike.
    """
    if field in mesh.cell_data:
        positions, values = mesh.cell_centers().points[:, 0], mesh.cell_data[field]
    elif field in mesh.point_data:
        positions, values = mesh.points[:, 0], mesh.point_data[field]
    else:
        raise KeyError(f"Missing field '{field}'")

    profile = {}
    for position, value in zip(positions, values):
        profile.setdefault(round(float(position), 9), []).append(float(value))

    x = sorted(profile)
    return x, [sum(profile[position])/len(profile[position]) for position in x]


def format_time(time):
    for limit, factor, unit in ((60.0, 1.0, "s"), (3600.0, 60.0, "min"), (86400.0, 3600.0, "h")):
        if time < limit:
            return f"{time/factor:.3g} {unit}"
    return f"{time/86400.0:.3g} d"


def read_params(input_file=INPUT_FILE):
    """Read the parameters of the input file into a dict with keys of the form Group.Key."""
    import configparser
    parser = configparser.ConfigParser(inline_comment_prefixes=("#",))
    parser.optionxform = str  # keep the case of the keys
    parser.read(input_file)
    return {f"{section}.{key}": value for section in parser.sections() for key, value in parser[section].items()}


def line_plot_title(executable):
    """Scheme, grid and maximum time step size of the run shown in the line plot."""
    params = read_params()
    cells = [int(n) for n in params["Grid.Cells"].split()]
    lower = [float(v) for v in params.get("Grid.LowerLeft", "0 0").split()]
    upper = [float(v) for v in params["Grid.UpperRight"].split()]
    cell_size = (upper[0] - lower[0])/cells[0]
    max_dt = float(params["TimeLoop.MaxTimeStepSize"])
    return (f"{scheme_of(executable)}, {' × '.join(map(str, cells))} cells ($h$ = {cell_size*1000:g} mm), "
            rf"maximum time step size $\Delta t_\mathrm{{max}}$ = {max_dt:g} s")


def create_line_plot(executable, image_file):
    pv, plt = require_plot_modules()
    results = representative_times(run_profiles(executable), NUM_PROFILES)

    fig, ax = plt.subplots(figsize=(8.2, 4.8), constrained_layout=True)
    for index, (time, vtu) in enumerate(results):
        mesh = pv.read(vtu)
        x, temperature = sorted_profile(mesh, TEMPERATURE_FIELD)
        exact_x, exact = sorted_profile(mesh, EXACT_FIELD)
        color = PROFILE_COLORS[index % len(PROFILE_COLORS)]
        # exact solution: thin black line below; numerical solution: filled crosses at the degrees of freedom
        ax.plot(exact_x, exact, linewidth=1.2, color="black", zorder=2)
        ax.plot(x, temperature, label=rf"$T$ numerical ($t$ = {format_time(time)})",
                linestyle="none", marker="X", markersize=8, color=color,
                markeredgewidth=0, zorder=3)

    ax.plot([], [], label=r"$T$ exact", linewidth=1.2, color="black")
    ax.set_xlabel("x [m]")
    ax.set_ylabel(r"Temperature $T$ [K]")
    ax.set_xlim(0, X_MAX)
    ax.set_ylim(*TEMPERATURE_LIMITS)
    ax.set_title(line_plot_title(executable))
    ax.grid(True, alpha=0.3)
    ax.legend()
    fig.savefig(image_file, dpi=200)
    plt.close(fig)


# ---------------------------------------------------------------- convergence

def run_convergence(executable):
    """Run the full space and time studies and print the tables."""
    data = {}
    for study in ("space", "time"):
        hs, dts, errors = run_study(executable, study, INPUT_FILE, full=True,
                                    name=f"{BASE_NAME}_{scheme_of(executable)}_{study}")
        sizes = hs if study == "space" else dts
        rates = compute_rates(sizes, errors)
        print(f"\n{study.capitalize()} convergence, {scheme_of(executable)}:")
        print_table(study, hs, dts, errors, rates)
        data[study] = (sizes, errors)
    return data


def add_reference_slope(ax, sizes, errors, order):
    """Dashed line with the expected order, anchored at the coarsest point and shifted below the data."""
    anchor_size, anchor_error = sizes[0], 0.5*errors[0]
    reference = [anchor_error*(size/anchor_size)**order for size in sizes]
    ax.plot(sizes, reference, color=REFERENCE_COLOR, linestyle="--", linewidth=1.2,
            label=f"order {order}")


def create_convergence_plot(results, title, norm_label, image_file):
    """One figure for one scheme family: space study (left) and time study (right)."""
    _, plt = require_plot_modules()
    fig, (ax_space, ax_time) = plt.subplots(1, 2, figsize=(9.6, 4.2), constrained_layout=True)

    for index, (scheme, data) in enumerate(results.items()):
        style = dict(color=SCHEME_COLORS[index], marker=SCHEME_MARKERS[index], markersize=8 - 2*index,
                     linewidth=2, label=scheme)
        ax_space.plot(*data["space"], **style)
        ax_time.plot(*data["time"], **style)

    first = next(iter(results.values()))
    add_reference_slope(ax_space, *first["space"], order=2)
    add_reference_slope(ax_time, *first["time"], order=1)

    ax_space.set_xlabel("cell size $h$ [m]")
    ax_space.set_title(r"space ($\Delta t \propto h^2$)")
    ax_time.set_xlabel(r"time step size $\Delta t$ [s]")
    ax_time.set_title("time (3200 cells)")
    for ax in (ax_space, ax_time):
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_ylabel(f"{norm_label} [K]")
        ax.grid(True, which="both", alpha=0.3)
        ax.legend()
    mantissa, exponent = f"{END_TIME:e}".split("e")
    time_label = rf"$10^{{{int(exponent)}}}$" if float(mantissa) == 1.0 else f"{END_TIME:g}"
    fig.suptitle(rf"{title} ($t$ = {time_label} s = {END_TIME/86400:.3g} d)")
    fig.savefig(image_file, dpi=200)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("executables", nargs="+")
    parser.add_argument("--output-dir", default=".", help="directory of the figures (default: current directory)")
    args = parser.parse_args()

    out_dir = Path(args.output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    cell_centered = [exe for exe in args.executables if scheme_of(exe) in ("tpfa", "mpfa")]
    box = [exe for exe in args.executables if scheme_of(exe) == "box"]

    print("Creating line plot...")
    create_line_plot((cell_centered or args.executables)[0], out_dir/f"{BASE_NAME}_lineplot.png")

    if box:
        results = {scheme_of(exe): run_convergence(exe) for exe in box}
        create_convergence_plot(results, "Box scheme", r"$L^2$ error", out_dir/"1pni_1d_conduction_convergence_box.png")
    if cell_centered:
        results = {scheme_of(exe): run_convergence(exe) for exe in cell_centered}
        create_convergence_plot(results, "Cell-centered schemes", r"discrete $L^2$ error",
                                out_dir/"1pni_1d_conduction_convergence_cc.png")


if __name__ == "__main__":
    main()
