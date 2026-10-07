#!/usr/bin/env python3
"""
Grid convergence test for the wellbore heat transport benchmark.

Options:
--variants: which discretizations of the 3D domain to compare (default: all)
--refinements: which refinements to run (default: 0 1)
--exe: executable of one variant as VARIANT=PATH (default: its name in the current directory)
--reuse: skip runs whose results already exist
--no-run: only evaluate and plot existing results

Arguments that are not recognized are forwarded to the simulation, e.g.

    python3 convergence_test.py --refinements 0 1 -TimeLoop.TEnd 86400
"""
import argparse
import configparser
import glob
import os
import re
import shutil
import subprocess
import sys
import time

import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
import numpy as np

from ogs_reference import OGS_DEPTH, OGS_LENGTH, OGS_OUTLET, OGS_PROFILE, OGS_T_END, OGS_TIME
from ramey import RameySolution

# ── Soil grid at refinement 0 ─────────────────────────────────────────────────
BASE_CELLS = ("10 10", "10 10", "1")
GRADING = ("-1.3 1.3", "-1.3 1.3", "1")

# ── Simulation variants (different discretizations of the 3D domain) ──────────
VARIANTS = {
    "tpfa": "test_wellbore_heat_transport_tpfa",
    "box": "test_wellbore_heat_transport_box",
}

# input file of the benchmark and directory for the runs and the figure, both
# relative to the current directory
INPUT_FILE = "params.input"
OUTPUT_DIR = "convergence"

# ── Plot style ────────────────────────────────────────────────────────────────
# One hue per variant, the refinements are shades of it, see shades()
VARIANT_COLORS = {"tpfa": "#0072B2", "box": "#D55E00"}
REFERENCE_COLOR = "#333333"
GUIDE_COLOR = "#9a9a9a"
OGS_COLOR = "#6A3D9A"


def shades(color, n):
    """n shades of color, from light (coarsest grid) to dark (finest grid)."""
    rgb = np.array(mcolors.to_rgb(color))
    if n == 1:
        return [tuple(rgb)]
    light = rgb + (1.0 - rgb) * 0.55  # mixed with white
    dark = rgb * 0.6                  # mixed with black
    return [tuple(light + (dark - light) * i / (n - 1)) for i in range(n)]


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

    coords = np.array(vertices)
    length = float(np.linalg.norm(coords.max(axis=0) - coords.min(axis=0)))
    first = elements[0]
    r_inner, r_outer = float(first[2]), float(first[3])
    return length, r_inner, r_outer, len(elements)


def collect_setup(input_file):
    """Everything the script needs to know about the simulation setup."""
    params = read_input(input_file)
    input_dir = os.path.dirname(os.path.abspath(input_file))

    grid_file = get(params, "Voids.Grid", "File")
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
        "injection_rate": get(params, "Voids.Problem", "InjectionRate", float),
        "injection_temperature": get(params, "Voids.Problem", "InjectionTemperature", float) - 273.15,
        "lambda_s": get(params, "Component", "SolidThermalConductivity", float),
        "rho_s": get(params, "Component", "SolidDensity", float),
        "c_p_s": get(params, "Component", "SolidHeatCapacity", float),
    }
    return setup


# ── Running the simulation ────────────────────────────────────────────────────
def run_simulation(executable, input_file, setup, refinement, run_dir, extra_args, reuse, tag):
    """Run one variant at one refinement in its own directory."""
    if reuse and glob.glob(os.path.join(run_dir, "*.csv")):
        print(f"[{tag}] reusing existing results in {run_dir}")
        return

    if os.path.exists(run_dir):
        shutil.rmtree(run_dir)
    os.makedirs(run_dir)

    command = [executable, os.path.abspath(input_file)]
    # the soil grid at refinement 0 comes from this script, not from params.input
    for direction, (cells, grading) in enumerate(zip(BASE_CELLS, GRADING)):
        command += [f"-Soil.Grid.Cells{direction}", cells]
        command += [f"-Soil.Grid.Grading{direction}", grading]
    command += ["-Soil.Grid.Refinement", str(refinement)]
    command += ["-Voids.Grid.Refinement", str(refinement)]
    # the grid file is given relative to the input file, make it absolute
    command += ["-Voids.Grid.File", setup["grid_file"]]
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
    times, temperatures = raw[:, 0] * 86400.0, raw[:, 1]  # time in s, temperature in °C
    # t = 0 would compare the initial condition with the Ramey solution, which
    # is not defined there either (the time function vanishes), so drop it
    mask = times > 0.0
    return times[mask], temperatures[mask]


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
def plot_solutions(variants, refinements, results, setup, ramey, path):
    """Outlet temperature over time and temperature profile at the final time."""
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 5.5))
    ax1e, ax2e = ax1.twinx(), ax2.twinx()
    # the OGS data belongs to the benchmark setup only
    show_ogs = (abs(setup["length"] - OGS_LENGTH) < 1e-6
                and abs(setup["t_end"] - OGS_T_END) < 1.0)
    # the shades of every variant, one per refinement
    colors = {variant: shades(VARIANT_COLORS[variant], len(refinements)) for variant in variants}

    times = np.linspace(0.0, setup["t_end"], 200)[1:]
    ax1.plot(times / 86400, ramey.outlet(times), color=REFERENCE_COLOR, linewidth=2,
             linestyle="--", label="Ramey (1962)")
    if show_ogs:
        mask = OGS_TIME > 0.0  # t = 0 is the initial condition, see load_outlet_temperature
        ax1.plot(OGS_TIME[mask] / 86400, OGS_OUTLET[mask], color=OGS_COLOR, linewidth=2,
                 label="OGS")
        ax1e.plot(OGS_TIME[mask] / 86400, OGS_OUTLET[mask] - ramey.outlet(OGS_TIME[mask]),
                  color=OGS_COLOR, linewidth=1.5, linestyle="--")
    for variant in variants:
        for i, refinement in enumerate(refinements):
            if "outlet_data" not in results[variant, refinement]:
                continue
            t, T = results[variant, refinement]["outlet_data"]
            color = colors[variant][i]
            ax1.plot(t / 86400, T, color=color, linewidth=2,
                     label=f"{variant} refinement {refinement}")
            ax1e.plot(t / 86400, T - ramey.outlet(t), color=color, linewidth=1.5,
                      linestyle="--")
    ax1.set_xlabel("Time (d)")
    ax1.set_ylabel("Outlet temperature (°C)")
    ax1.set_title("Outlet temperature over time")

    depths = np.linspace(0.0, setup["length"], 200)
    ax2.plot(depths, ramey(depths, setup["t_end"]), color=REFERENCE_COLOR, linewidth=2,
             linestyle="--", label="Ramey (1962)")
    if show_ogs:
        ax2.plot(OGS_DEPTH, OGS_PROFILE, color=OGS_COLOR, linewidth=2, label="OGS")
        ax2e.plot(OGS_DEPTH, OGS_PROFILE - ramey(OGS_DEPTH, setup["t_end"]),
                  color=OGS_COLOR, linewidth=1.5, linestyle="--")
    for variant in variants:
        for i, refinement in enumerate(refinements):
            if "profile_data" not in results[variant, refinement]:
                continue
            x, T = results[variant, refinement]["profile_data"]
            color = colors[variant][i]
            ax2.plot(x, T, color=color, linewidth=2, label=f"{variant} refinement {refinement}")
            ax2e.plot(x, T - ramey(x, setup["t_end"]), color=color, linewidth=1.5,
                      linestyle="--")
    ax2.set_xlabel("Distance along borehole (m)")
    ax2.set_ylabel("Fluid temperature (°C)")
    ax2.set_xlim(0, setup["length"])
    ax2.set_title(f"Temperature profile at t = {setup['t_end'] / 86400:.3g} d")

    for ax, error_ax in ((ax1, ax1e), (ax2, ax2e)):
        ax.grid(color=GUIDE_COLOR, linewidth=0.5, alpha=0.4)
        ax.set_axisbelow(True)
        ax.spines["top"].set_visible(False)
        error_ax.set_ylabel("Absolute error vs Ramey (°C, dashed)")
        error_ax.spines["top"].set_visible(False)
        # the legend holds the solution curves only; it is put on the twin
        # because that one is drawn on top
        error_ax.legend(*ax.get_legend_handles_labels(), frameon=False, fontsize=9)

    fig.tight_layout()
    fig.savefig(path, dpi=150)
    print(f"Saved: {path}")
    return fig


# ── Main ──────────────────────────────────────────────────────────────────────
def main():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--refinements", type=int, nargs="+", default=[0, 1],
                        help="refinements to run (default: 0 1)")
    parser.add_argument("--variants", nargs="+", choices=list(VARIANTS), default=list(VARIANTS),
                        help="discretizations of the 3D domain to compare (default: all)")
    parser.add_argument("--exe", action="append", default=[], metavar="VARIANT=PATH",
                        help="path to the executable of one variant "
                             "(default: its name in the current directory)")
    parser.add_argument("--reuse", action="store_true",
                        help="skip runs whose results already exist")
    parser.add_argument("--no-run", action="store_true",
                        help="only evaluate and plot existing results")
    args, extra_args = parser.parse_known_args()

    input_file = os.path.abspath(INPUT_FILE)
    if not os.path.exists(input_file):
        sys.exit(f"error: no {INPUT_FILE} in the current directory — "
                 "run the script from the build directory")

    setup = collect_setup(input_file)
    try:
        ramey = RameySolution(
            injection_temperature=setup["injection_temperature"],
            injection_rate=setup["injection_rate"],
            length=setup["length"],
            r_inner=setup["r_inner"],
            r_outer=setup["r_outer"],
            lambda_s=setup["lambda_s"],
            rho_s=setup["rho_s"],
            c_p_s=setup["c_p_s"],
        )
    except ValueError as e:
        sys.exit(f"error: {e}")
    print(f"Input file: {input_file}")
    print(f"Borehole: L = {setup['length']} m, r_inner = {setup['r_inner']} m, "
          f"r_outer = {setup['r_outer']} m, {setup['num_segments']} elements at refinement 0")
    print(f"Soil grid at refinement 0: Cells = {' | '.join(BASE_CELLS)}, "
          f"Grading = {' | '.join(GRADING)}")
    print(f"Ramey: Re = {ramey.Re:.1f}, Pr = {ramey.Pr:.3f}, Nu = {ramey.Nu:.4f}, "
          f"h = {ramey.h:.4f} W/m²K, U = {ramey.U:.4f} W/m²K")

    # every run gets its own directory below this one, next to the figure
    output_dir = os.path.abspath(OUTPUT_DIR)
    os.makedirs(output_dir, exist_ok=True)

    refinements = sorted(set(args.refinements))
    variants = [v for v in VARIANTS if v in args.variants]  # keep the order of VARIANTS

    executables = {variant: VARIANTS[variant] for variant in variants}
    for override in args.exe:
        variant, _, path = override.partition("=")
        if variant not in executables:
            sys.exit(f"error: --exe {override}: expected one of "
                     f"{', '.join(variants)} before the '='")
        executables[variant] = path
    executables = {variant: os.path.abspath(path) for variant, path in executables.items()}
    if not args.no_run:
        for variant, path in executables.items():
            if not os.access(path, os.X_OK):
                sys.exit(f"error: '{path}' is not an executable, "
                         f"use --exe {variant}=<path> to point to it")

    results = {}
    for variant in variants:
        for refinement in refinements:
            run_dir = os.path.join(output_dir, f"{variant}_refinement_{refinement}")
            if not args.no_run:
                run_simulation(executables[variant], input_file, setup, refinement,
                               run_dir, extra_args, args.reuse, f"{variant} refinement {refinement}")
            run_results = load_results(run_dir)
            if not run_results:
                sys.exit(f"error: no results found for {variant} refinement {refinement} in {run_dir}")
            results[variant, refinement] = run_results

    if not any("profile_data" in r for r in results.values()):
        print("\nnote: python3-vtk not available or no vtp output — "
              "the temperature profile along the borehole is skipped")

    plot_solutions(variants, refinements, results, setup, ramey,
                   os.path.join(output_dir, "convergence_temperatures.png"))

    plt.show()


if __name__ == "__main__":
    main()
