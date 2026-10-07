#!/usr/bin/env python3
"""
Grid convergence test for the wellbore heat transport benchmark, including the
kernel (distributed source) coupling method.

The surface coupling variants (tpfa, box) run on the graded tensor-product grid of
params.input, which resolves the near field around the wellbore. The kernel variant
runs on the equidistant grid of params_kernel.input, because the kernel coupling
manager locates bulk elements by index arithmetic. Its kernel width stays fixed
while the grid is refined, so the kernel runs converge to the solution of the kernel
model with that width. That solution depends on the width (the flux scaling factor
assumes a steady radial temperature profile within the kernel), so several widths
can be compared with --kernel-width-factors.

The refinement is split into a radial part (the horizontal directions of the soil grid)
and an axial part (the vertical direction of the soil grid and the wellbore grid), so
both can be studied separately. A refinement level is a pair (radial, axial). A
refinement halves the cells of its directions, on the graded grid too, so the grids
of all levels are nested and the radial grid does not depend on the axial refinement.

Options:
--variants: which coupling variants to compare (default: all)
--refinements: refinements in both directions at once (default: 0 1)
--radial: radial refinements, overrides --refinements for the radial direction
--axial: axial refinements, overrides --refinements for the axial direction
    Lists of equal length are paired element by element, e.g. --radial 0 1
    --axial 1 2 runs (0, 1) and (1, 2). Lists of different length are combined
    to every pair, e.g. --radial 2 --axial 0 1 2 runs (2, 0), (2, 1) and (2, 2),
    and --radial 0 1 2 --axial 1 2 runs all six pairs.
--kernel-width-factors: kernel widths of the kernel variant, as multiples of the outer
    wellbore radius (default: MixedDimension.KernelWidthFactor of params_kernel.input).
    Every width is a series of runs of its own.
--exe: executable of one variant as VARIANT=PATH (default: its name in the current directory)
--reuse: skip runs whose results already exist
--no-run: only evaluate and plot existing results
--all: plot every run found in the output directory (implies --no-run; the
    refinement and kernel width options are ignored, --variants still filters)

Arguments that are not recognized are forwarded to the simulation, e.g.

    python3 convergence_test_kernel.py --refinements 0 1 -TimeLoop.TEnd 86400
    python3 convergence_test_kernel.py --radial 0 1 2 --axial 1
    python3 convergence_test_kernel.py --variants kernel --kernel-width-factors 2.5 5 --radial 2 3
    python3 convergence_test_kernel.py --all
"""
import argparse
import configparser
import glob
import itertools
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
# The graded tensor-product grid of the surface coupling methods, per direction
# (x, y, z), between the Soil.Grid.Positions of the input file. A refinement halves
# every cell of this grid (see soil_grid_args).
BASE_CELLS = ("10 10", "10 10", "2")
GRADING = ("-1.3 1.3", "-1.3 1.3", "1")

# The equidistant grid of the kernel method, as cell counts per direction. The kernel
# coupling manager rejects refined bulk grids, so the script doubles the cell counts
# per refinement itself instead of passing Soil.Grid.Refinement.
KERNEL_BASE_CELLS = (20, 20, 10)

# the soil grid directions refined by the radial and the axial refinement
RADIAL, AXIAL = "radial", "axial"
DIRECTIONS = (RADIAL, RADIAL, AXIAL)

# ── Simulation variants (different couplings of the 1D and 3D domain) ─────────
VARIANTS = {
    "tpfa": "test_wellbore_heat_transport_tpfa",
    "box": "test_wellbore_heat_transport_box",
    "kernel": "test_wellbore_heat_transport_kernel",
}

# input file per variant, relative to the current directory
INPUT_FILES = {
    "tpfa": "params.input",
    "box": "params.input",
    "kernel": "params_kernel.input",
}

# directory for the runs and the figure, relative to the current directory
OUTPUT_DIR = "convergence_kernel"

# the variant whose runs are parametrized by the kernel width
KERNEL = "kernel"

# ── Plot style ────────────────────────────────────────────────────────────────
# One hue per series and axial refinement, the radial refinements are shades of it,
# see shades() and level_colors()
VARIANT_COLORS = {"tpfa": "#0072B2", "box": "#D55E00", "kernel": "#009E73"}
# hues of further kernel widths and axial refinements, the first of a variant gets
# the variant hue
EXTRA_COLORS = ("#CC79A7", "#E69F00", "#56B4E9", "#882255", "#999933", "#44AA99",
                "#AA4499", "#661100")
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


# everything that has to be equal for the variants to be comparable with each
# other and with the Ramey solution; the soil grid is deliberately not part of it
COMPARABLE_KEYS = ("length", "r_inner", "r_outer", "t_end", "injection_rate",
                   "injection_temperature", "lambda_s", "rho_s", "c_p_s")


def common_setup(setups):
    """The setup shared by all variants, or an error if their input files disagree."""
    reference_variant, reference = next(iter(setups.items()))
    for variant, setup in setups.items():
        differing = [key for key in COMPARABLE_KEYS if setup[key] != reference[key]]
        if differing:
            sys.exit(f"error: {INPUT_FILES[variant]} and {INPUT_FILES[reference_variant]} "
                     f"disagree on {', '.join(differing)}, so the variants would not "
                     "describe the same problem")
    return reference


# ── Refinement levels ─────────────────────────────────────────────────────────
def refinement_levels(radial, axial):
    """Combine the radial and axial refinements to (radial, axial) levels.

    Lists of equal length are paired element by element, otherwise every
    radial refinement is combined with every axial refinement.
    """
    if min(radial + axial) < 0:
        sys.exit("error: refinements have to be non-negative")
    if len(radial) == len(axial):
        levels = zip(radial, axial)
    else:
        levels = itertools.product(radial, axial)
    # drop duplicates, keep the order from coarse to fine
    return sorted(set(levels), key=lambda level: (sum(level), level))


def level_name(level):
    """Label of a refinement level in the output and the plot."""
    radial, axial = level
    return f"radial {radial} axial {axial}"


# ── Series of runs ────────────────────────────────────────────────────────────
# A series is a variant together with its kernel width factor (None for the surface
# coupling variants), refined over all refinement levels.
def series_name(series):
    """Label of a series in the output and the plot."""
    variant, width = series
    return variant if width is None else f"{variant} width {width:g}"


def level_colors(series_levels):
    """The color of every (series, level).

    The runs of one series with the same axial refinement form a group of one hue,
    their radial refinements are shades of it. The first group of every variant
    gets the variant hue, all further groups get one of EXTRA_COLORS.
    """
    colors, used_variants, extra_count = {}, set(), 0
    for series, levels in series_levels.items():
        variant, _ = series
        for axial in sorted({axial for _, axial in levels}):
            group = [level for level in levels if level[1] == axial]
            if variant not in used_variants:
                color = VARIANT_COLORS[variant]
                used_variants.add(variant)
            else:
                color = EXTRA_COLORS[extra_count % len(EXTRA_COLORS)]
                extra_count += 1
            # the levels are ordered from coarse to fine, so are the shades
            for level, shade in zip(group, shades(color, len(group))):
                colors[series, level] = shade
    return colors


def level_dir(series, level):
    """Directory name of the run of one series at one refinement level."""
    variant, width = series
    radial, axial = level
    prefix = variant if width is None else f"{variant}_width_{width:g}"
    return f"{prefix}_radial_{radial}_axial_{axial}"


# the inverse of level_dir
RUN_DIR_PATTERN = re.compile(r"^(?P<variant>[a-z]+)(?:_width_(?P<width>[^_]+))?"
                             r"_radial_(?P<radial>\d+)_axial_(?P<axial>\d+)$")


def discover_runs(output_dir, variants):
    """All runs in output_dir as {series: [levels]}, for the given variants only."""
    runs = {}
    if not os.path.isdir(output_dir):
        return runs
    for name in os.listdir(output_dir):
        match = RUN_DIR_PATTERN.match(name)
        if (not match or match["variant"] not in variants
                or not os.path.isdir(os.path.join(output_dir, name))):
            continue
        variant, width = match["variant"], match["width"]
        # only the kernel variant is parametrized by the width
        if (width is None) == (variant == KERNEL):
            continue
        try:
            width = None if width is None else float(width)
        except ValueError:
            continue
        level = (int(match["radial"]), int(match["axial"]))
        runs.setdefault((variant, width), []).append(level)
    # order the series like VARIANTS and from narrow to wide kernels
    order = list(VARIANTS)
    return {series: sorted(runs[series], key=lambda level: (sum(level), level))
            for series in sorted(runs, key=lambda s: (order.index(s[0]), s[1] or 0.0))}


# ── Running the simulation ────────────────────────────────────────────────────
def graded_vertices(lower, upper, cells, grading):
    """Vertices of a graded interval without the upper one, placed like DuMux places them.

    The cell sizes are a geometric sequence, growing from lower to upper for a positive
    grading factor and shrinking for a negative one. A factor and its inverse are equivalent.
    """
    ratio = max(abs(grading), 1.0 / abs(grading))
    sizes = (ratio if grading > 0.0 else 1.0 / ratio) ** np.arange(cells)
    return lower + (upper - lower) * np.append(0.0, np.cumsum(sizes[:-1])) / sizes.sum()


def soil_grid_args(variant, level, params):
    """Command line arguments for the soil grid of one variant at one refinement level."""
    refinement = dict(zip((RADIAL, AXIAL), level))
    if variant == "kernel":
        # the kernel coupling manager does not support refined bulk grids
        cells = (count * 2**refinement[direction]
                 for count, direction in zip(KERNEL_BASE_CELLS, DIRECTIONS))
        return ["-Soil.Grid.Cells", " ".join(str(count) for count in cells)]

    # A refinement halves the cells like Soil.Grid.Refinement does, but DuMux refines all
    # directions at once and a halved graded grid is no geometric sequence anymore. So
    # every cell of the grid at refinement 0 becomes an interval of the tensor-product
    # grid of its own, with 2^refinement equidistant cells.
    args = []
    for direction, (cells, grading, name) in enumerate(zip(BASE_CELLS, GRADING, DIRECTIONS)):
        positions = get(params, "Soil.Grid", f"Positions{direction}").split()
        if len(positions) != len(cells.split()) + 1:
            sys.exit(f"error: [Soil.Grid] Positions{direction} of the input file does not "
                     "match the intervals of BASE_CELLS")
        vertices = []
        for lower, upper, count, factor in zip(positions, positions[1:],
                                               cells.split(), grading.split()):
            vertices += list(graded_vertices(float(lower), float(upper), int(count), float(factor)))
        vertices.append(float(positions[-1]))
        intervals = len(vertices) - 1
        args += [f"-Soil.Grid.Positions{direction}", " ".join(f"{v:.17g}" for v in vertices)]
        args += [f"-Soil.Grid.Cells{direction}", " ".join([str(2**refinement[name])] * intervals)]
        args += [f"-Soil.Grid.Grading{direction}", " ".join(["1"] * intervals)]
    return args + ["-Soil.Grid.Refinement", "0"]


def run_simulation(executable, input_file, setup, series, level, run_dir,
                   extra_args, reuse, tag):
    """Run one series at one refinement level in its own directory."""
    variant, width = series
    if reuse and glob.glob(os.path.join(run_dir, "*.csv")):
        print(f"[{tag}] reusing existing results in {run_dir}")
        return

    if os.path.exists(run_dir):
        shutil.rmtree(run_dir)
    os.makedirs(run_dir)

    command = [executable, os.path.abspath(input_file)]
    # the soil grid at refinement 0 comes from this script, not from the input file
    command += soil_grid_args(variant, level, setup["params"])
    # the wellbore runs along the axis, so it follows the axial refinement
    command += ["-Voids.Grid.Refinement", str(level[1])]
    # the grid file is given relative to the input file, make it absolute
    command += ["-Voids.Grid.File", setup["grid_file"]]
    if width is not None:
        command += ["-MixedDimension.KernelWidthFactor", f"{width:g}"]
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
def plot_solutions(all_series, levels, results, setup, ramey, path):
    """Outlet temperature over time and temperature profile at the final time."""
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 5.5))
    ax1e, ax2e = ax1.twinx(), ax2.twinx()
    # the OGS data belongs to the benchmark setup only
    show_ogs = (abs(setup["length"] - OGS_LENGTH) < 1e-6
                and abs(setup["t_end"] - OGS_T_END) < 1.0)
    # the levels with results of every series, and their colors
    series_levels = {series: [level for level in levels if (series, level) in results]
                     for series in all_series}
    colors = level_colors(series_levels)

    times = np.linspace(0.0, setup["t_end"], 200)[1:]
    ax1.plot(times / 86400, ramey.outlet(times), color=REFERENCE_COLOR, linewidth=2,
             linestyle="--", label="Ramey (1962)")
    if show_ogs:
        mask = OGS_TIME > 0.0  # t = 0 is the initial condition, see load_outlet_temperature
        ax1.plot(OGS_TIME[mask] / 86400, OGS_OUTLET[mask], color=OGS_COLOR, linewidth=2,
                 label="OGS")
        ax1e.plot(OGS_TIME[mask] / 86400, OGS_OUTLET[mask] - ramey.outlet(OGS_TIME[mask]),
                  color=OGS_COLOR, linewidth=1.5, linestyle="--")
    for series in all_series:
        for level in series_levels[series]:
            if "outlet_data" not in results[series, level]:
                continue
            t, T = results[series, level]["outlet_data"]
            color = colors[series, level]
            ax1.plot(t / 86400, T, color=color, linewidth=2,
                     label=f"{series_name(series)} {level_name(level)}")
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
    for series in all_series:
        for level in series_levels[series]:
            if "profile_data" not in results[series, level]:
                continue
            x, T = results[series, level]["profile_data"]
            color = colors[series, level]
            ax2.plot(x, T, color=color, linewidth=2,
                     label=f"{series_name(series)} {level_name(level)}")
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
                        help="refinements in both directions (default: 0 1)")
    parser.add_argument("--radial", type=int, nargs="+",
                        help="radial refinements of the soil grid (default: --refinements)")
    parser.add_argument("--axial", type=int, nargs="+",
                        help="axial refinements of the soil and the wellbore grid "
                             "(default: --refinements)")
    parser.add_argument("--variants", nargs="+", choices=list(VARIANTS), default=list(VARIANTS),
                        help="coupling variants to compare (default: all)")
    parser.add_argument("--kernel-width-factors", type=float, nargs="+", metavar="FACTOR",
                        help="kernel widths of the kernel variant as multiples of the outer "
                             "wellbore radius (default: MixedDimension.KernelWidthFactor "
                             "of params_kernel.input)")
    parser.add_argument("--exe", action="append", default=[], metavar="VARIANT=PATH",
                        help="path to the executable of one variant "
                             "(default: its name in the current directory)")
    parser.add_argument("--reuse", action="store_true",
                        help="skip runs whose results already exist")
    parser.add_argument("--no-run", action="store_true",
                        help="only evaluate and plot existing results")
    parser.add_argument("--all", action="store_true",
                        help="plot every run found in the output directory (implies --no-run)")
    args, extra_args = parser.parse_known_args()

    levels = refinement_levels(args.radial or args.refinements,
                               args.axial or args.refinements)
    variants = [v for v in VARIANTS if v in args.variants]  # keep the order of VARIANTS

    # every run gets its own directory below this one, next to the figure
    output_dir = os.path.abspath(OUTPUT_DIR)

    if args.all:
        args.no_run = True
        runs = discover_runs(output_dir, variants)
        if not runs:
            sys.exit(f"error: no runs of {', '.join(variants)} found in {output_dir}")
        variants = [v for v in variants if any(series[0] == v for series in runs)]
        levels = sorted({level for series_levels in runs.values() for level in series_levels},
                        key=lambda level: (sum(level), level))

    input_files = {}
    for variant in variants:
        path = os.path.abspath(INPUT_FILES[variant])
        if not os.path.exists(path):
            sys.exit(f"error: no {INPUT_FILES[variant]} in the current directory — "
                     "run the script from the build directory")
        input_files[variant] = path

    setups = {variant: collect_setup(path) for variant, path in input_files.items()}
    setup = common_setup(setups)

    # the kernel width is set by the script, one series per width
    if any("KernelWidthFactor" in arg for arg in extra_args):
        sys.exit("error: set the kernel width with --kernel-width-factors")
    all_series = list(runs) if args.all else []
    for variant in ([] if args.all else variants):
        if variant != KERNEL:
            all_series.append((variant, None))
            continue
        widths = args.kernel_width_factors or [
            get(setups[variant]["params"], "MixedDimension", "KernelWidthFactor", float)]
        if min(widths) <= 0.0:
            sys.exit("error: kernel width factors have to be positive")
        # drop duplicates, keep the order from narrow to wide
        all_series += [(variant, width) for width in sorted(set(widths))]
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
    for variant in variants:
        print(f"Input file ({variant}): {input_files[variant]}")
    print(f"Borehole: L = {setup['length']} m, r_inner = {setup['r_inner']} m, "
          f"r_outer = {setup['r_outer']} m, {setup['num_segments']} elements at refinement 0")
    print(f"Soil grid at refinement 0: graded Cells = {' | '.join(BASE_CELLS)}, "
          f"Grading = {' | '.join(GRADING)}, cells halved per refinement")
    if "kernel" in variants:
        print("Soil grid at refinement 0 (kernel): equidistant Cells = "
              f"{' '.join(str(count) for count in KERNEL_BASE_CELLS)}, doubled per refinement")
    for variant, width in all_series:
        if width is not None:
            print(f"Kernel width ({variant}): factor {width:g} = {width * setup['r_outer']:.4g} m")
    print("Refinement levels: " + ", ".join(f"({level_name(level)})" for level in levels))
    print(f"Ramey: Re = {ramey.Re:.1f}, Pr = {ramey.Pr:.3f}, Nu = {ramey.Nu:.4f}, "
          f"h = {ramey.h:.4f} W/m²K, U = {ramey.U:.4f} W/m²K")

    os.makedirs(output_dir, exist_ok=True)

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
    for series in all_series:
        variant = series[0]
        for level in (runs[series] if args.all else levels):
            run_dir = os.path.join(output_dir, level_dir(series, level))
            tag = f"{series_name(series)} {level_name(level)}"
            if not args.no_run:
                run_simulation(executables[variant], input_files[variant], setups[variant],
                               series, level, run_dir, extra_args, args.reuse, tag)
            run_results = load_results(run_dir)
            if not run_results:
                if args.all:
                    print(f"note: skipping {tag}, no results in {run_dir}")
                    continue
                sys.exit(f"error: no results found for {tag} in {run_dir}")
            results[series, level] = run_results

    if not any("profile_data" in r for r in results.values()):
        print("\nnote: python3-vtk not available or no vtp output — "
              "the temperature profile along the borehole is skipped")

    plot_solutions(all_series, levels, results, setup, ramey,
                   os.path.join(output_dir, "convergence_temperatures.png"))

    plt.show()


if __name__ == "__main__":
    main()
