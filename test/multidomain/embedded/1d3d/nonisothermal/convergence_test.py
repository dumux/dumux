#!/usr/bin/env python3
"""
Grid convergence test for the wellbore heat transport benchmark.

The simulation is run several times, each time on a finer grid: the level-0 soil
grid is defined by ``BASE_CELLS``/``GRADING`` below (the cells and grading given
in params.input are ignored), and both domains are refined with the DuMux
parameters ``Soil.Grid.Refinement`` and ``Voids.Grid.Refinement``, i.e. every
level halves the mesh size in each direction.

For every level the outlet temperature over time and the fluid temperature
profile along the borehole at the final time are plotted together with the
Ramey (1962) analytical solution (same solution as in analytical_solution.py),
the OGS reference of the benchmark, and the absolute error against Ramey on a
second y axis.

Usage (run it from the build directory, where the executable lives)::

    python3 convergence_test.py --levels 0 1 2

Arguments that are not recognized are forwarded to the simulation, which is
handy for cheap trial runs::

    python3 convergence_test.py --levels 0 1 -TimeLoop.TEnd 86400

Note that each refinement level multiplies the number of soil cells by 8, so
levels beyond 1-2 quickly become expensive.
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

import matplotlib.pyplot as plt
import numpy as np

# ── Fixed physical properties ─────────────────────────────────────────────────
# Everything that cannot be read from params.input (fluid properties are
# hard-coded in simpleh2o.hh, the undisturbed formation temperature in
# problem_soil.hh).  Values are identical to analytical_solution.py.
T_D = 55.0            # undisturbed formation temperature, °C
RHO_F = 1000.0        # fluid density, kg/m³
C_P_F = 4190.0        # fluid specific heat, J/kg/K
MU_F = 1.14e-3        # dynamic viscosity, Pa·s
LAMBDA_F = 0.59       # fluid thermal conductivity, W/m/K
LAMBDA_G = 0.73       # grout thermal conductivity, W/m/K
LAMBDA_PI = 1.3       # pipe wall thermal conductivity, W/m/K
T_PI = 0.00587        # pipe wall thickness, m

# ── Soil grid at refinement level 0 ───────────────────────────────────────────
# One entry per coordinate direction, holding the cells / grading factors of
# every interval of the corresponding Soil.Grid.Positions in params.input (the
# default grid is graded with two intervals in x and y and one in z).  The
# grading factors belong to the cell numbers, so they are set here as well.
# Higher levels are created by Soil.Grid.Refinement from this grid.
BASE_CELLS = ("10 10", "10 10", "10")
GRADING = ("-1.3 1.3", "-1.3 1.3", "1")

# ── Plot style (Okabe-Ito, colorblind safe; grey is reserved for guides) ──────
LEVEL_COLORS = ["#0072B2", "#D55E00", "#009E73", "#CC79A7", "#56B4E9", "#E69F00"]
REFERENCE_COLOR = "#333333"
GUIDE_COLOR = "#9a9a9a"
OGS_COLOR = "#6A3D9A"

# ── OGS reference solution ────────────────────────────────────────────────────
# Result of the OGS benchmark that this test reproduces (same data as in
# analytical_solution.py): outlet temperature over five days and the fluid
# temperature along the 30 m borehole at the final time.  It only describes this
# one setup, so it is plotted only if params.input asks for the same one.
OGS_LENGTH = 30.0        # borehole length of the benchmark, m
OGS_T_END = 5 * 86400.0  # final time of the benchmark, s

OGS_TIME = np.linspace(0.0, OGS_T_END, 31)  # s
OGS_OUTLET = np.array([  # outlet temperature, °C
    20.        , 25.67074983, 26.67948325, 26.54800871, 26.27755879, 26.03968084,
    25.84599915, 25.68828698, 25.55846907, 25.45030474, 25.35907264, 25.28119084,
    25.21392656, 25.15518214, 25.10333715, 25.05713156, 25.01557889, 24.97790140,
    24.94348150, 24.91182514, 24.88253415, 24.85528518, 24.82981365, 24.80590137,
    24.78336707, 24.76205898, 24.74184908, 24.72262861, 24.70430449, 24.68679654,
    24.67003519,
])

OGS_DEPTH = np.linspace(0.0, OGS_LENGTH, 101)  # m
OGS_PROFILE = np.array([  # fluid temperature at t = 5 d, °C
    20.        , 20.0499748 , 20.0999496 , 20.1499244 , 20.19959088, 20.24910318,
    20.29861549, 20.34812923, 20.39764582, 20.44716241, 20.49667901, 20.5457412 ,
    20.59480339, 20.64386559, 20.69293331, 20.74200381, 20.7910743 , 20.83999301,
    20.88860817, 20.93722332, 20.98583847, 21.03446673, 21.08309499, 21.13172325,
    21.18004725, 21.22821912, 21.276391  , 21.32456885, 21.37275866, 21.42094847,
    21.46913828, 21.5168706 , 21.56460292, 21.61233525, 21.66008277, 21.70783788,
    21.755593  , 21.80319523, 21.85049171, 21.89778819, 21.94508467, 21.99240882,
    22.03973297, 22.08705711, 22.1340747 , 22.18093901, 22.22780331, 22.27467848,
    22.32157535, 22.36847222, 22.4153691 , 22.46180488, 22.50824065, 22.55467643,
    22.6011372 , 22.64761047, 22.69408374, 22.74040287, 22.78641374, 22.8324246 ,
    22.87843546, 22.92448876, 22.97054207, 23.01659537, 23.06233949, 23.10792902,
    23.15351855, 23.19912389, 23.24476083, 23.29039778, 23.33603473, 23.38120648,
    23.42637824, 23.47154999, 23.51675669, 23.56198086, 23.60720503, 23.65227365,
    23.69703115, 23.74178866, 23.78654617, 23.83136112, 23.87617606, 23.92099101,
    23.96549383, 24.00984059, 24.05418734, 24.09855493, 24.14296419, 24.18737345,
    24.2317827 , 24.27572213, 24.31966155, 24.36360098, 24.40758592, 24.45159361,
    24.4956013 , 24.5394493 , 24.58297793, 24.62650656, 24.67003519,
])


# ── Ramey (1962) analytical solution ──────────────────────────────────────────
def outlet_temp(T_d, T_i, z, X):
    """Fluid temperature at distance z along the borehole."""
    return T_d + (T_i - T_d) * np.exp(-z / X)


def coefficient_x(q, rho_f, c_p_f, lambda_re, r_ref, U, f_t):
    """Characteristic length X of the exponential temperature decay."""
    return (q * rho_f * c_p_f) * (lambda_re + r_ref * U * f_t) / (2 * np.pi * r_ref * U * lambda_re)


def dimensionless_time(lambda_re, delta_t, rho_re, c_p_re, r_b):
    return lambda_re * delta_t / (rho_re * c_p_re * r_b**2)


def time_function(t_D):
    if t_D > 1.5:
        return (0.4063 + 0.5 * np.log(t_D)) * (1 + 0.6 / t_D)
    return 1.1281 * np.sqrt(t_D) * (1 - 0.3 * np.sqrt(t_D))


def nusselt_number(Re, Pr, r_pi, length):
    """Nusselt number of the pipe flow (laminar / transitional / turbulent)."""
    if Re < 2300:
        return 4.364
    entrance = 1 + (r_pi / length) ** (2 / 3)
    if Re > 10000:
        xi = (1.8 * np.log(Re) - 1.5) ** -2
        return (xi / 8 * Re * Pr) / (1 + 12.7 * np.sqrt(xi / 8) * (Pr ** (2 / 3) - 1)) * entrance
    gamma = (Re - 2300) / (10000 - 2300)
    xi_t = (1.8 * np.log(10000) - 1.5) ** -2
    Nu_t = (xi_t / 8 * 1e4 * Pr) / (1 + 12.7 * np.sqrt(xi_t / 8) * (Pr ** (2 / 3) - 1)) * entrance
    return (1 - gamma) * 4.364 + gamma * Nu_t


class RameySolution:
    """Ramey solution T(z, t) for the setup described by an input file."""

    def __init__(self, setup):
        self.T_i = setup["injection_temperature"]
        self.q = setup["injection_rate"]
        self.length = setup["length"]
        self.r_pi = setup["r_inner"]
        self.r_b = setup["r_outer"]
        self.lambda_re = setup["lambda_re"]
        self.rho_re = setup["rho_re"]
        self.c_p_re = setup["c_p_re"]

        v = self.q / (np.pi * self.r_pi**2)
        self.Pr = MU_F * C_P_F / LAMBDA_F
        self.Re = RHO_F * v * (2 * self.r_pi) / MU_F
        self.Nu = nusselt_number(self.Re, self.Pr, self.r_pi, self.length)
        self.h = LAMBDA_F * self.Nu / (2 * self.r_pi)

        # overall heat transfer coefficient, referred to the inner pipe radius
        r_po = self.r_pi + T_PI
        self.U = 1 / (1 / self.h + self.r_pi * (np.log(r_po / self.r_pi) / LAMBDA_PI
                                                + np.log(self.r_b / r_po) / LAMBDA_G))

    def __call__(self, z, t):
        """Temperature in °C at distance z [m] along the borehole and time t [s]."""
        z = np.asarray(z, dtype=float)
        f_t = time_function(dimensionless_time(self.lambda_re, t, self.rho_re, self.c_p_re, self.r_b))
        X = coefficient_x(self.q, RHO_F, C_P_F, self.lambda_re, self.r_pi, self.U, f_t)
        return outlet_temp(T_D, self.T_i, z, X)

    def outlet(self, times):
        """Outlet temperature in °C for an array of times [s]."""
        return np.array([float(self(self.length, t)) for t in times])


# ── Input file handling ───────────────────────────────────────────────────────
def read_input(path):
    """Read a DuMux input file into a nested dict {group: {key: value}}."""
    parser = configparser.ConfigParser(inline_comment_prefixes=("#",), strict=False)
    parser.optionxform = str  # keep the case of the keys
    with open(path) as f:
        parser.read_string("[__root__]\n" + f.read())
    return {section: dict(parser[section]) for section in parser.sections()}


def get(params, group, key, default=None, cast=str):
    try:
        return cast(params[group][key])
    except (KeyError, ValueError):
        return default


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

    grid_file = get(params, "Voids.Grid", "File", "./grids/pipe.dgf")
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
        "t_end": get(params, "TimeLoop", "TEnd", 432000.0, float),
        "injection_rate": get(params, "Voids.Problem", "InjectionRate", 2e-4, float),
        "injection_temperature": get(params, "Voids.Problem", "InjectionTemperature", 293.15, float) - 273.15,
        "lambda_re": get(params, "Component", "SolidThermalConductivity", 2.78018, float),
        "rho_re": get(params, "Component", "SolidDensity", 1800.0, float),
        "c_p_re": get(params, "Component", "SolidHeatCapacity", 1778.0, float),
    }
    return setup


def num_soil_cells(level):
    """Total number of cells of the 3D grid at the given refinement level."""
    cells = 8**level
    for entry in BASE_CELLS:
        cells *= sum(int(c) for c in entry.split())
    return cells


# ── Running the simulation ────────────────────────────────────────────────────
def find_executable(explicit, script_dir):
    if explicit:
        if not os.access(explicit, os.X_OK):
            sys.exit(f"error: '{explicit}' is not an executable")
        return os.path.abspath(explicit)

    candidates = []
    # the build directory of the dune module this script lives in
    module_root = script_dir
    while module_root != "/" and not os.path.exists(os.path.join(module_root, "dune.module")):
        module_root = os.path.dirname(module_root)
    if os.path.exists(os.path.join(module_root, "dune.module")):
        rel = os.path.relpath(script_dir, module_root)
        for build_dir in sorted(glob.glob(os.path.join(module_root, "build*"))):
            candidates.append(os.path.join(build_dir, rel))
    candidates = [os.getcwd(), script_dir] + candidates

    for directory in candidates:
        matches = sorted(glob.glob(os.path.join(directory, "test_wellbore_heat_transport*")))
        matches = [m for m in matches if os.path.isfile(m) and os.access(m, os.X_OK)]
        if matches:
            return os.path.abspath(matches[0])

    sys.exit("error: no test_wellbore_heat_transport* executable found, use --exe to point to it")


def run_simulation(executable, input_file, setup, level, refine, run_dir, extra_args, reuse):
    """Run one refinement level in its own directory, return the run directory."""
    if reuse and glob.glob(os.path.join(run_dir, "*.csv")):
        print(f"[level {level}] reusing existing results in {run_dir}")
        return

    if os.path.exists(run_dir):
        shutil.rmtree(run_dir)
    os.makedirs(run_dir)

    command = [executable, os.path.abspath(input_file)]
    # the level-0 soil grid comes from this script, not from params.input
    for direction, (cells, grading) in enumerate(zip(BASE_CELLS, GRADING)):
        command += [f"-Soil.Grid.Cells{direction}", cells]
        command += [f"-Soil.Grid.Grading{direction}", grading]
    if refine in ("both", "soil"):
        command += ["-Soil.Grid.Refinement", str(level)]
    if refine in ("both", "voids"):
        command += ["-Voids.Grid.Refinement", str(level)]
    # the grid file is given relative to the input file, make it absolute
    command += ["-Voids.Grid.File", setup["grid_file"]]
    command += extra_args

    print(f"[level {level}] {' '.join(command)}", flush=True)
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
    print(f"[level {level}] finished after {time.time() - start:.1f} s", flush=True)

    if process.returncode != 0:
        sys.exit(f"error: the simulation of level {level} failed, see {log_file}")


# ── Reading the simulation results ────────────────────────────────────────────
def load_outlet_temperature(run_dir):
    """Outlet temperature over time from the csv written by the 1D problem.

    The file name follows Vtk.OutputName of the input file, so results of an
    older run carry a different prefix — any csv in the run directory is used.
    """
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
    """Temperature profile along the borehole at the final time from the 1D vtp.

    Like the csv, the vtp files are taken from the run directory regardless of
    the prefix their name was written with (only the 1D domain writes vtp).
    """
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


# ── Results of one refinement level ───────────────────────────────────────────
def load_results(run_dir):
    """Outlet temperature over time and temperature profile of one level."""
    results = {}

    times, temperatures = load_outlet_temperature(run_dir)
    if times is not None:
        results["outlet_data"] = (times, temperatures)

    depth, profile = load_profile(run_dir)
    if depth is not None:
        results["profile_data"] = (depth, profile)

    return results


# ── Output ────────────────────────────────────────────────────────────────────
def plot_solutions(levels, results, setup, ramey, path):
    """Outlet temperature over time and temperature profile at the final time.

    Both subplots show the Ramey solution and, for the benchmark setup, the OGS
    reference, and carry the absolute error against Ramey on a second y axis
    (dashed), in the same style as analytical_solution.py.
    """
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 5.5))
    ax1e, ax2e = ax1.twinx(), ax2.twinx()
    # the OGS data belongs to the benchmark setup only
    show_ogs = (abs(setup["length"] - OGS_LENGTH) < 1e-6
                and abs(setup["t_end"] - OGS_T_END) < 1.0)

    times = np.linspace(0.0, setup["t_end"], 200)[1:]
    ax1.plot(times / 86400, ramey.outlet(times), color=REFERENCE_COLOR, linewidth=2,
             linestyle="--", label="Ramey (1962)")
    if show_ogs:
        mask = OGS_TIME > 0.0  # t = 0 is the initial condition, see load_outlet_temperature
        ax1.plot(OGS_TIME[mask] / 86400, OGS_OUTLET[mask], color=OGS_COLOR, linewidth=2,
                 label="OGS")
        ax1e.plot(OGS_TIME[mask] / 86400, OGS_OUTLET[mask] - ramey.outlet(OGS_TIME[mask]),
                  color=OGS_COLOR, linewidth=1.5, linestyle="--", label="error OGS")
    for i, level in enumerate(levels):
        if "outlet_data" not in results[i]:
            continue
        t, T = results[i]["outlet_data"]
        color = LEVEL_COLORS[i % len(LEVEL_COLORS)]
        ax1.plot(t / 86400, T, color=color, linewidth=2, label=f"refinement {level}")
        ax1e.plot(t / 86400, T - ramey.outlet(t), color=color, linewidth=1.5,
                  linestyle="--", label=f"error refinement {level}")
    ax1.set_xlabel("Time (d)")
    ax1.set_ylabel("Outlet temperature (°C)")
    ax1.set_title("Outlet temperature over time")

    depths = np.linspace(0.0, setup["length"], 200)
    ax2.plot(depths, ramey(depths, setup["t_end"]), color=REFERENCE_COLOR, linewidth=2,
             linestyle="--", label="Ramey (1962)")
    if show_ogs:
        ax2.plot(OGS_DEPTH, OGS_PROFILE, color=OGS_COLOR, linewidth=2, label="OGS")
        ax2e.plot(OGS_DEPTH, OGS_PROFILE - ramey(OGS_DEPTH, setup["t_end"]),
                  color=OGS_COLOR, linewidth=1.5, linestyle="--", label="error OGS")
    for i, level in enumerate(levels):
        if "profile_data" not in results[i]:
            continue
        x, T = results[i]["profile_data"]
        color = LEVEL_COLORS[i % len(LEVEL_COLORS)]
        ax2.plot(x, T, color=color, linewidth=2, label=f"refinement {level}")
        ax2e.plot(x, T - ramey(x, setup["t_end"]), color=color, linewidth=1.5,
                  linestyle="--", label=f"error refinement {level}")
    ax2.set_xlabel("Distance along borehole (m)")
    ax2.set_ylabel("Fluid temperature (°C)")
    ax2.set_xlim(0, setup["length"])
    ax2.set_title(f"Temperature profile at t = {setup['t_end'] / 86400:.3g} d")

    for ax, error_ax in ((ax1, ax1e), (ax2, ax2e)):
        ax.grid(color=GUIDE_COLOR, linewidth=0.5, alpha=0.4)
        ax.set_axisbelow(True)
        ax.spines["top"].set_visible(False)
        error_ax.set_ylabel("Absolute error vs Ramey (°C)")
        error_ax.spines["top"].set_visible(False)
        # one legend per subplot, holding the curves of both axes; it is put on
        # the twin because that one is drawn on top
        handles, labels = ax.get_legend_handles_labels()
        error_handles, error_labels = error_ax.get_legend_handles_labels()
        error_ax.legend(handles + error_handles, labels + error_labels,
                        frameon=False, fontsize=9)

    fig.tight_layout()
    fig.savefig(path, dpi=150)
    print(f"Saved: {path}")
    return fig


# ── Main ──────────────────────────────────────────────────────────────────────
def main():
    script_dir = os.path.dirname(os.path.abspath(__file__))
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--levels", type=int, nargs="+", default=[0, 1],
                        help="refinement levels to run (default: 0 1)")
    parser.add_argument("--exe", default=None,
                        help="path to the simulation executable (default: search build dirs)")
    parser.add_argument("--input", default=None,
                        help="input file (default: params.input next to the script or in the cwd)")
    parser.add_argument("--refine", choices=["both", "soil", "voids"], default="both",
                        help="which domain to refine (default: both)")
    parser.add_argument("--output-dir", default="convergence",
                        help="directory for the runs and the plots (default: convergence)")
    parser.add_argument("--reuse", action="store_true",
                        help="skip levels whose output already exists")
    parser.add_argument("--no-run", action="store_true",
                        help="only evaluate and plot existing results")
    parser.add_argument("--no-show", action="store_true",
                        help="save the plots without opening a window")
    args, extra_args = parser.parse_known_args()

    input_file = args.input
    if input_file is None:
        for candidate in ("params.input", os.path.join(script_dir, "params.input")):
            if os.path.exists(candidate):
                input_file = candidate
                break
        else:
            sys.exit("error: no params.input found, use --input to point to it")
    input_file = os.path.abspath(input_file)

    setup = collect_setup(input_file)
    ramey = RameySolution(setup)
    print(f"Input file: {input_file}")
    print(f"Borehole: L = {setup['length']} m, r_inner = {setup['r_inner']} m, "
          f"r_outer = {setup['r_outer']} m, {setup['num_segments']} elements at level 0")
    print(f"Soil grid at level 0: Cells = {' | '.join(BASE_CELLS)}, "
          f"Grading = {' | '.join(GRADING)} ({num_soil_cells(0)} cells)")
    print(f"Ramey: Re = {ramey.Re:.1f}, Pr = {ramey.Pr:.3f}, Nu = {ramey.Nu:.4f}, "
          f"h = {ramey.h:.4f} W/m²K, U = {ramey.U:.4f} W/m²K")

    output_dir = os.path.abspath(args.output_dir)
    os.makedirs(output_dir, exist_ok=True)

    levels = sorted(set(args.levels))
    executable = None if args.no_run else find_executable(args.exe, script_dir)

    results = []
    for level in levels:
        run_dir = os.path.join(output_dir, f"refinement_{level}")
        if not args.no_run:
            run_simulation(executable, input_file, setup, level, args.refine,
                           run_dir, extra_args, args.reuse)
        level_results = load_results(run_dir)
        if not level_results:
            sys.exit(f"error: no results found for level {level} in {run_dir}")
        results.append(level_results)

    if not any("profile_data" in r for r in results):
        print("\nnote: python3-vtk not available or no vtp output — "
              "the temperature profile along the borehole is skipped")

    plot_solutions(levels, results, setup, ramey,
                   os.path.join(output_dir, "convergence_temperatures.png"))

    if not args.no_show:
        plt.show()


if __name__ == "__main__":
    main()
