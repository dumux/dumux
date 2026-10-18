#!/usr/bin/env python3
# SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

"""
Elder benchmark: grid and time-step study of the low-Rayleigh Elder problem (Ra = 60) against the pseudo-spectral
reference solution, the solution without buoyancy (Ra = 0) against the exact solution, the sharpness of the
source boundary data, and the steady states at Ra = 400.

Usage: python3 run_benchmark.py ./test_1p2c_elder_box [output directory of the figures] [--jobs N] [--max-cells N]

Each run is stored as <name>.json with the commit and a hash of the uncommitted changes of the source tree, the
hashes of the executable and of params.input, the arguments, the return code and the results. A run is reused
when the executable, params.input and the arguments are unchanged. Pseudo-spectral reference solutions are cached
in reference_*.npz files.
"""

import hashlib
import io
import json
import math
import os
import re
import subprocess
import sys
import time
import xml.etree.ElementTree as ET
from concurrent.futures import ThreadPoolExecutor

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

plt.rcParams.update(
    {
        "font.size": 13,
        "axes.titlesize": 13,
        "legend.fontsize": 12,
        "xtick.labelsize": 12,
        "ytick.labelsize": 12,
    }
)
import numpy as np
from PIL import Image

SOURCE_DIR = os.path.dirname(os.path.realpath(__file__))
sys.path.insert(0, SOURCE_DIR)
sys.dont_write_bytecode = True  # keep the source directory clean
from elder_pseudospectral import ElderPseudoSpectral  # noqa: E402

HEIGHT, LENGTH = 150.0, 300.0  # [m], half domain
TIMES = [0.03, 0.1, 0.3, 1.0]  # dimensionless output times
DELTA = 0.05  # transition length of the smoothed source in units of the height
REFERENCE_RESOLUTION = (129, 41, 2.5e-5)  # collocation points in x and z, time step
SPATIAL_TIME_STEP = 1.25e-3  # DIRK3 time step of the grid studies

NUMBER = r"([-+0-9.eE]+|nan|-?inf)"
RESULT = re.compile(
    rf"^t = {NUMBER}: mass {NUMBER}, Sherwood {NUMBER}, bottom {NUMBER}"
    rf"(?:, relative L2 error without buoyancy {NUMBER})?",
    re.M,
)


def sha256_of_file(filename):
    with open(filename, "rb") as f:
        return hashlib.sha256(f.read()).hexdigest()


def source_state():
    """Commit and hash of the uncommitted changes (tracked and untracked files) of the source tree"""

    def git(*args):
        return subprocess.run(
            ["git", "-C", SOURCE_DIR, *args], capture_output=True, text=True, check=True
        ).stdout

    changes = hashlib.sha256(git("diff", "HEAD").encode())
    for name in sorted(git("ls-files", "--others", "--exclude-standard").split()):
        changes.update(name.encode())
        changes.update(open(os.path.join(SOURCE_DIR, name), "rb").read())
    return {"commit": git("rev-parse", "HEAD").strip(), "uncommitted": changes.hexdigest()}


class Runner:
    def __init__(self, executable, jobs):
        self.executable = os.path.abspath(executable)
        sources = [os.path.join(SOURCE_DIR, f) for f in ("main.cc", "problem.hh", "properties.hh")]
        if os.path.getmtime(self.executable) < max(os.path.getmtime(f) for f in sources):
            sys.exit(f"{self.executable} is older than its sources, rebuild it")
        self.state = source_state()
        self.identity = {
            "executable": sha256_of_file(self.executable),
            "params": sha256_of_file("params.input"),
        }
        self.jobs = jobs
        self.env = dict(
            os.environ,
            DUMUX_NUM_THREADS="1",
            VECLIB_MAXIMUM_THREADS="1",
            OMP_NUM_THREADS="1",
            OPENBLAS_NUM_THREADS="1",
        )

    def _run(self, name, args):
        record_file = f"{name}.json"
        if os.path.exists(record_file):
            record = json.load(open(record_file))
            if (
                record.get("identity") == self.identity
                and record.get("args") == args
                and record.get("returncode") == 0
                and complete(record)
            ):
                return record
        start = time.time()
        process = subprocess.run(
            [
                self.executable,
                "params.input",
                "-Problem.Name",
                name,
                "-Vtk.Precision",
                "Float64",
                *args,
            ],
            capture_output=True,
            text=True,
            env=self.env,
        )
        results = [
            [float(v) if v is not None else None for v in m.groups()]
            for m in RESULT.finditer(process.stdout)
        ]
        record = {
            "name": name,
            "args": args,
            "identity": self.identity,
            "source": self.state,
            "returncode": process.returncode,
            "seconds": time.time() - start,
            "results": results,
        }
        with open(record_file, "w") as f:
            json.dump(record, f, indent=1)
        with open(f"{name}.log", "w") as f:
            f.write(process.stdout[-200000:] + process.stderr)
        print(f"{name}: return code {process.returncode}, {record['seconds']:.0f} s", flush=True)
        return record

    def run(self, cases):
        """Run the cases {name: args} in parallel, the most expensive first, and return the checked records"""

        def cost(item):
            args = item[1]
            nx, nz = map(int, args[args.index("-Grid.Cells") + 1].split())
            dt = float(args[args.index("-TimeLoop.DimensionlessTimeStepSize") + 1])
            return nx * nz / dt

        ordered = sorted(cases.items(), key=cost, reverse=True)
        with ThreadPoolExecutor(self.jobs) as pool:
            records = dict(
                zip([name for name, _ in ordered], pool.map(lambda item: self._run(*item), ordered))
            )
        for name, record in records.items():
            check_record(record, self.identity)
        return records


def output_times(args):
    return [float(t) for t in args[args.index("-TimeLoop.DimensionlessOutputTimes") + 1].split()]


def complete(record):
    """One result and one VTU file per output time"""
    n = len(output_times(record["args"]))
    return len(record["results"]) == n and all(
        os.path.exists(f"{record['name']}-{i:05d}.vtu") for i in range(1, n + 1)
    )


def check_record(record, identity):
    """Refuse runs without source record, failed runs, incomplete output and non-finite results"""
    for key in ("source", "identity", "args", "returncode", "results"):
        if key not in record:
            sys.exit(f"{record.get('name')}: incomplete run record")
    if record["identity"] != identity:
        sys.exit(f"{record['name']}: produced by another executable or input file")
    if record["returncode"] != 0:
        sys.exit(
            f"{record['name']}: the run failed (Newton's method did not converge or an assertion failed)"
        )
    values = [v for row in record["results"] for v in row if v is not None]
    if not complete(record) or not all(math.isfinite(v) for v in values):
        sys.exit(f"{record['name']}: incomplete or non-finite results")


def reference(rayleigh, delta, times, resolution=REFERENCE_RESOLUTION, initial=None, label=""):
    """Pseudo-spectral solution at the given times, cached with the hash of the solver"""
    nx, nz, dt = resolution
    filename = f"reference_ra{rayleigh:g}_delta{delta:g}_{nx}x{nz}_dt{dt:g}{label}.npz"
    solver = ElderPseudoSpectral(rayleigh, delta, nx, nz)
    solver_hash = sha256_of_file(os.path.join(SOURCE_DIR, "elder_pseudospectral.py"))
    if os.path.exists(filename):
        data = np.load(filename)
        if str(data["solver_hash"]) == solver_hash and np.allclose(data["times"], times):
            return solver, data["snapshots"]
    print(f"computing the pseudo-spectral reference {filename}", flush=True)
    snapshots, _ = solver.run(times, dt, initial=initial(solver) if initial else None)
    snapshots = np.array(snapshots)
    if not np.all(np.isfinite(snapshots)):
        sys.exit(f"non-finite pseudo-spectral solution {filename}")
    np.savez(
        filename,
        snapshots=snapshots,
        times=np.array(times),
        solver_hash=solver_hash,
        rayleigh=rayleigh,
        delta=delta,
        nx=nx,
        nz=nz,
        dt=dt,
    )
    return solver, snapshots


def read_vtu(filename):
    """Points in units of the height and point data of an ASCII VTU file"""
    piece = ET.parse(filename).getroot().find("UnstructuredGrid/Piece")

    def to_array(array):
        return np.array(array.text.split(), dtype=float).reshape(
            -1, int(array.get("NumberOfComponents", "1"))
        )

    points = to_array(piece.find("Points/DataArray"))[:, :2] / HEIGHT
    return points, {array.get("Name"): to_array(array) for array in piece.find("PointData")}


class StructuredField:
    """Vertex values of a VTU output of a structured grid with the nodal quadrature weights of the box scheme"""

    def __init__(self, filename):
        points, data = read_vtu(filename)
        self.x, ix = np.unique(np.round(points[:, 0], 10), return_inverse=True)
        self.z, iz = np.unique(np.round(points[:, 1], 10), return_inverse=True)
        self.c = np.empty((len(self.x), len(self.z)))
        self.c[ix, iz] = data["X^Solute_liq"][:, 0]
        self.w = np.empty_like(self.c)
        self.w[ix, iz] = data["velocity_liq (m/s)"][:, 1]
        hx, hz = np.gradient(self.x), np.gradient(self.z)
        hx[[0, -1]] *= 0.5
        hz[[0, -1]] *= 0.5
        self.weights = np.outer(hx, hz)

    def relative_l2(self, difference, reference):
        """L2 norm of the difference relative to that of the reference, with the control volumes as quadrature"""
        return math.sqrt(np.sum(difference**2 * self.weights) / np.sum(reference**2 * self.weights))


def field(name, index):
    return StructuredField(f"{name}-{index + 1:05d}.vtu")


def count_plumes(vertical_velocity_at_mid_depth):
    """Number of downward plumes in the full domain obtained by mirroring at the symmetry axis"""
    w = np.concatenate([vertical_velocity_at_mid_depth[::-1], vertical_velocity_at_mid_depth[1:]])
    down = w < -1e-3 * np.abs(w).max()
    return int(np.sum(down[1:] & ~down[:-1]) + down[0])


def rates(errors):
    return ["-"] + [f"{math.log2(a / b):.2f}" for a, b in zip(errors[:-1], errors[1:])]


def sci(value, digits=2):
    """Number in the notation of the README tables"""
    exponent = int(math.floor(math.log10(abs(value))))
    return f"${{{value / 10**exponent:.{digits}f} \\times 10^{{{exponent}}}}}$"


def exact(value):
    """Exact decimal value, such as a time step, in the notation of the README tables"""
    exponent = int(math.floor(math.log10(abs(value)) + 1e-12))
    mantissa = f"{value / 10**exponent:.6g}"
    return (
        f"$10^{{{exponent}}}$" if mantissa == "1" else f"${{{mantissa} \\times 10^{{{exponent}}}}}$"
    )


def with_rates(errors):
    return [sci(e) + ("" if r == "-" else f" ({r})") for e, r in zip(errors, rates(errors))]


def grid_label(n):
    return f"${{{n} \\times {n // 2}}}$"


def print_table(header, rows):
    print("| " + " | ".join(header) + " |")
    print("|" + "------|" * len(header))
    for row in rows:
        print("| " + " | ".join(row) + " |")


def grid_args(n, *extra):
    return ["-Grid.Cells", f"{n} {n // 2}", *extra]


def main():
    args = [a for a in sys.argv[1:]]
    jobs, max_cells = 8, 512
    for option in ("--jobs", "--max-cells"):
        if option in args:
            i = args.index(option)
            value = int(args[i + 1])
            del args[i : i + 2]
            jobs, max_cells = (value, max_cells) if option == "--jobs" else (jobs, value)
    if not args:
        sys.exit(__doc__)
    try:
        lock = os.open("run_benchmark.lock", os.O_CREAT | os.O_EXCL)
    except FileExistsError:
        sys.exit("run_benchmark.lock exists: another run_benchmark.py works in this directory")
    try:
        run_studies(args, jobs, max_cells)
    finally:
        os.close(lock)
        os.remove("run_benchmark.lock")


def run_studies(args, jobs, max_cells):
    runner = Runner(args[0], jobs)
    image_dir = args[1] if len(args) > 1 else "."
    grids = [n for n in (64, 128, 256, 512) if n <= max_cells]
    dirk = [
        "-TimeLoop.Method",
        "DIRK3",
        "-TimeLoop.DimensionlessTimeStepSize",
        f"{SPATIAL_TIME_STEP}",
    ]
    times_arg = ["-TimeLoop.DimensionlessOutputTimes", " ".join(map(str, TIMES))]

    # all DuMux runs
    cases = {}
    for weight in (1.0, 0.5):
        for n in grids:
            cases[f"elder_grid_{n}_w{weight:g}"] = grid_args(
                n, "-Flux.UpwindWeight", f"{weight}", *dirk, *times_arg
            )
    for n in grids:
        cases[f"elder_nobuoyancy_{n}"] = grid_args(
            n, "-Problem.DensityDifference", "0", *dirk, *times_arg
        )
        cases[f"elder_delta0.005_{n}"] = grid_args(
            n,
            "-Flux.UpwindWeight",
            "0.5",
            "-Problem.SourceTransitionLength",
            f"{0.005 * HEIGHT}",
            *dirk,
            *times_arg,
        )
        cases[f"elder_sharp_{n}"] = grid_args(
            n,
            "-Flux.UpwindWeight",
            "0.5",
            "-Problem.SourceTransitionLength",
            "0",
            *dirk,
            *times_arg,
        )
    time_steps = {
        "ImplicitEuler": [1e-2, 5e-3, 2.5e-3, 1.25e-3, 6.25e-4, 3.125e-4],
        "DIRK3": [1e-2, 5e-3, 2.5e-3, 1.25e-3],
    }
    for method, steps in time_steps.items():
        for dt in steps:
            cases[f"elder_time_{method}_{dt:g}"] = grid_args(
                64,
                "-Flux.UpwindWeight",
                "0.5",
                "-TimeLoop.Method",
                method,
                "-TimeLoop.DimensionlessTimeStepSize",
                f"{dt}",
                *times_arg,
            )
    cases["elder_time_reference"] = grid_args(
        64,
        "-Flux.UpwindWeight",
        "0.5",
        "-TimeLoop.Method",
        "DIRK3",
        "-TimeLoop.DimensionlessTimeStepSize",
        "3.125e-4",
        *times_arg,
    )
    ra400_times = [0.1, 0.3, 0.6, 1.0]
    ra400 = [(n, dt) for n in (64, 128, 256) if n <= max_cells for dt in (2.5e-3, 1.25e-3)]
    for n, dt in ra400:
        cases[f"elder_ra400_{n}_dt{dt:g}"] = grid_args(
            n,
            "-Problem.DensityDifference",
            "200",
            "-Flux.UpwindWeight",
            "1",
            "-TimeLoop.Method",
            "DIRK3",
            "-TimeLoop.DimensionlessTimeStepSize",
            f"{dt}",
            "-TimeLoop.DimensionlessOutputTimes",
            " ".join(map(str, ra400_times)),
        )
    records = runner.run(cases)

    # pseudo-spectral references
    solver, snapshots = reference(60.0, DELTA, TIMES)
    diag = [solver.diagnostics(a) for a in snapshots]

    def reference_error(name, index):
        f = field(name, index)
        reference_values = solver.evaluate(snapshots[index], f.x, f.z)
        return f.relative_l2(f.c - reference_values, reference_values)

    time_header = [f"$\\varepsilon$, $t = {t:g}$" for t in TIMES]
    print("\nReference (pseudo-spectral, delta = 0.05, 129 x 41, dt = 2.5e-5)")
    print_table(
        ["$t$", "$M$", "$Sh$", "$Sh_b$"],
        [
            [f"{t:g}", f"{d['mass']:.8f}", f"{d['sherwood']:.8f}", f"{d['bottom']:.8f}"]
            for t, d in zip(TIMES, diag)
        ],
    )

    # grid study at Ra = 60
    grid_errors = {}
    for weight in (0.5, 1.0):
        errors = np.array(
            [
                [reference_error(f"elder_grid_{n}_w{weight:g}", i) for i in range(len(TIMES))]
                for n in grids
            ]
        )
        grid_errors[weight] = errors
        print(
            f"\nGrid study, Ra = 60, upwind weight {weight:g}: relative L2 error of X against the reference (rate)"
        )
        columns = [with_rates(errors[:, i]) for i in range(len(TIMES))]
        print_table(
            ["Grid", *time_header],
            [[grid_label(n), *(c[k] for c in columns)] for k, n in enumerate(grids)],
        )
        for quantity, column, key in (
            ("solute mass $M$", 1, "mass"),
            ("Sherwood number $Sh$", 2, "sherwood"),
        ):
            print(f"\nRelative error of the {quantity}, upwind weight {weight:g}")
            rows = []
            for n in grids:
                res = records[f"elder_grid_{n}_w{weight:g}"]["results"]
                rows.append(
                    [
                        grid_label(n),
                        *(
                            sci(abs(res[i][column] - d[key]) / d[key], 1)
                            for i, d in enumerate(diag)
                        ),
                    ]
                )
            print_table(["Grid", *(f"$t = {t:g}$" for t in TIMES)], rows)

    # time-step study on 64 x 32
    print(
        "\nTime-step study on 64 x 32 cells: relative L2 difference of X to DIRK3 with dt = 3.125e-4, "
        "maximum over the output times"
    )
    time_errors = {}
    ref_fields = [field("elder_time_reference", i) for i in range(len(TIMES))]
    for method, steps in time_steps.items():
        errors = [
            max(
                f.relative_l2(field(f"elder_time_{method}_{dt:g}", i).c - f.c, f.c)
                for i, f in enumerate(ref_fields)
            )
            for dt in steps
        ]
        time_errors[method] = (steps, errors)
    all_steps = time_steps["ImplicitEuler"]
    columns = {
        method: dict(zip(steps, with_rates(time_errors[method][1])))
        for method, steps in time_steps.items()
    }
    print_table(
        ["$\\Delta t$", "implicit Euler", "DIRK3"],
        [
            [exact(dt), columns["ImplicitEuler"].get(dt, ""), columns["DIRK3"].get(dt, "")]
            for dt in all_steps
        ],
    )

    # discrete solute balance of implicit Euler, M^{n+1} - M^n = dt L/H (Sh - Sh_b)^{n+1}, with the boundary fluxes
    # from the balance of the boundary control volumes
    history = np.genfromtxt(
        "elder_time_ImplicitEuler_0.01_diagnostics.csv", delimiter=",", names=True
    )
    balance = (
        np.diff(history["mass"])
        - np.diff(history["time"]) * LENGTH / HEIGHT * (history["sherwood"] - history["bottom"])[1:]
    )
    print(
        f"\nSolute balance of implicit Euler (64 x 32, dt = 0.01): max |M^(n+1) - M^n - dt L/H (Sh - Sh_b)| = "
        f"{np.abs(balance).max():.1e}, relative to M(1) = {history['mass'][-1]:.4f}"
    )

    # Ra = 0 against the exact solution
    print("\nRa = 0: relative L2 error against the exact solution (rate)")
    nobuoyancy_errors = np.array(
        [[row[4] for row in records[f"elder_nobuoyancy_{n}"]["results"]] for n in grids]
    )
    columns = [with_rates(nobuoyancy_errors[:, i]) for i in range(len(TIMES))]
    print_table(
        ["Grid", *time_header],
        [[grid_label(n), *(c[k] for c in columns)] for k, n in enumerate(grids)],
    )

    # sharpness of the source boundary data: differences to the next finer grid at the common vertices
    variants = {
        "$\\delta = 0.05 H$": "elder_grid_{}_w0.5",
        "$\\delta = 0.005 H$": "elder_delta0.005_{}",
        "step": "elder_sharp_{}",
    }
    differences = {}
    for label, case in variants.items():
        l2, maximum = [], []
        for coarse, fine in zip(grids[:-1], grids[1:]):
            fc, ff = field(case.format(coarse), 2), field(case.format(fine), 2)
            difference = fc.c - ff.c[::2, ::2]
            l2.append(fc.relative_l2(difference, ff.c[::2, ::2]))
            maximum.append(np.abs(difference).max())
        differences[label] = (l2, maximum)
    for norm, k in (("relative L2", 0), ("maximum", 1)):
        print(
            f"\nSharpness of the source, t = 0.3, upwind weight 0.5: {norm} difference of X "
            "to the next finer grid (rate)"
        )
        columns = [with_rates(differences[label][k]) for label in variants]
        print_table(
            ["Grids", *variants],
            [
                [f"{grid_label(coarse)}, {grid_label(fine)}", *(c[m] for c in columns)]
                for m, (coarse, fine) in enumerate(zip(grids[:-1], grids[1:]))
            ],
        )

    # Ra = 400: steady states of the reference and DuMux
    def initial_s1(s):
        x, z = np.meshgrid(s.x, s.z[1:-1], indexing="ij")
        return s.cos_coefficients(np.where(x < 1.0, (1.0 + np.cos(np.pi * x)) / 2.0, 0.0))

    steady_times = [2.0]
    solver_400, s2 = reference(
        400.0, DELTA, steady_times, resolution=(129, 41, 5e-5), label="_zero"
    )
    _, s1 = reference(
        400.0, DELTA, steady_times, resolution=(129, 41, 5e-5), initial=initial_s1, label="_s1"
    )
    states = {"S1": s1[-1], "S2": s2[-1]}
    words = {1: "one plume", 2: "two plumes", 3: "three plumes"}
    print("\nRa = 400, steady states of the pseudo-spectral solution (129 x 41)")
    xs = np.linspace(0.0, 2.0, 801)
    for label, a in states.items():
        _, qz = solver_400.evaluate_velocity(a, xs, [0.5])
        print(
            f"{label}: {words.get(count_plumes(qz[:, 0]))}, Sh = {solver_400.diagnostics(a)['sherwood']:.4f}"
        )
    print("\nRa = 400, DuMux states at t = 1 from X = 0, full upwinding, DIRK3")
    plumes = {}
    for n, dt in ra400:
        f = field(f"elder_ra400_{n}_dt{dt:g}", len(ra400_times) - 1)
        plumes[(n, dt)] = count_plumes(f.w[:, np.argmin(np.abs(f.z - 0.5))])
    steps = sorted({dt for _, dt in ra400}, reverse=True)
    print_table(
        ["Grid", *(f"$\\Delta t = {exact(dt)[1:-1]}$" for dt in steps)],
        [
            [
                grid_label(n),
                *(
                    f"{words.get(plumes[(n, dt)])}, $Sh = {records[f'elder_ra400_{n}_dt{dt:g}']['results'][-1][2]:.4f}$"
                    for dt in steps
                ),
            ]
            for n in sorted({n for n, _ in ra400})
        ],
    )

    # figures
    plot_fields(
        image_dir,
        f"elder_grid_{grids[-2] if len(grids) > 1 else grids[-1]}_w0.5",
        solver,
        snapshots,
    )
    plot_convergence(image_dir, grids, grid_errors, nobuoyancy_errors, time_errors)
    plot_ra400(image_dir, ra400, plumes, solver_400, states)


def save_figure(fig, filename):
    """PNG with an adaptive palette of 256 colors, which keeps the file small"""
    buffer = io.BytesIO()
    fig.savefig(buffer, format="png", dpi=150)
    buffer.seek(0)
    Image.open(buffer).convert("RGB").quantize(colors=256).save(filename, optimize=True)


def plot_fields(image_dir, name, solver, snapshots):
    """DuMux mass fraction and isolines compared with the isolines of the pseudo-spectral solution"""
    fig, axes = plt.subplots(
        2, 2, figsize=(12, 6.6), layout="constrained", sharex=True, sharey=True
    )
    levels = np.linspace(0.1, 0.9, 9)
    for ax, i, t in zip(axes.flat, range(len(TIMES)), TIMES):
        f = field(name, i)
        image = ax.contourf(f.x, f.z, f.c.T, levels=np.linspace(0.0, 1.0, 21), cmap="viridis")
        ax.contour(f.x, f.z, f.c.T, levels=levels, colors="w", linewidths=1.6)
        xs, zs = np.linspace(0.0, 2.0, 401), np.linspace(0.0, 1.0, 201)
        ax.contour(
            xs,
            zs,
            solver.evaluate(snapshots[i], xs, zs).T,
            levels=levels,
            colors="k",
            linewidths=0.7,
            linestyles="--",
        )
        ax.set_aspect("equal")
        ax.set_xticks([0.0, 0.5, 1.0, 1.5, 2.0])
        ax.set_title(f"$t = {t:g}$ ({t * 200:g} years)")
    for ax in axes[-1]:
        ax.set_xlabel("$x / H$")
    for ax in axes[:, 0]:
        ax.set_ylabel("$z / H$")
    fig.colorbar(image, ax=axes, label="solute mass fraction $X$", ticks=np.linspace(0, 1, 6))
    save_figure(fig, f"{image_dir}/elder_concentration.png")


def slope_triangle(ax, x0, y0, order, corner_below, factor=2.0):
    """Triangle with the slope `order` in a log-log plot whose horizontal leg spans x0 to factor x0"""
    x1, y1 = x0 * factor, y0 * factor**order
    corner = (x1, y0) if (order > 0) == corner_below else (x0, y1)
    ax.plot([x0, x1, corner[0], x0], [y0, y1, corner[1], y0], color="k", lw=0.8)
    ax.text(
        math.sqrt(x0 * x1),
        corner[1] / 1.4 if corner_below else corner[1] * 1.4,
        "1",
        ha="center",
        va="center",
        fontsize=12,
    )
    right = corner[0] == x1
    ax.text(
        corner[0] * (1.08 if right else 1 / 1.08),
        math.sqrt(y0 * y1),
        f"{abs(order):g}",
        ha="left" if right else "right",
        va="center",
        fontsize=12,
    )


def plot_convergence(image_dir, grids, grid_errors, nobuoyancy_errors, time_errors):
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5.6), layout="constrained")
    h = np.array([n // 2 for n in grids])
    markers = {0.03: "o", 0.1: "s", 0.3: "^", 1.0: "D"}
    for i, t in enumerate(TIMES):
        ax1.loglog(
            h,
            grid_errors[0.5][:, i],
            marker=markers[t],
            color=f"C{i}",
            label=f"$t = {t:g}$, weight 0.5",
        )
        ax1.loglog(
            h,
            grid_errors[1.0][:, i],
            marker=markers[t],
            color=f"C{i}",
            ls=":",
            mfc="none",
            label=f"$t = {t:g}$, weight 1",
        )
    ax1.loglog(h, nobuoyancy_errors[:, 2], "k-", marker="x", label="$Ra = 0$, $t = 0.3$")
    k2 = max(len(h) - 2, 0)
    slope_triangle(ax1, h[k2], 0.4 * nobuoyancy_errors[k2, 2], -2, corner_below=True)
    k1 = max(len(h) - 3, 0)
    slope_triangle(ax1, h[k1], 1.6 * grid_errors[1.0][k1, 2], -1, corner_below=False)
    ax1.set_xticks(h, [str(n) for n in h])
    ax1.minorticks_off()
    ax1.set_xlabel("cells per height $H$")
    ax1.set_ylabel("relative $L^2$ error of $X$")
    fig.legend(*ax1.get_legend_handles_labels(), loc="outside lower left", ncol=5)
    for method, color, style in (("ImplicitEuler", "C0", "o-"), ("DIRK3", "C1", "s--")):
        steps, errors = time_errors[method]
        ax2.loglog(
            steps,
            errors,
            style,
            color=color,
            label="implicit Euler" if method == "ImplicitEuler" else "DIRK3",
        )
    # the triangles show the classical orders of the methods
    for method, order, scale, below in (("ImplicitEuler", 1, 0.4, True), ("DIRK3", 3, 3.0, False)):
        steps, errors = time_errors[method]
        j = min(len(steps) - 2, 2)
        slope_triangle(ax2, steps[j + 1], scale * errors[j + 1], order, corner_below=below)
    ax2.set_xlabel(r"time step $\Delta t$ $[\phi H^2 / D_\mathrm{eff}]$")
    ax2.set_ylabel("relative $L^2$ time error ($64 \\times 32$ cells)")
    ax2.legend()
    save_figure(fig, f"{image_dir}/elder_convergence.png")


def plot_ra400(image_dir, runs, plumes, solver, states):
    """Steady states of the pseudo-spectral solution and one DuMux state with one and with two plumes"""
    fig, axes = plt.subplots(
        2, 2, figsize=(12, 6.6), layout="constrained", sharex=True, sharey=True
    )
    xs, zs = np.linspace(0.0, 2.0, 401), np.linspace(0.0, 1.0, 201)
    panels = [
        (f"pseudo-spectral {k}", lambda a=a: (xs, zs, solver.evaluate(a, xs, zs)))
        for k, a in states.items()
    ]
    for count in (1, 2):
        matches = [run for run in runs if plumes[run] == count]
        if matches:
            n, dt = max(matches)
            panels.append(
                (
                    f"DuMux ${n} \\times {n // 2}$, $\\Delta t = {dt:g}$, $t = 1$",
                    lambda n=n, dt=dt: (lambda f: (f.x, f.z, f.c))(
                        field(f"elder_ra400_{n}_dt{dt:g}", 3)
                    ),
                )
            )
    for ax, (title, data) in zip(axes.flat, panels):
        x, z, c = data()
        image = ax.contourf(x, z, c.T, levels=np.linspace(0.0, 1.0, 21), cmap="viridis")
        ax.contour(x, z, c.T, levels=np.linspace(0.1, 0.9, 9), colors="w", linewidths=0.8)
        ax.set_aspect("equal")
        ax.set_xticks([0.0, 0.5, 1.0, 1.5, 2.0])
        ax.set_title(title)
    for ax in axes[-1]:
        ax.set_xlabel("$x / H$")
    for ax in axes[:, 0]:
        ax.set_ylabel("$z / H$")
    fig.colorbar(image, ax=axes, label="solute mass fraction $X$", ticks=np.linspace(0, 1, 6))
    save_figure(fig, f"{image_dir}/elder_ra400.png")


if __name__ == "__main__":
    main()
