#!/usr/bin/env python3
# SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later
"""
Run the Taylor-Green vortex benchmark: spatial and temporal convergence studies
and the kinetic energy decay, see README.md.

The script has to be executed from the build directory of this test,
where the executables and parameter files are located.
"""

import argparse
import csv
import math
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path


SCHEMES = ("pq1bubble", "pq1bubblehybrid", "pq2hybrid")
TIME_SCHEMES = ("ImplicitEuler", "CrankNicolson", "DIRK3")

# expected convergence orders under grid refinement (velocity L2, velocity H1, pressure L2)
EXPECTED_SPATIAL_ORDER = {
    "pq1bubble": {"velocityL2": 2.0, "velocityH1": 1.0, "pressureL2": 2.0, "pressureH1": 1.0},
    "pq1bubblehybrid": {"velocityL2": 2.0, "velocityH1": 1.0, "pressureL2": 2.0, "pressureH1": 1.0},
    "pq2hybrid": {"velocityL2": 3.0, "velocityH1": 2.0, "pressureL2": 2.0, "pressureH1": 1.0},
}
EXPECTED_TEMPORAL_ORDER = {"ImplicitEuler": 1.0, "CrankNicolson": 2.0, "DIRK3": 3.0}
ORDER_TOLERANCE = 0.3

# number of cells per direction on the coarsest grid (coarser grids do not resolve
# the vortex well enough for the stationary Navier-Stokes problem to converge)
BASE_CELLS = {2: 16, 3: 6}
# default number of refinement levels (full study, test mode)
DEFAULT_LEVELS = {2: (4, 3), 3: (3, 2)}
ERROR_KEYS = ("velocityL2", "velocityH1", "pressureL2", "pressureH1")

# plot labels (notation as in README.md)
SCHEME_LABELS = {
    "pq1bubble": "PQ1Bubble",
    "pq1bubblehybrid": "hybrid PQ1Bubble",
    "pq2hybrid": "hybrid PQ2",
}
TIME_SCHEME_LABELS = {
    "ImplicitEuler": "implicit Euler",
    "CrankNicolson": "Crank-Nicolson",
    "DIRK3": "DIRK3",
}
ERROR_LABELS = {
    "velocityL2": r"$\|\mathbf{u} - \mathbf{u}_h\|_{L^2(\Omega)}$",
    "velocityH1": r"$\|\mathbf{u} - \mathbf{u}_h\|_{H^1(\Omega)}$",
    "pressureL2": r"$\|p - p_h\|_{L^2(\Omega)}$",
    "pressureH1": r"$\|p - p_h\|_{H^1(\Omega)}$",
}


@dataclass
class Run:
    """Result of a single simulation run (last row of the error file)"""

    label: str
    h: float
    dt: float
    errors: dict
    history: list


def executable(dim, scheme, multistage=False):
    """Name of the executable for the given dimension and momentum scheme"""
    suffix = "_multistage" if multistage else ""
    return f"test_ff_navierstokes_taylorgreen_{dim}d_{scheme}{suffix}"


def run_simulation(exe, dim, name, args):
    """Run a simulation and return the rows of the error file"""
    if not Path(exe).is_file():
        raise FileNotFoundError(f"Executable {exe} not found. Build it first.")

    command = [f"./{exe}", f"params_{dim}d.input", "-Problem.Name", name]
    command += ["-Problem.EnableVtkOutput", "false", "-Problem.EnableGravity", "false"] + args
    print("+ " + " ".join(command), flush=True)
    subprocess.run(command, check=True, stdout=subprocess.DEVNULL)

    with open(f"{name}_errors.csv", newline="") as errorFile:
        reader = csv.DictReader(errorFile)
        rows = [{key: float(value) for key, value in row.items()} for row in reader]
    if not rows:
        raise RuntimeError(f"No errors written by {exe}")
    return rows


def to_run(label, rows):
    """Convert the rows of an error file to a run"""
    last = rows[-1]
    return Run(
        label=label,
        h=last["h"],
        dt=last["dt"],
        errors={key: last[key] for key in ERROR_KEYS},
        history=rows,
    )


def rates(runs, key, step):
    """Experimental orders of convergence between consecutive runs"""
    result = []
    for coarse, fine in zip(runs[:-1], runs[1:]):
        errorRatio = coarse.errors[key] / fine.errors[key]
        stepRatio = step(coarse) / step(fine)
        result.append(math.log(errorRatio) / math.log(stepRatio))
    return result


def write_table(fileName, title, runs, stepName, step):
    """Write a Markdown table with errors and convergence rates"""
    rateTable = {key: [float("nan")] + rates(runs, key, step) for key in ERROR_KEYS}
    header = f"| {stepName} | " + " | ".join(f"{key} | EOC" for key in ERROR_KEYS) + " |"
    lines = [f"### {title}", "", header, "|" + "---|" * (1 + 2 * len(ERROR_KEYS))]
    for i, run in enumerate(runs):
        cells = [f"{step(run):.4e}"]
        for key in ERROR_KEYS:
            cells += [f"{run.errors[key]:.4e}", f"{rateTable[key][i]:.2f}"]
        lines.append("| " + " | ".join(cells) + " |")
    lines.append("")

    text = "\n".join(lines)
    print(text)
    with open(fileName, "a") as tableFile:
        tableFile.write(text + "\n")


def import_pyplot():
    """Import matplotlib (if available) with a LaTeX-like style"""
    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        print("matplotlib not available, skipping plots")
        return None

    plt.rcParams.update(
        {
            "font.family": "serif",
            "font.serif": ["cmr10", "DejaVu Serif"],
            "mathtext.fontset": "cm",
            "axes.formatter.use_mathtext": True,
            "font.size": 11,
            "axes.labelsize": 12,
            "legend.fontsize": 9,
            "lines.markersize": 5,
        }
    )
    return plt


def format_order(order):
    """Format a convergence order for a LaTeX exponent"""
    return f"{order:g}"


def add_reference_slope(ax, steps, errors, order, stepSymbol, color):
    """Add a dashed line with the expected slope through the finest data point"""
    stepRange = [min(steps), max(steps)]
    finest = steps.index(min(steps))
    reference = [errors[finest] * (h / steps[finest]) ** order for h in stepRange]
    ax.loglog(stepRange, reference, "--", color=color, linewidth=1.0, alpha=0.8)
    ax.annotate(
        rf"$\mathcal{{O}}({stepSymbol}^{{{format_order(order)}}})$",
        xy=(stepRange[1], reference[1]),
        xytext=(4, 0),
        textcoords="offset points",
        va="center",
        fontsize=9,
        color=color,
    )


def plot_convergence(fileName, results, stepSymbol, stepUnit, step, title):
    """
    Plot all error norms over the step size together with the expected slopes

    results maps a label to a tuple (runs, expected orders)
    """
    plt = import_pyplot()
    if plt is None:
        return

    fig, axes = plt.subplots(2, 2, figsize=(10, 8), constrained_layout=True)
    colors = plt.rcParams["axes.prop_cycle"].by_key()["color"]
    for ax, key in zip(axes.flat, ERROR_KEYS):
        for (label, (runs, expected)), color in zip(results.items(), colors):
            steps = [step(r) for r in runs]
            errors = [r.errors[key] for r in runs]
            ax.loglog(steps, errors, "o-", color=color, label=label)
            if key in expected:
                add_reference_slope(ax, steps, errors, expected[key], stepSymbol, color)
        ax.set_xlabel(rf"${stepSymbol}$ [{stepUnit}]")
        ax.set_ylabel(ERROR_LABELS[key])
        ax.grid(True, which="both", alpha=0.3)
        ax.margins(x=0.15)

    # one legend for all panels, with a proxy entry for the reference slopes
    handles, labels = axes.flat[0].get_legend_handles_labels()
    handles.append(plt.Line2D([], [], color="gray", linestyle="--", linewidth=1.0))
    labels.append("expected order")
    axes.flat[0].legend(handles, labels, loc="lower right")
    fig.suptitle(title)
    fig.savefig(fileName, dpi=200)
    plt.close(fig)


def check_rates(label, runs, expected, step):
    """Check that the rate of the finest refinement matches the expected order"""
    success = True
    for key, order in expected.items():
        rate = rates(runs, key, step)[-1]
        if rate < order - ORDER_TOLERANCE:
            print(f"FAILED: {label}: {key} rate {rate:.2f} < expected {order} - {ORDER_TOLERANCE}")
            success = False
        else:
            print(f"passed: {label}: {key} rate {rate:.2f} (expected {order})")
    return success


def spatial_study(dim, schemes, levels, test):
    """Stationary problem under uniform grid refinement"""
    results = {}
    success = True
    tableFile = f"taylorgreen_spatial_{dim}d.md"
    Path(tableFile).unlink(missing_ok=True)
    for scheme in schemes:
        runs = []
        for level in range(levels):
            cells = BASE_CELLS[dim] * 2**level
            name = f"taylorgreen_spatial_{dim}d_{scheme}_{level}"
            rows = run_simulation(
                executable(dim, scheme),
                dim,
                name,
                ["-Problem.IsStationary", "true", "-Grid.Cells", " ".join([str(cells)] * dim)],
            )
            runs.append(to_run(name, rows))

        results[SCHEME_LABELS[scheme]] = (runs, EXPECTED_SPATIAL_ORDER[scheme])
        write_table(tableFile, f"{dim}D {scheme} (stationary)", runs, "h", lambda r: r.h)
        if test:
            expected = EXPECTED_SPATIAL_ORDER[scheme]
            success &= check_rates(f"{dim}D {scheme}", runs, expected, lambda r: r.h)

    if not test:
        plot_convergence(
            f"taylorgreen_spatial_{dim}d.png",
            results,
            "h",
            "m",
            lambda r: r.h,
            f"{dim}D Taylor-Green vortex, stationary, spatial convergence",
        )
    return success


def temporal_study(dim, schemes, levels, test):
    """Instationary problem on a fine grid under refinement of the time step size"""
    results = {}
    success = True
    tableFile = f"taylorgreen_temporal_{dim}d.md"
    Path(tableFile).unlink(missing_ok=True)
    # the spatial error has to be small compared to the temporal error
    cells = BASE_CELLS[dim] * (2 if test else 4)
    tEnd = 0.4
    for scheme in schemes:
        # higher-order methods require the multi-stage executable, otherwise use implicit Euler
        exe = executable(dim, scheme, multistage=True)
        timeSchemes = TIME_SCHEMES
        if not Path(exe).is_file():
            exe = executable(dim, scheme)
            timeSchemes = ("ImplicitEuler",)
        if not Path(exe).is_file():
            print(f"Skipping temporal study for {dim}D {scheme}: {exe} not found")
            continue

        for timeScheme in timeSchemes:
            runs = []
            for level in range(levels):
                dt = tEnd / 2 ** (level + 1)
                name = f"taylorgreen_temporal_{dim}d_{scheme}_{timeScheme}_{level}"
                rows = run_simulation(
                    exe,
                    dim,
                    name,
                    [
                        "-Problem.IsStationary", "false",
                        "-Grid.Cells", " ".join([str(cells)] * dim),
                        "-TimeLoop.Scheme", timeScheme,
                        "-TimeLoop.TEnd", str(tEnd),
                        "-TimeLoop.DtInitial", str(dt),
                        "-TimeLoop.MaxTimeStepSize", str(dt),
                    ],
                )
                runs.append(to_run(name, rows))

            label = f"{scheme} {timeScheme}"
            expected = {"velocityL2": EXPECTED_TEMPORAL_ORDER[timeScheme]}
            results[f"{SCHEME_LABELS[scheme]}, {TIME_SCHEME_LABELS[timeScheme]}"] = (runs, expected)
            write_table(tableFile, f"{dim}D {label} (t = {tEnd})", runs, "dt", lambda r: r.dt)
            if test:
                success &= check_rates(f"{dim}D {label}", runs, expected, lambda r: r.dt)

    if results and not test:
        plot_convergence(
            f"taylorgreen_temporal_{dim}d.png",
            results,
            r"\Delta t",
            "s",
            lambda r: r.dt,
            f"{dim}D Taylor-Green vortex, temporal convergence at $t = {tEnd}$ s",
        )
    return success


def energy_study(dim, schemes):
    """Kinetic energy decay of the instationary problem"""
    results = {}
    for scheme in schemes:
        name = f"taylorgreen_energy_{dim}d_{scheme}"
        rows = run_simulation(
            executable(dim, scheme), dim, name, ["-Problem.IsStationary", "false"]
        )
        results[scheme] = rows

        # dissipation rate -dE/dt by finite differences compared to the analytical rate
        last, previous = rows[-1], rows[-2]
        dissipation = -(last["kineticEnergy"] - previous["kineticEnergy"]) / last["dt"]
        exactDissipation = (
            -(last["kineticEnergyExact"] - previous["kineticEnergyExact"]) / last["dt"]
        )
        energyRatio = last["kineticEnergy"] / last["kineticEnergyExact"]
        print(
            f"{dim}D {scheme}: E(T)/E_exact(T) = {energyRatio:.6f}, "
            f"dissipation rate {dissipation:.6e} (exact {exactDissipation:.6e})"
        )

    plt = import_pyplot()
    if plt is None:
        return True

    fig, (axEnergy, axError) = plt.subplots(1, 2, figsize=(10, 4), constrained_layout=True)
    colors = plt.rcParams["axes.prop_cycle"].by_key()["color"]
    exact = next(iter(results.values()))
    energy0 = exact[0]["kineticEnergyExact"]
    times = [r["t"] for r in exact]
    axEnergy.plot(
        times,
        [r["kineticEnergyExact"] / energy0 for r in exact],
        "k--",
        label=r"exact, $F^2(t)$",
    )
    for (scheme, rows), color in zip(results.items(), colors):
        t = [r["t"] for r in rows]
        label = SCHEME_LABELS[scheme]
        axEnergy.plot(t, [r["kineticEnergy"] / energy0 for r in rows], color=color, label=label)
        relativeError = [
            (r["kineticEnergy"] - r["kineticEnergyExact"]) / r["kineticEnergyExact"] for r in rows
        ]
        axError.plot(t, relativeError, "o-", color=color, label=label)

    axEnergy.set_xlabel(r"$t$ [s]")
    axEnergy.set_ylabel(r"$E_h(t) / E(0)$")
    axEnergy.legend()
    axError.set_xlabel(r"$t$ [s]")
    axError.set_ylabel(r"$\left(E_h(t) - E(t)\right) / E(t)$")
    axError.legend()
    for ax in (axEnergy, axError):
        ax.grid(True, alpha=0.3)

    dt = exact[1]["dt"] if len(exact) > 1 else exact[0]["dt"]
    fig.suptitle(
        rf"{dim}D Taylor-Green vortex, kinetic energy "
        rf"$E(t) = \frac{{1}}{{2}}\int_\Omega \rho \|\mathbf{{u}}\|^2 \, \mathrm{{d}}x$, "
        rf"$\Delta t = {dt:g}$ s"
    )
    fig.savefig(f"taylorgreen_energy_{dim}d.png", dpi=200)
    plt.close(fig)
    return True


def check_solution():
    """Verify symbolically that the analytical solutions solve the Navier-Stokes equations"""
    import sympy as sp

    x, y, z, t, k, nu, rho, u0 = sp.symbols("x y z t k nu rho U_0", positive=True)

    def residuals(velocity, pressure, coords):
        dim = len(coords)
        div = sp.simplify(sum(sp.diff(velocity[i], coords[i]) for i in range(dim)))
        momentum = [
            sp.simplify(
                rho * sp.diff(velocity[i], t)
                + rho * sum(velocity[j] * sp.diff(velocity[i], coords[j]) for j in range(dim))
                - rho * nu * sum(sp.diff(velocity[i], c, 2) for c in coords)
                + sp.diff(pressure, coords[i])
            )
            for i in range(dim)
        ]
        return [div] + momentum

    decay2d = sp.exp(-2 * nu * k**2 * t)
    velocity2d = [
        u0 * sp.sin(k * x) * sp.cos(k * y) * decay2d,
        -u0 * sp.cos(k * x) * sp.sin(k * y) * decay2d,
    ]
    pressure2d = rho * u0**2 / 4 * (sp.cos(2 * k * x) + sp.cos(2 * k * y)) * decay2d**2

    decay3d = sp.exp(-3 * nu * k**2 * t)
    scale = 4 * sp.sqrt(2) / (3 * sp.sqrt(3)) * u0

    def component(a, b, c):
        return (
            scale
            * (
                sp.sin(k * a - 5 * sp.pi / 6) * sp.cos(k * b - sp.pi / 6) * sp.sin(k * c)
                - sp.cos(k * c - 5 * sp.pi / 6) * sp.sin(k * a - sp.pi / 6) * sp.sin(k * b)
            )
            * decay3d
        )

    velocity3d = [component(x, y, z), component(y, z, x), component(z, x, y)]
    pressure3d = -rho / 2 * sum(v**2 for v in velocity3d)

    success = True
    for label, velocity, pressure, coords in (
        ("2D", velocity2d, pressure2d, [x, y]),
        ("3D", velocity3d, pressure3d, [x, y, z]),
    ):
        result = residuals(velocity, pressure, coords)
        print(f"{label}: divergence = {result[0]}, momentum residual = {result[1:]}")
        success &= all(r == 0 for r in result)
    return success


def main():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("--dim", type=int, nargs="+", choices=(2, 3), default=[2, 3])
    parser.add_argument("--schemes", nargs="+", choices=SCHEMES, default=list(SCHEMES))
    parser.add_argument(
        "--study",
        nargs="+",
        choices=("spatial", "temporal", "energy", "all"),
        default=["all"],
    )
    parser.add_argument("--levels", type=int, default=None, help="number of refinement levels")
    parser.add_argument(
        "--test", action="store_true", help="short study checking the convergence rates"
    )
    parser.add_argument(
        "--check-solution",
        action="store_true",
        help="verify the analytical solutions (requires sympy)",
    )
    args = parser.parse_args()

    if args.check_solution:
        sys.exit(0 if check_solution() else 1)

    studies = {"spatial", "temporal", "energy"} if "all" in args.study else set(args.study)
    if args.levels is not None and args.levels < 2:
        parser.error("At least two refinement levels are needed to compute convergence rates")

    success = True
    for dim in args.dim:
        levels = args.levels if args.levels is not None else DEFAULT_LEVELS[dim][args.test]
        if "spatial" in studies:
            success &= spatial_study(dim, args.schemes, levels, args.test)
        if "temporal" in studies:
            success &= temporal_study(dim, args.schemes, levels, args.test)
        if "energy" in studies:
            success &= energy_study(dim, args.schemes)

    sys.exit(0 if success else 1)


if __name__ == "__main__":
    main()
