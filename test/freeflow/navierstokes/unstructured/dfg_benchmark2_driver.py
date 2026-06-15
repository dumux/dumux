#!/usr/bin/env python3
# SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later
"""
DFG flow-around-cylinder benchmark 2D-2 (Re = 100) driver.

Runs the unstructured cylinder Navier-Stokes test long enough to reach the oscillatory
steady state (Karman vortex street), then evaluates the benchmark quantities from the
recorded drag/lift/pressure-difference time series and compares them against the FeatFlow
reference values.

Reference (see benchmarkreferences.md):
  https://wwwold.mathematik.tu-dortmund.de/~featflow/en/benchmarks/cfdbenchmarking/flow/dfg_benchmark2_re100.html

This is intentionally NOT a ctest: it runs many time steps and takes minutes, but it makes
the benchmark result reproducible. By default it uses the most efficient solver available
here: the Trilinos Amesos2 parallel direct solver.

Usage (from the build directory containing the test binary + params.input + the .msh):
  ./dfg_benchmark2_driver.py                      # run + analyze (Trilinos, 2 ranks)
  ./dfg_benchmark2_driver.py --ranks 4 --tend 10  # more ranks, longer
  ./dfg_benchmark2_driver.py --solver simple      # use the matrix-free SIMPLE solver
  ./dfg_benchmark2_driver.py --analyze-only f.csv # re-analyze an existing time series
"""
import argparse
import csv
import os
import subprocess
import sys

# FeatFlow DFG 2D-2 (Re=100) reference intervals
REFERENCE = {
    "cD_max":    (3.2200, 3.2400),
    "cL_max":    (0.9900, 1.0100),
    "Strouhal":  (0.2950, 0.3050),
    "dP":        (2.4600, 2.5000),
}
# characteristic quantities of the canonical DFG 2D-2 benchmark
D_CYL = 0.1     # cylinder diameter
U_MEAN = 1.0    # mean inflow velocity (= 2/3 * U_max; U_max = 1.5).
                # Re = U_MEAN * D_CYL / nu = 100 for nu = 1e-3. At this scale the dimensional
                # pressure difference dP is directly comparable to the reference (it scales as U^2).


def solver_args(solver):
    """Return the linear-solver command-line arguments for the chosen solver."""
    if solver == "trilinos":
        return ["-LinearSolver.UseTrilinos", "true"]
    if solver == "simple":
        return [
            "-LinearSolver.UseSimpleSolver", "true",
            "-LinearSolver.Type", "fgmres",
            "-LinearSolver.MaxIterations", "500",
            "-LinearSolver.ResidualReduction", "1e-7",
            "-LinearSolver.Preconditioner.SchurReduction", "1e-2",
            "-LinearSolver.Preconditioner.SchurIterations", "30",
            "-LinearSolver.Preconditioner.PressureRelaxation", "1.0",
        ]
    if solver == "umfpack":
        return []  # sequential direct solver (main.cc default branch)
    raise ValueError(f"unknown solver '{solver}'")


def run_simulation(args):
    name = "dfg_benchmark2_re100"
    csv_file = f"{name}_benchmark_indicators.csv"
    cmd = []
    if args.solver != "umfpack" and args.ranks > 1:
        cmd += [args.mpiexec, "-np", str(args.ranks)]
    cmd += [
        os.path.abspath(args.exe), "params.input",
        "-Problem.Name", name,
        "-Problem.EnableInertiaTerms", "true",
        "-Problem.Instationary", "true",
        "-Problem.InflowMaxVelocity", str(args.inflow),
        "-Component.LiquidDynamicViscosity", str(args.viscosity),
        "-TimeLoop.DtInitial", str(args.dt),
        "-TimeLoop.MaxTimeStepSize", str(args.dt),
        "-TimeLoop.TEnd", str(args.tend),
        "-Newton.MaxSteps", "20",
    ]
    cmd += solver_args(args.solver)
    # DuMux assembly is multithreaded (element coloring); with several MPI ranks we cap the
    # threads per rank via DUMUX_NUM_THREADS so that ranks x threads does not oversubscribe the
    # cores (and, conversely, so the otherwise single-threaded assembly does use several cores).
    env = os.environ.copy()
    env["DUMUX_NUM_THREADS"] = str(args.threads)
    print(f"Running (DUMUX_NUM_THREADS={args.threads}):\n  " + " ".join(cmd) + "\n", flush=True)
    subprocess.run(cmd, check=True, env=env)
    return csv_file


def read_series(csv_file):
    t, cD, cL, dP = [], [], [], []
    with open(csv_file, newline="") as f:
        for row in csv.DictReader(f):
            t.append(float(row["time"]))
            cD.append(float(row["cDrag"]))
            cL.append(float(row["cLift"]))
            dP.append(float(row["pDiff"]))
    return t, cD, cL, dP


def local_maxima(t, y, min_value=None):
    """Indices of strict local maxima (optionally above a threshold value)."""
    idx = []
    for i in range(1, len(y) - 1):
        if y[i] > y[i - 1] and y[i] >= y[i + 1]:
            if min_value is None or y[i] > min_value:
                idx.append(i)
    return idx


def interpolate(t, y, t_query):
    """Linear interpolation of y(t) at t_query."""
    if t_query <= t[0]:
        return y[0]
    if t_query >= t[-1]:
        return y[-1]
    for i in range(1, len(t)):
        if t[i] >= t_query:
            w = (t_query - t[i - 1]) / (t[i] - t[i - 1])
            return y[i - 1] + w * (y[i] - y[i - 1])
    return y[-1]


def analyze(t, cD, cL, dP, skip_fraction):
    """Extract the benchmark quantities from the (transient-trimmed) time series."""
    n = len(t)
    if n < 5:
        raise RuntimeError("time series too short to analyze")

    # trim the initial transient
    t0_cut = t[0] + skip_fraction * (t[-1] - t[0])
    start = next((i for i, ti in enumerate(t) if ti >= t0_cut), 0)
    ts, cDs, cLs, dPs = t[start:], cD[start:], cL[start:], dP[start:]

    amplitude = 0.5 * (max(cLs) - min(cLs))

    # saturation check: compare the lift amplitude in the first vs the second half of the
    # analysis window. While the Karman street is still developing the amplitude grows
    # exponentially, so a much larger late amplitude means the limit cycle is NOT yet reached.
    half = len(cLs) // 2
    amp_early = (max(cLs[:half]) - min(cLs[:half])) if half > 1 else 0.0
    amp_late = (max(cLs[half:]) - min(cLs[half:])) if half > 1 else 0.0
    saturated = amp_late <= 1.15 * amp_early if amp_early > 0 else False

    # shedding frequency from the lift-coefficient maxima
    peaks = local_maxima(ts, cLs, min_value=0.0)
    period = strouhal = None
    if len(peaks) >= 2:
        periods = [ts[peaks[k + 1]] - ts[peaks[k]] for k in range(len(peaks) - 1)]
        period = sum(periods) / len(periods)
        if period > 0:
            strouhal = (1.0 / period) * D_CYL / U_MEAN

    # max drag/lift over the last full period (or the whole analysis window as fallback)
    if len(peaks) >= 2:
        lo, hi = peaks[-2], peaks[-1] + 1
    else:
        lo, hi = 0, len(ts)
    cD_max = max(cDs[lo:hi])
    cL_max = max(cLs[lo:hi])

    # pressure difference at t(cL_max) + T/2 (benchmark definition); fallback: at t(cL_max)
    i_clmax = lo + max(range(hi - lo), key=lambda k: cLs[lo + k])
    t_clmax = ts[i_clmax]
    t_dp = t_clmax + 0.5 * period if period else t_clmax
    dP_bench = interpolate(ts, dPs, t_dp)

    return {
        "n_steps": n,
        "analysis_window": (ts[0], ts[-1]),
        "n_peaks": len(peaks),
        "cL_amplitude": amplitude,
        "amp_early": amp_early,
        "amp_late": amp_late,
        "saturated": saturated,
        "period": period,
        "Strouhal": strouhal,
        "cD_max": cD_max,
        "cL_max": cL_max,
        "dP": dP_bench,
    }


def report(res, tol):
    print("\n================ DFG 2D-2 (Re=100) benchmark result ================")
    print(f"time steps analyzed : {res['n_steps']}  "
          f"(window t = {res['analysis_window'][0]:.3f} .. {res['analysis_window'][1]:.3f})")
    print(f"lift maxima found   : {res['n_peaks']}")
    print(f"lift amplitude      : {res['cL_amplitude']:.4f}")
    if res["period"]:
        print(f"shedding period T   : {res['period']:.4f} s  (f = {1.0/res['period']:.4f} Hz)")
    print()
    if res["n_peaks"] < 2 or res["cL_amplitude"] < 1e-3:
        print("WARNING: no developed oscillation detected -- run longer (increase --tend) "
              "and/or reduce --dt. The flow has not reached the oscillatory steady state yet.")
        return False

    if not res["saturated"]:
        print(f"WARNING: the oscillation is NOT saturated -- the lift amplitude is still growing "
              f"(early {res['amp_early']:.3f} -> late {res['amp_late']:.3f} over the analysis "
              f"window). The Karman street has not reached its limit cycle; increase --tend "
              f"(the frequency/Strouhal already locks in early, but cD_max/cL_max/dP will keep "
              f"rising until saturation). Reported amplitudes are lower bounds.\n")

    ok = True
    header = f"{'quantity':<10}{'computed':>14}{'reference range':>22}{'status':>10}"
    print(header)
    print("-" * len(header))
    for key in ("cD_max", "cL_max", "Strouhal", "dP"):
        val = res[key]
        lo, hi = REFERENCE[key]
        if val is None:
            print(f"{key:<10}{'n/a':>14}{f'[{lo}, {hi}]':>22}{'SKIP':>10}")
            continue
        # accept the reference interval widened by a relative tolerance (the 1st-order
        # implicit-Euler time scheme and this mesh resolution introduce some deviation)
        wlo, whi = lo * (1 - tol), hi * (1 + tol)
        passed = wlo <= val <= whi
        ok = ok and passed
        print(f"{key:<10}{val:>14.4f}{f'[{lo}, {hi}]':>22}{'PASS' if passed else 'FAIL':>10}")
    print("-" * len(header))
    print(f"(reference intervals widened by tol = {tol:.0%} for the comparison)")
    return ok


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--exe",
                   default="./test_ff_stokes_dfg_benchmark_stationary_pq1bubble_box_trilinos",
                   help="path to the test binary (default: the Trilinos build target)")
    p.add_argument("--solver", choices=["trilinos", "simple", "umfpack"], default="trilinos",
                   help="linear solver (default: trilinos parallel direct = most efficient)")
    p.add_argument("--ranks", type=int, default=2, help="number of MPI ranks (default: 2)")
    p.add_argument("--threads", type=int, default=2,
                   help="DUMUX_NUM_THREADS per rank for the multithreaded assembly (default: 2)")
    p.add_argument("--mpiexec", default=os.environ.get("MPIEXEC", "mpirun"))
    p.add_argument("--dt", type=float, default=0.01, help="time step size (default: 0.01)")
    p.add_argument("--tend", type=float, default=8.0, help="end time (default: 8.0)")
    p.add_argument("--inflow", type=float, default=1.5,
                   help="peak inflow velocity U_max (default 1.5 -> U_mean=1.0, canonical DFG 2D-2)")
    p.add_argument("--viscosity", type=float, default=1e-3,
                   help="dynamic viscosity (default 1e-3 -> Re=100 at U_mean=1.0)")
    p.add_argument("--skip-fraction", type=float, default=0.5,
                   help="fraction of the time series discarded as transient (default: 0.5)")
    p.add_argument("--tol", type=float, default=0.05,
                   help="relative tolerance widening the reference intervals (default: 0.05)")
    p.add_argument("--analyze-only", metavar="CSV", default=None,
                   help="skip the run and analyze an existing time-series CSV")
    args = p.parse_args()

    csv_file = args.analyze_only if args.analyze_only else run_simulation(args)
    t, cD, cL, dP = read_series(csv_file)
    res = analyze(t, cD, cL, dP, args.skip_fraction)
    ok = report(res, args.tol)
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
