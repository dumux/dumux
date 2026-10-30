#!/usr/bin/env python3
# SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

"""
Runs the 1pni conduction test with refined cells (space study) or refined time steps (time study)
and checks the convergence rate of the temperature error at the end time with respect to the
analytical solution. The executable appends the errors to <name>-errors.csv.

Space study: the time step size is fixed and scaled with the squared cell size,
             so the time discretization error stays a small, constant fraction of the error.
Time study:  the grid is fine enough for the spatial error to be negligible,
             the time step size is fixed and halved.

usage: convergencetest.py <executable> [input file] --study space|time [--full] [runtime arguments]

With --full, more refinements are run (used for the benchmark documentation).
"""

import os
import sys
import subprocess
from math import log

LENGTH = 5.0  # [m] domain length in x-direction
END_TIME = 1e5  # [s]

# space study: dt = SPACE_DT_COARSEST*(SPACE_CELLS[0]/cells)^2, all dt divide the end time
SPACE_CELLS = {"test": [25, 50, 100, 200], "full": [25, 50, 100, 200, 400]}
SPACE_DT_COARSEST = 2000.0

# time study: fine grid, dt halved
TIME_CELLS = 3200
TIME_DTS = {"test": [2e4, 1e4, 5e3, 2.5e3], "full": [2e4, 1e4, 5e3, 2.5e3, 1250.0, 625.0]}

EXPECTED_RATE = {"space": 2.0, "time": 1.0}
TOLERANCE = 0.2


def study_setup(study, full=False):
    """Return a list of (cells, time step size) for the given study."""
    level = "full" if full else "test"
    if study == "space":
        cells = SPACE_CELLS[level]
        return [(n, SPACE_DT_COARSEST*(cells[0]/n)**2) for n in cells]
    if study == "time":
        return [(TIME_CELLS, dt) for dt in TIME_DTS[level]]
    raise ValueError(f"Unknown study '{study}', use 'space' or 'time'")


def run_study(executable, study, input_file="params.input", full=False, name=None, extra_args=()):
    """Run the executable for each (cells, dt) and return the cell sizes, time step sizes and errors."""
    exe = executable if os.path.isabs(executable) else "./" + executable
    name = name or f"{os.path.basename(executable)}_{study}"
    error_file = name + "-errors.csv"
    if os.path.exists(error_file):
        os.remove(error_file)

    for cells, dt in study_setup(study, full):
        subprocess.run(
            [exe, input_file, "-Problem.Name", name,
             "-Grid.Cells", f"{cells} 1",
             "-TimeLoop.TEnd", str(END_TIME),
             "-TimeLoop.DtInitial", str(dt), "-TimeLoop.MaxTimeStepSize", str(dt),
             "-Problem.OutputInterval", "1000000000",
             *extra_args],
            check=True,
            stdout=subprocess.DEVNULL,
        )

    with open(error_file) as f:
        rows = [[float(v) for v in line.split(",")] for line in f if line.strip()]
    return [row[0] for row in rows], [row[1] for row in rows], [row[2] for row in rows]


def compute_rates(sizes, errors):
    """Convergence rates between successive refinements."""
    return [log(errors[i]/errors[i + 1])/log(sizes[i]/sizes[i + 1]) for i in range(len(errors) - 1)]


def print_table(study, hs, dts, errors, rates):
    print(f"  {'cell size [m]':>14}  {'time step [s]':>14}  {'error [K]':>12}  {'rate':>6}")
    for i, (h, dt, error) in enumerate(zip(hs, dts, errors)):
        rate = f"{rates[i - 1]:.2f}" if i > 0 else "-"
        print(f"  {h:>14.6g}  {dt:>14.6g}  {error:>12.4e}  {rate:>6}")


if __name__ == "__main__":
    args = sys.argv[1:]
    if not args or "--study" not in args:
        sys.stderr.write(__doc__)
        sys.exit(1)

    i = args.index("--study")
    study = args[i + 1]
    del args[i:i + 2]
    full = "--full" in args
    if full:
        args.remove("--full")

    executable, extra_args = args[0], args[1:]
    input_file = "params.input"
    if extra_args and not extra_args[0].startswith("-"):
        input_file, extra_args = extra_args[0], extra_args[1:]

    hs, dts, errors = run_study(executable, study, input_file, full, extra_args=extra_args)
    sizes = hs if study == "space" else dts
    rates = compute_rates(sizes, errors)

    print(f"\n{study.capitalize()} convergence of the temperature error at t = {END_TIME:g} s for {executable}:")
    print_table(study, hs, dts, errors, rates)

    mean_rate = sum(rates)/len(rates)
    expected = EXPECTED_RATE[study]
    print(f"Mean convergence rate: {mean_rate:.2f} (expected {expected:g})")
    if not expected - TOLERANCE <= mean_rate <= expected + TOLERANCE:
        sys.stderr.write(f"Convergence rate {mean_rate:.2f} not close enough to {expected:g}!\n")
        sys.exit(1)
