#!/usr/bin/env python3
# SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later
"""Convergence of the VE model to the similarity solution of Nordbotten and Celia (2006).

The similarity solution neglects buoyancy (gravity number Gamma -> 0) and assumes a sharp interface
(capillary fringe -> 0). The script refines the radial cell size and the maximum time step size together
at Gamma = 1.4e-4 and an entry pressure of 1 Pa, for which the gas plume distance error is dominated by the discretization.

Usage (run from the build directory of the test, which contains params.input):
  python3 convergencetest.py <executable>
"""

import glob
import os
import subprocess
import sys
from math import log

# the paper case of params.input: 120 m^3/day for 10000 days
BASE_RATE = 1.388888889e-3
BASE_END_TIME = 8.64e8
BASE_MAX_TIME_STEP = 8.64e5
BASE_INITIAL_TIME_STEP = 1e3
RADIAL_EXTENT = 2000.0 - 0.1

GRID_STUDY_RATE_FACTOR = 1000
GRID_STUDY_ENTRY_PRESSURE = 1.0
GRID_STUDY_CELLS = [400, 800, 1600, 3200]
VERTICAL_CELLS = 30
EXPECTED_ORDER = 1.0


def run_case(executable, name, rate_factor=1, radial_cells=400, entry_pressure=None, initial_time_step=None, time_step_factor=1.0):
    """Run one case and return the gravity number and the final relative gas plume distance error."""
    exe = executable if os.path.isabs(executable) else "./" + executable
    command = [
        exe, "params.input",
        "-Problem.Name", name,
        "-Problem.InjectionRate", str(BASE_RATE*rate_factor),
        "-TimeLoop.TEnd", str(BASE_END_TIME/rate_factor),
        "-TimeLoop.MaxTimeStepSize", str(BASE_MAX_TIME_STEP/rate_factor*time_step_factor),
        "-TimeLoop.DtInitial", str(initial_time_step if initial_time_step else BASE_INITIAL_TIME_STEP/rate_factor),
        "-Grid.Cells", f"{radial_cells} {VERTICAL_CELLS}",
        "-Benchmark.MaxRelativeGasPlumeDistanceError", "1",
    ]
    if entry_pressure is not None:
        command += ["-SpatialParams.BrooksCoreyPcEntry", str(entry_pressure)]

    print(f"Starting simulation {name}: "
          f"grid={radial_cells} {VERTICAL_CELLS} (radial, vertical cells), "
          f"radial cell size={RADIAL_EXTENT/radial_cells:.6g} m, "
          f"maximum time step={BASE_MAX_TIME_STEP/rate_factor*time_step_factor:.6g} s, "
          f"injection rate={BASE_RATE*rate_factor:.6g} m^3/s",
          flush=True)
    output = subprocess.check_output(command, text=True)
    gravity_number, error = None, None
    for line in output.split("\n"):
        if line.startswith("Mobility ratio: "):
            gravity_number = float(line.split()[-1])
        if line.startswith("Relative gas plume distance error"):
            error = float(line.split()[-1])
    print(f"Finished simulation {name}: gravity number={gravity_number}, gas plume distance error={error}", flush=True)
    return gravity_number, error


def remove_outputs(name):
    for pattern in (f"{name}-*.vtu", f"{name}.pvd", f"fine_{name}-*.vtu", f"fine_{name}.pvd"):
        for f in glob.glob(pattern):
            os.remove(f)


def grid_study(executable):
    """Return the radial cell sizes and the gas plume distance errors for a small gravity number and a thin capillary fringe.

    The maximum time step size is refined together with the radial cell size.
    """
    cell_sizes, errors = [], []
    for cells in GRID_STUDY_CELLS:
        name = f"gridstudy_{cells}"
        _, error = run_case(executable, name, rate_factor=GRID_STUDY_RATE_FACTOR, radial_cells=cells,
                            entry_pressure=GRID_STUDY_ENTRY_PRESSURE, initial_time_step=0.01,
                            time_step_factor=GRID_STUDY_CELLS[0]/cells)
        remove_outputs(name)
        cell_sizes.append(RADIAL_EXTENT/cells)
        errors.append(error)
    return cell_sizes, errors


def compute_rates(sizes, errors):
    """Compute convergence rates from successive refinements."""
    return [log(errors[i]/errors[i + 1])/log(sizes[i]/sizes[i + 1]) for i in range(len(errors) - 1)]


def print_table(label, sizes, errors, rates):
    print(f"  {label:>12}  {'error':>12}  {'rate':>6}")
    for i, (h, e) in enumerate(zip(sizes, errors)):
        rate_str = f"{rates[i - 1]:.2f}" if i > 0 else "   -"
        print(f"  {h:>12.4g}  {e:>12.4e}  {rate_str:>6}")


if __name__ == "__main__":
    if len(sys.argv) < 2:
        sys.stderr.write("Please provide the test executable name as argument\n")
        sys.exit(1)

    testname = str(sys.argv[1])

    cell_sizes, errors = grid_study(testname)
    rates = compute_rates(cell_sizes, errors)
    print("\nGas plume distance error for decreasing radial cell size and time step size (Gamma = 1.4e-4, entry pressure 1 Pa):")
    print_table("dr [m]", cell_sizes, errors, rates)

    mean_rate = sum(rates)/len(rates)
    print(f"Mean convergence rate: {mean_rate:.2f}")

    tol = 0.2
    if not (EXPECTED_ORDER - tol <= mean_rate <= EXPECTED_ORDER + tol):
        sys.stderr.write("*"*70 + "\n")
        sys.stderr.write(f"Convergence rate {mean_rate:.2f} not close enough to {EXPECTED_ORDER}! Test failed.\n")
        sys.stderr.write("*"*70 + "\n")
        sys.exit(1)

    sys.exit(0)
