#!/usr/bin/env python3
# SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later
"""The sheared Kapton membrane of Wong and Pellegrino at the experimental shear displacements.

Runs the relaxation at the given shear displacements and reports, in millimetres, the crest
count along the mid-line against the counts of the experiment, the crest spacing along the
mid-line and the amplitude against the averages the experiment reports as functions of the
displacement, the wavelength across the crests against the bending-tension balance, the
wrinkle orientation, the shear resultant against the tension-field value and the strain
energy. The exponent of the wavelength in the load is fitted across the displacements.

The wrinkled sheet has many metastable states, one per crest count, so the state reached
depends on where the relaxation starts. By default every displacement starts from its own
tension-field seed at the wavelength of the balance; `--seed-crests` seeds a given number of
crests along the mid-line instead, one run per count, so that the states can be compared by
their energy; `--path` runs the displacements in sequence from the first seed, the way the
experiment loads the sheet.

Usage:
  python3 wrinkling_benchmark.py <executable> [--displacements 0.5 1 2 3] [--cells 144 48]
                                 [--seed-crests 13 15 17] [--path]
"""

import argparse
import os
import subprocess
import sys
from math import log, exp, sqrt

HERE = os.path.dirname(os.path.abspath(__file__))

WIDTH_MM = 128.0
LENGTH_MM = 380.0
THICKNESS_MM = 0.025
SCALE = WIDTH_MM/THICKNESS_MM**2

# Wong and Pellegrino 2006, Part I: crests counted along y = 64 mm (Fig. 4) while loading,
# and the averages over the central region (Fig. 8, t = 0.025 mm, read off the figure to
# about 0.3 mm and 0.005 mm): the crest spacing along x, which the paper calls the
# wavelength, and the amplitude
EXPERIMENT_CRESTS = {0.5: 15, 1.0: 17, 2.0: 18, 3.0: 19}
EXPERIMENT_SPACING = [(0.25, 32.3), (0.35, 29.2), (0.45, 27.0), (0.75, 25.0),
                      (1.05, 24.0), (1.5, 22.0), (2.25, 21.1), (2.85, 20.5)]
EXPERIMENT_AMPLITUDE = [(0.25, 0.150), (0.42, 0.200), (0.75, 0.250), (1.35, 0.297),
                        (1.65, 0.330), (2.25, 0.365), (2.85, 0.382)]


def interpolate(table, x):
    """Log-log interpolation within the measured range, None outside it"""
    if x < table[0][0] or x > table[-1][0]:
        return None
    for (x0, y0), (x1, y1) in zip(table, table[1:]):
        if x0 <= x <= x1:
            s = (log(x) - log(x0))/(log(x1) - log(x0))
            return exp(log(y0) + s*(log(y1) - log(y0)))
    return None


def fit_exponent(xs, ys):
    lx, ly = [log(x) for x in xs], [log(y) for y in ys]
    mx, my = sum(lx)/len(lx), sum(ly)/len(ly)
    return sum((a - mx)*(b - my) for a, b in zip(lx, ly))/sum((a - mx)**2 for a in lx)


def read_summary(path):
    with open(path) as f:
        lines = f.read().splitlines()
    names = lines[0].lstrip("#").split()
    return [dict(zip(names, (float(v) for v in line.split())))
            for line in lines[1:] if line.strip() and not line.startswith("#")]


def run(executable, outdir, cells, loads, extra):
    os.makedirs(outdir, exist_ok=True)
    command = [
        os.path.abspath(executable), os.path.join(HERE, "params.input"),
        "-Grid.Cells", f"{cells[0]} {cells[1]}",
        "-Problem.LoadSteps", " ".join(f"{load:.8g}" for load in loads),
        "-Problem.Name", "benchmark",
    ] + extra
    with open(os.path.join(outdir, "run.log"), "w") as log_file:
        log_file.write(" ".join(command) + "\n")
        log_file.flush()
        subprocess.call(command, cwd=outdir, stdout=log_file, stderr=subprocess.STDOUT)
    return read_summary(os.path.join(outdir, "benchmark_wrinkles.dat"))


def fmt(value, width, digits):
    return f"{value:>{width}.{digits}f}" if value is not None else f"{'-':>{width}}"


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("executable")
    parser.add_argument("--displacements", nargs="+", type=float, default=[0.5, 1.0, 2.0, 3.0],
                        help="shear displacements in mm")
    parser.add_argument("--cells", nargs=2, type=int, default=[144, 48])
    parser.add_argument("--seed-crests", nargs="+", type=int, default=[],
                        help="seed this many crests along the mid-line instead of the balance")
    parser.add_argument("--path", action="store_true",
                        help="run the displacements in sequence from the first seed")
    parser.add_argument("--prestrain", type=float, default=0.0,
                        help="strain across the sheet locked in before the shear, physical")
    parser.add_argument("--outdir", default="wrinkling_benchmark")
    parser.add_argument("--extra", nargs=argparse.REMAINDER, default=[],
                        help="further arguments passed to the executable")
    args = parser.parse_args()

    seeds = [(None, [])]
    if args.seed_crests:
        # the seed wave runs at 45 degrees, so n crests along the mid-line of length L are a
        # wavelength of L/(n sqrt 2) across the crests, in the scaled units of the executable
        seeds = [(n, ["-Problem.SeedWavelength", f"{LENGTH_MM/WIDTH_MM/(n*sqrt(2.0)):.8g}"])
                 for n in args.seed_crests]

    if args.prestrain != 0.0:
        # the executable takes strains in the scaled units, (W/t)^2 times the physical
        args.extra += ["-Problem.PreStrain", f"{args.prestrain*(WIDTH_MM/THICKNESS_MM)**2:.8g}"]

    rows = []
    for count, seed_args in seeds:
        tag = f"_n{count}" if count is not None else ""
        if args.path:
            outdir = os.path.join(os.path.abspath(args.outdir), f"path{tag}")
            loads = [d*SCALE for d in args.displacements]
            for row in run(args.executable, outdir, args.cells, loads, seed_args + args.extra):
                rows.append((count, row))
        else:
            for displacement in args.displacements:
                outdir = os.path.join(os.path.abspath(args.outdir), f"delta_{displacement:g}mm{tag}")
                for row in run(args.executable, outdir, args.cells, [displacement*SCALE], seed_args + args.extra):
                    rows.append((count, row))

    print(f"{'delta [mm]':>10} {'seed':>5} {'crests':>6} {'exp.':>5} {'spacing':>8} {'exp.':>6}"
          f" {'A [mm]':>7} {'exp.':>6} {'lambda':>7} {'balance':>8} {'normal':>7} {'Nxy/TF':>7}"
          f" {'energy':>13} {'steps':>6} {'steady':>6}")
    for count, row in rows:
        delta = row["load"]/SCALE
        print(f"{delta:>10.3g} {str(count) if count is not None else '-':>5}"
              f" {int(row['numCrests']):>6d} {str(EXPERIMENT_CRESTS.get(round(delta, 3), '-')):>5}"
              f" {row['crestSpacing']*WIDTH_MM:>8.1f} {fmt(interpolate(EXPERIMENT_SPACING, delta), 6, 1)}"
              f" {row['amplitude']*THICKNESS_MM:>7.3f} {fmt(interpolate(EXPERIMENT_AMPLITUDE, delta), 6, 3)}"
              f" {row['wavelength']*WIDTH_MM:>7.1f} {row['predictedWavelength']*WIDTH_MM:>8.1f}"
              f" {row['normalAngle']:>7.1f} {row['centralNxy']/row['tensionFieldNxy']:>7.3f}"
              f" {row['bendingEnergy'] + row['membraneEnergy']:>13.6e}"
              f" {int(row['relaxationSteps']):>6d} {int(row['steadyState']):>6d}")
    print("spacing and amplitude in mm along the mid-line, experiment read off Fig. 8 of Wong and Pellegrino;"
          " lambda across the crests against the bending-tension balance")
    by_seed = {}
    for count, row in rows:
        by_seed.setdefault(count, []).append(row)
    for count, group in by_seed.items():
        if len(group) > 1:
            print(f"seed {count if count is not None else 'balance'}: fitted exponent"
                  f" d log(lambda)/d log(delta) = "
                  f"{fit_exponent([r['load'] for r in group], [r['wavelength'] for r in group]):+.3f}"
                  "  (balance -0.250)")


if __name__ == "__main__":
    sys.exit(main())
