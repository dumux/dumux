#!/usr/bin/env python3
# SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

"""
Run two commands on a sequence of refinement levels and check that the mean convergence
rates of the errors of the first command are not lower than those of the second by more
than a tolerance.

Each executable has to write a CSV file with a header row and one data row, with the mesh size
in the first column (see test/testing/convergencetest.hh). The rate of an error between two
levels is log(e_i/e_{i+1})/log(h_i/h_{i+1}).

Example:
  compare_rates.py --command-a "./test_a params.input -Problem.Name a"
                   --command-b "./test_b params.input -Problem.Name b"
                   --errors-a a-errors.csv --errors-b b-errors.csv
                   --levels "6 12" "12 24" --level-args "-Vessel.Grid.Cells {0} -Tissue.Grid.Cells '{1} {1} {0}'"
                   --tolerance 0.05
"""

import argparse
import csv
import math
import shlex
import subprocess
import sys


def readErrors(fileName):
    """Read the header and the single data row of an error file"""
    with open(fileName, newline="") as file:
        rows = list(csv.reader(file))
    if len(rows) != 2:
        sys.exit(f"{fileName} must have a header row and one data row")
    return dict(zip(rows[0], [float(x) for x in rows[1]]))


def runLevels(command, errorFile, levels, levelArgs):
    """Run the command for each level and return the errors of each level"""
    result = []
    for level in levels:
        call = shlex.split(command) + shlex.split(levelArgs.format(*level.split()))
        subprocess.run(call, check=True, stdout=subprocess.DEVNULL)
        result.append(readErrors(errorFile))
    return result


def meanRates(rows):
    """The mean over the levels of the rates of each error"""
    rates = {}
    for name in rows[0]:
        if name == "h":
            continue
        values = [
            math.log(r0[name] / r1[name]) / math.log(r0["h"] / r1["h"])
            for r0, r1 in zip(rows, rows[1:])
        ]
        rates[name] = sum(values) / len(values)
    return rates


def main():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("--command-a", required=True)
    parser.add_argument("--command-b", required=True)
    parser.add_argument("--errors-a", required=True)
    parser.add_argument("--errors-b", required=True)
    parser.add_argument("--levels", nargs="+", required=True)
    parser.add_argument("--level-args", required=True)
    parser.add_argument("--tolerance", type=float, default=0.0)
    parser.add_argument("--max-error-ratio", type=float, default=None,
                        help="Also check that the error of the first command on the finest level is at most this multiple of the second")
    parser.add_argument("--min-rate", type=float, default=None,
                        help="Also check that the mean rates of both commands are at least this value")
    args = parser.parse_args()

    rowsA = runLevels(args.command_a, args.errors_a, args.levels, args.level_args)
    rowsB = runLevels(args.command_b, args.errors_b, args.levels, args.level_args)
    ratesA, ratesB = meanRates(rowsA), meanRates(rowsB)

    success = True
    for name, rateA in ratesA.items():
        ok = rateA >= ratesB[name] - args.tolerance
        print(f"{name}: mean rate {rateA:.3f} (a), {ratesB[name]:.3f} (b), {'ok' if ok else 'a is lower than b'}")
        success = success and ok
        if args.min_rate is not None:
            rateOk = min(rateA, ratesB[name]) >= args.min_rate
            print(f"{name}: mean rates at least {args.min_rate}: {'ok' if rateOk else 'too low'}")
            success = success and rateOk
        if args.max_error_ratio is not None:
            ratio = rowsA[-1][name] / rowsB[-1][name]
            ratioOk = ratio <= args.max_error_ratio
            print(f"{name}: error on the finest level {ratio:.2f} times that of b, {'ok' if ratioOk else 'too large'}")
            success = success and ratioOk

    if not success:
        sys.exit("A check of the rates or errors failed")


if __name__ == "__main__":
    main()
