#!/usr/bin/env python3
# SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

"""
Run two commands and check that the errors they write to their CSV files agree.

Each executable has to write a CSV file with a header row and one data row (see
test/testing/convergencetest.hh). All columns have to agree up to the given relative tolerance.

Example:
  compare_errors.py --command-a "./test_a params.input -Problem.Name a"
                    --command-b "./test_b params.input -Problem.Name b"
                    --errors-a a-errors.csv --errors-b b-errors.csv --tolerance 1e-10
"""

import argparse
import csv
import shlex
import subprocess
import sys


def readErrors(fileName):
    """Read the header and the single data row of an error file"""
    with open(fileName, newline="") as file:
        rows = list(csv.reader(file))
    if len(rows) != 2:
        sys.exit(f"{fileName} must have a header row and one data row")
    return rows[0], [float(x) for x in rows[1]]


def main():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("--command-a", required=True)
    parser.add_argument("--command-b", required=True)
    parser.add_argument("--errors-a", required=True)
    parser.add_argument("--errors-b", required=True)
    parser.add_argument("--tolerance", type=float, default=1e-10)
    parser.add_argument("--columns", nargs="+", default=None,
                        help="Compare only these columns (default: all, which have to be the same in both files)")
    args = parser.parse_args()

    for command in (args.command_a, args.command_b):
        subprocess.run(shlex.split(command), check=True, stdout=subprocess.DEVNULL)

    headerA, valuesA = readErrors(args.errors_a)
    headerB, valuesB = readErrors(args.errors_b)
    if args.columns is None and headerA != headerB:
        sys.exit(f"Different columns: {headerA} and {headerB}")

    columnsA = dict(zip(headerA, valuesA))
    columnsB = dict(zip(headerB, valuesB))
    success = True
    for name in (args.columns if args.columns is not None else headerA):
        a, b = columnsA[name], columnsB[name]
        relDiff = abs(a - b) / max(abs(a), abs(b)) if a != b else 0.0
        print(f"{name}: {a:.16e} {b:.16e} relative difference {relDiff:.2e}")
        success = success and relDiff <= args.tolerance

    if not success:
        sys.exit(f"The errors differ by more than the relative tolerance {args.tolerance}")


if __name__ == "__main__":
    main()
