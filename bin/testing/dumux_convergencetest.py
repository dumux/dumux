#!/usr/bin/env python3
# SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

"""
Run an executable on a sequence of refinement levels and check the convergence rates of its errors.

After each run, the executable has to write a CSV file with a header row and one data row.
The first column is the mesh size h, the other columns are named errors. The rate of an
error e between two levels is log(e_i/e_{i+1})/log(h_i/h_{i+1}). A check passes without a
rate if its error is below the round-off value on all levels, for a solution that the
discretization represents exactly.

Example:
  dumux_convergencetest.py --command "./test_name params.input"
                           --levels 0 1 2 --level-args "-Grid.Refinement {}"
                           --errors test_name-errors.csv
                           --check "L2(p) mean 1.9" --check "L2(v) finest 1.9 2.1"
"""

import argparse
import csv
import math
import os
import shlex
import subprocess
import sys

STATISTICS = ("mean", "finest", "all")


def parseCheck(check):
    """Parse a check of the form '<error> <mean|finest|all> <min> [<max>]'"""
    tokens = shlex.split(check)
    if len(tokens) not in (3, 4) or tokens[1] not in STATISTICS:
        raise argparse.ArgumentTypeError(
            f"invalid check '{check}', expected '<error> <{'|'.join(STATISTICS)}> <min> [<max>]'"
        )
    return {
        "error": tokens[0],
        "statistic": tokens[1],
        "min": float(tokens[2]),
        "max": float(tokens[3]) if len(tokens) == 4 else math.inf,
    }


def readCmdParameters():
    """Read the command line parameters"""
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument(
        "-c",
        "--command",
        required=True,
        help="The executable and the arguments shared by all levels as a single string",
    )
    parser.add_argument(
        "-l",
        "--levels",
        nargs="+",
        required=True,
        help=(
            "One entry per refinement level, each entry one or more space separated values"
            " inserted into the level arguments"
        ),
    )
    parser.add_argument(
        "-a",
        "--level-args",
        required=True,
        help=(
            "The arguments of one level as a single string, with the values of the level"
            " inserted at {} or {0}, {1}, ..."
        ),
    )
    parser.add_argument(
        "-e",
        "--errors",
        required=True,
        help=(
            "The CSV file the executable writes: header row, then h and the errors of the run."
            " The values of the level are inserted as in the level arguments."
        ),
    )
    parser.add_argument(
        "-k",
        "--check",
        type=parseCheck,
        action="append",
        default=[],
        help=(
            "'<error> <mean|finest|all> <min> [<max>]': the mean rate, the rate between the two"
            " finest levels, or every rate of the error has to lie in [min, max]"
        ),
    )
    parser.add_argument(
        "-r",
        "--round-off",
        type=float,
        default=0.0,
        help=(
            "A check passes without a rate if its error is below this value on all levels"
            " (default: 0)"
        ),
    )
    return parser.parse_args()


def runLevels(command, levels, levelArgs, errors):
    """Run the command for each level and return the values of the error file per level"""
    rows = []
    for level in levels:
        try:
            call = shlex.split(command) + shlex.split(levelArgs.format(*level.split()))
            errorFile = errors.format(*level.split())
        except (IndexError, KeyError) as error:
            sys.exit(
                "A placeholder of the level arguments or the error file"
                f" has no value in level '{level}': {error}"
            )
        if os.path.exists(errorFile):
            os.remove(errorFile)
        print("Running " + shlex.join(call), flush=True)
        try:
            result = subprocess.run(call, check=False)
        except OSError as error:
            sys.exit(f"Could not run {call[0]}: {error}")
        if result.returncode != 0:
            sys.stderr.write(f"Run for level '{level}' failed with exit code {result.returncode}\n")
            sys.exit(result.returncode)
        try:
            with open(errorFile, newline="", encoding="utf-8") as csvFile:
                reader = csv.reader(csvFile, skipinitialspace=True)
                header, values = next(reader), next(reader)
        except (OSError, StopIteration):
            sys.exit(
                f"Run for level '{level}' did not write the error file {errorFile}"
                " (header row and data row)"
            )
        if header[0] != "h":
            sys.exit(f"The first column of {errorFile} has to be h, not {header[0]}")
        rows.append(dict(zip(header, map(float, values))))
    return rows


def convergenceRates(h, errors):
    """The rates between consecutive levels, NaN for a non-positive error or an unchanged h"""
    rates = []
    for e0, e1, h0, h1 in zip(errors, errors[1:], h, h[1:]):
        defined = e0 > 0.0 and e1 > 0.0 and h0 > 0.0 and h1 > 0.0 and h0 != h1
        rates.append(math.log(e0 / e1) / math.log(h0 / h1) if defined else math.nan)
    return rates


def printTable(levels, h, errors, rates):
    """Print the errors and rates of all levels"""
    print(f"\n{'level':>12}  {'h':>10}" + "".join(f"  {name:>14}  {'rate':>5}" for name in errors))
    for i, level in enumerate(levels):
        line = f"{level:>12}  {h[i]:>10.4e}"
        for name, values in errors.items():
            rate = f"{rates[name][i - 1]:.2f}" if i > 0 else "-"
            line += f"  {values[i]:>14.4e}  {rate:>5}"
        print(line)
    print()


def checkRates(checks, errors, rates, roundOff):
    """Apply each check to the rates, print the outcome and return whether all passed"""
    passed = True
    for check in checks:
        name, statistic = check["error"], check["statistic"]
        if name not in rates:
            sys.stderr.write(f"Error '{name}' is not in the error file ({', '.join(rates)})\n")
            passed = False
            continue
        if all(e < roundOff for e in errors[name]):
            print(f"{name}: below the round-off value {roundOff} on all levels")
            continue
        if statistic == "mean":
            values = [sum(rates[name]) / len(rates[name])]
        elif statistic == "finest":
            values = [rates[name][-1]]
        else:
            values = rates[name]
        ok = all(check["min"] <= v <= check["max"] for v in values)
        bounds = f"[{check['min']}, {check['max']}]"
        print(
            f"{name}: {statistic} rate {', '.join(f'{v:.3f}' for v in values)}"
            f" {'in' if ok else 'NOT in'} {bounds}"
        )
        passed = passed and ok
    return passed


def runConvergenceTest():
    """Run all levels, print the convergence table and check the rates"""
    args = readCmdParameters()
    if len(args.levels) < 2:
        sys.exit("At least two levels are needed for a convergence rate")

    rows = runLevels(args.command, args.levels, args.level_args, args.errors)
    for level, row in zip(args.levels, rows):
        nonFinite = [
            n
            for n in ["h"] + [c["error"] for c in args.check]
            if not math.isfinite(row.get(n, 0.0))
        ]
        if nonFinite:
            sys.exit(f"Non-finite value of {', '.join(nonFinite)} for level '{level}': {row}")
    h = [row["h"] for row in rows]
    errors = {name: [row[name] for row in rows] for name in rows[0] if name != "h"}
    rates = {name: convergenceRates(h, values) for name, values in errors.items()}

    printTable(args.levels, h, errors, rates)
    if not checkRates(args.check, errors, rates, args.round_off):
        sys.stderr.write("At least one convergence rate is out of its bounds\n")
        sys.exit(1)
    sys.exit(0)


if __name__ == "__main__":
    runConvergenceTest()
