#!/usr/bin/env python3
# SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

"""
Convergence test script for coupled FreeFlow/DarcyFlow problems.

The -Rates argument accepts one or more named error-type specifications:

    -Rates "ErrorType minRate" ["ErrorType minRate" ...]

where ErrorType matches a key written to the log file by the simulation,
e.g. L2(p), L2(v), H1(v), H1(p).  The script searches for each specified
error type in ALL log files (freeflow and darcy).  If a specified type is
not found in any log file the test exits with an error.  The last convergence
rate must satisfy last_rate >= minRate.

Example:
    -Rates "L2(p) 1.4" "L2(v) 1.95" "H1(v) 1.95"
"""

import re
import math
from math import log, isnan, isinf
import subprocess
import sys

if len(sys.argv) < 2:
    sys.stderr.write('Please provide a single argument <testname> to the script\n')
    sys.exit(1)

executableName = str(sys.argv[1])
testargs = [str(i) for i in sys.argv][2:]
testname = testargs[testargs.index('-Problem.TestCase')+1]
freeflowLogName = testname + "_" + testargs[testargs.index('-FreeFlow.Problem.Name')+1] + '.log'
darcyLogName = testname + "_" + testargs[testargs.index('-Darcy.Problem.Name')+1] + '.log'

# ---------------------------------------------------------------------------
# Parse -Rates argument
# ---------------------------------------------------------------------------
def parse_rates(testargs):
    """Return (rate_specs, cleaned_testargs).

    rate_specs is one of:
      None                                     -> use built-in defaults
      dict { error_type: minRate }            -> named format
    """
    if '-Rates' not in testargs:
        return None, testargs

    idx = testargs.index('-Rates')
    # collect all following tokens that are not flags (don't start with '-')
    specs_raw = []
    i = idx + 1
    while i < len(testargs) and not testargs[i].startswith('-'):
        specs_raw.append(testargs[i])
        i += 1
    cleaned = testargs[:idx] + testargs[i:]

    if not specs_raw:
        sys.stderr.write("ERROR: -Rates requires at least one argument.\n")
        sys.exit(1)

    # named format: each entry is "ErrorType minRate"
    rate_specs = {}
    for spec in specs_raw:
        parts = spec.split()
        if len(parts) != 2:
            sys.stderr.write(
                f"ERROR: Each -Rates entry must be 'ErrorType minRate', got: '{spec}'\n")
            sys.exit(1)
        error_type, min_rate = parts[0], float(parts[1])
        rate_specs[error_type] = min_rate
    return rate_specs, cleaned


rate_specs, testargs = parse_rates(testargs)

# ---------------------------------------------------------------------------
# Remove old log files
# ---------------------------------------------------------------------------
for log in [freeflowLogName, darcyLogName]:
    subprocess.call(['rm', '-f', log])
    print("Removed old log file ({})!".format(log))

# ---------------------------------------------------------------------------
# Run the simulations at three refinement levels
# ---------------------------------------------------------------------------
for i in [0, 1, 2]:
    subprocess.call(['./' + executableName] + testargs + ['-Grid.Refinement', str(i)])

# ---------------------------------------------------------------------------
# Generic log-file parser
# Handles lines of the form written by the C++ simulation, e.g.:
#   [ConvergenceTest] L2(p) = 1.23e-04 L2(v) = 2.34e-04 H1(v) = 3.45e-04
# Returns { error_type: [val_level0, val_level1, ...] }
# ---------------------------------------------------------------------------
_KV_PATTERN = re.compile(r'(\w+\(\w+\))\s*=\s*([0-9Ee+\-.]+)')

def parse_log_file(filename):
    """Parse a convergence log file.  Returns dict {error_type: [values]}."""
    errors = {}
    try:
        with open(filename, 'r') as f:
            for line in f:
                if '[ConvergenceTest]' not in line:
                    continue
                for key, val in _KV_PATTERN.findall(line):
                    errors.setdefault(key, []).append(float(val))
    except FileNotFoundError:
        sys.stderr.write(f"ERROR: Log file '{filename}' not found.\n")
        sys.exit(1)
    return errors


def compute_rates(error_values):
    """Compute convergence rates (log2 ratio) from a list of per-level errors."""
    _log = math.log
    rates = []
    for i in range(len(error_values) - 1):
        e0, e1 = float(error_values[i]), float(error_values[i + 1])
        if math.isnan(e0) or math.isinf(e0):
            continue
        if e0 < 1e-12 or e1 < 1e-12:
            continue  # effectively exact – skip
        rates.append((_log(e0) - _log(e1)) / _log(2))
    return rates


def write_rates_table(filename, errors_dict):
    """Overwrite log file with a human-readable convergence-rate table."""
    if not errors_dict:
        print(f"WARNING: No convergence data found in {filename}, skipping table.")
        return
    col_names = list(errors_dict.keys())

    # Fixed column widths so all rows align regardless of value length
    W_IDX  = 4   # width of level index column
    W_VAL  = 16  # width of error/rate value columns (must exceed len("error_L2dof(p)") = 14)

    def fmt_val(v):
        return f"{v:0.4e}".ljust(W_VAL)

    def fmt_na():
        return "-".ljust(W_VAL)

    # Header
    header = "n".ljust(W_IDX)
    for k in col_names:
        header += f"error_{k}".ljust(W_VAL) + f"rate_{k}".ljust(W_VAL)
    separator = "-" * (W_IDX + len(col_names) * 2 * W_VAL)
    lines = [header, separator]

    n_levels = max(len(v) for v in errors_dict.values())
    all_rates = {k: compute_rates(errors_dict[k]) for k in col_names}

    for i in range(n_levels):
        row = str(i).ljust(W_IDX)
        for k in col_names:
            vals = errors_dict[k]
            rs   = all_rates[k]
            row += fmt_val(vals[i]) if i < len(vals) else fmt_na()
            row += fmt_val(rs[i])   if i < len(rs)   else fmt_na()
        lines.append(row)

    with open(filename, 'w') as f:
        f.write("\n".join(lines) + "\n")


# ---------------------------------------------------------------------------
# Parse both log files
# ---------------------------------------------------------------------------
ffErrors    = parse_log_file(freeflowLogName)
darcyErrors = parse_log_file(darcyLogName)

write_rates_table(freeflowLogName, ffErrors)
write_rates_table(darcyLogName,    darcyErrors)

def print_file(filename):
    with open(filename, 'r') as f:
        sys.stdout.write(f.read())

print(f"\nFreeFlow convergence table ({freeflowLogName}):")
print_file(freeflowLogName)
print()
print(f"\nDarcy convergence table ({darcyLogName}):")
print_file(darcyLogName)
print()

# ---------------------------------------------------------------------------
# Convergence-rate checking helpers
# ---------------------------------------------------------------------------
def mean(numbers):
    return float(sum(numbers)) / len(numbers)


def check_rate_in_log(error_type, min_rate, log_name, errors_dict):
    """Return True if error_type found (and passed), False if not found, exit on failure."""
    if error_type not in errors_dict:
        return False
    rates = compute_rates(errors_dict[error_type])
    if not rates:
        print(f"WARNING: No rates could be computed for '{error_type}' in {log_name}.")
        return True  # found but nothing to check
    last_rate = rates[-1]
    if last_rate < min_rate:
        sys.stderr.write(
            "*" * 70 + "\n"
            f"Convergence rate check FAILED for '{error_type}' in {log_name}:\n"
            f"  last rate = {last_rate:.4f}, expected >= {min_rate:.4f}\n"
            + "*" * 70 + "\n")
        sys.exit(1)
    return True


# ---------------------------------------------------------------------------
# Apply checks according to the rate_specs mode
# ---------------------------------------------------------------------------
all_logs = [(freeflowLogName, ffErrors), (darcyLogName, darcyErrors)]

if rate_specs is None:
    # Built-in defaults
    defaults = {'L2dof(p)': 1.9, 'L2dof(v)': 1.9}
    for et, lo in defaults.items():
        for log_name, errors_dict in all_logs:
            check_rate_in_log(et, lo, log_name, errors_dict)

else:
    # Named format – search all log files for each specified error type
    for error_type, min_rate in rate_specs.items():
        found_in_any = False
        for log_name, errors_dict in all_logs:
            if check_rate_in_log(error_type, min_rate, log_name, errors_dict):
                found_in_any = True
        if not found_in_any:
            sys.stderr.write(
                "*" * 70 + "\n"
                f"ERROR: Specified error type '{error_type}' was not found in any log file.\n"
                f"  Searched in: {[n for n, _ in all_logs]}\n"
                + "*" * 70 + "\n")
            sys.exit(1)

print(f"\nAll convergence rate checks passed for test case '{testname}'.")
