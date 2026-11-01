#!/usr/bin/env python3
# SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

"""
Convergence test script for coupled FreeFlow/DarcyFlow problems.

The -Rates argument accepts one or more named error-type specifications:

    -Rates "ErrorType minRate maxRate" ["ErrorType minRate maxRate" ...]

where ErrorType matches a key written to the log file by the simulation,
e.g. L2(p), L2(v), H1(v), H1(p).  The script searches for each specified
error type in ALL log files (freeflow and darcy), unless the error type is
prefixed with "freeflow:" or "darcy:", which restricts the search to that log
file.  If a specified type is not found in any log file the test exits with an
error.  Every observed
convergence rate must satisfy minRate <= rate <= maxRate.  The test also fails
if a simulation fails or an error is not finite.

Example:
    -Rates "freeflow:L2(p) 1.4 1.6" "darcy:L2(p) 1.9 2.1" "L2(v) 1.95 2.1" "H1(v) 0.95 1.1"
"""

import re
import math
import subprocess
import sys

if len(sys.argv) < 2:
    sys.stderr.write('Please provide a single argument <testname> to the script\n')
    sys.exit(1)

executableName = str(sys.argv[1])
testargs = [str(i) for i in sys.argv][2:]
testname = testargs[testargs.index('-Problem.TestCase')+1]
def argument(args, key, default):
    return args[args.index(key)+1] if key in args else default

freeflowLogName = testname + "_" + argument(testargs, '-FreeFlow.Problem.Name', 'freeFlow') + '.log'
darcyLogName = testname + "_" + argument(testargs, '-Darcy.Problem.Name', 'darcy') + '.log'

# ---------------------------------------------------------------------------
# Parse -Rates argument
# ---------------------------------------------------------------------------
def parse_rates(testargs):
    """Return (rate_specs, cleaned_testargs).

    rate_specs is one of:
      None                                     -> use built-in defaults
      dict { error_type: (minRate, maxRate) } -> named format
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

    # named format: each entry is "ErrorType minRate maxRate"
    rate_specs = {}
    for spec in specs_raw:
        parts = spec.split()
        if len(parts) != 3:
            sys.stderr.write(
                f"ERROR: Each -Rates entry must be 'ErrorType minRate maxRate', got: '{spec}'\n")
            sys.exit(1)
        error_type, min_rate, max_rate = parts[0], float(parts[1]), float(parts[2])
        rate_specs[error_type] = (min_rate, max_rate)
    return rate_specs, cleaned


rate_specs, testargs = parse_rates(testargs)

# ---------------------------------------------------------------------------
# Remove old log files
# ---------------------------------------------------------------------------
for logName in [freeflowLogName, darcyLogName]:
    subprocess.call(['rm', '-f', logName])
    print("Removed old log file ({})!".format(logName))

# ---------------------------------------------------------------------------
# Run the simulations at three refinement levels
# ---------------------------------------------------------------------------
for i in [0, 1, 2]:
    subprocess.run(['./' + executableName] + testargs + ['-Grid.Refinement', str(i)], check=True)

# ---------------------------------------------------------------------------
# Generic log-file parser
# Handles lines of the form written by the C++ simulation, e.g.:
#   [ConvergenceTest] L2(p) = 1.23e-04 L2(v) = 2.34e-04 H1(v) = 3.45e-04
# Returns { error_type: [val_level0, val_level1, ...] }
# ---------------------------------------------------------------------------
_KV_PATTERN = re.compile(r'(\w+\(\w+\))\s*=\s*(\S+)')

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
        if not (math.isfinite(e0) and math.isfinite(e1)):
            rates.append(math.nan)
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
def check_rate_in_log(error_type, min_rate, max_rate, log_name, errors_dict):
    """Return True if error_type found (and passed), False if not found, exit on failure."""
    if error_type not in errors_dict:
        return False
    errors = errors_dict[error_type]
    if len(errors) != 3:
        sys.stderr.write(f"ERROR: Expected 3 values of '{error_type}' in {log_name}, found {len(errors)}.\n")
        sys.exit(1)
    rates = compute_rates(errors)
    if not rates:
        if all(e < 1e-12 for e in errors):
            print(f"'{error_type}' in {log_name} is at round-off level on all grids.")
            return True
        sys.stderr.write(f"ERROR: No rates could be computed for '{error_type}' in {log_name}: {errors}\n")
        sys.exit(1)
    if not all(min_rate <= rate <= max_rate for rate in rates):
        sys.stderr.write(
            "*" * 70 + "\n"
            f"Convergence rate check FAILED for '{error_type}' in {log_name}:\n"
            f"  rates = {', '.join(f'{r:.4f}' for r in rates)}, expected in [{min_rate:.4f}, {max_rate:.4f}]\n"
            + "*" * 70 + "\n")
        sys.exit(1)
    return True


# ---------------------------------------------------------------------------
# Apply checks according to the rate_specs mode
# ---------------------------------------------------------------------------
all_logs = [(freeflowLogName, ffErrors), (darcyLogName, darcyErrors)]

if rate_specs is None:
    # Built-in defaults
    defaults = {'L2dof(p)': (1.9, 2.1), 'L2dof(v)': (1.9, 2.1)}
    for et, (lo, hi) in defaults.items():
        for log_name, errors_dict in all_logs:
            check_rate_in_log(et, lo, hi, log_name, errors_dict)

else:
    # Named format – search all log files for each specified error type
    for error_type, (min_rate, max_rate) in rate_specs.items():
        logs = all_logs
        if ':' in error_type:
            domain, error_type = error_type.split(':', 1)
            logs = [{'freeflow': all_logs[0], 'darcy': all_logs[1]}[domain]]
        found_in_any = False
        for log_name, errors_dict in logs:
            if check_rate_in_log(error_type, min_rate, max_rate, log_name, errors_dict):
                found_in_any = True
        if not found_in_any:
            sys.stderr.write(
                "*" * 70 + "\n"
                f"ERROR: Specified error type '{error_type}' was not found in any log file.\n"
                f"  Searched in: {[n for n, _ in logs]}\n"
                + "*" * 70 + "\n")
            sys.exit(1)

print(f"\nAll convergence rate checks passed for test case '{testname}'.")
