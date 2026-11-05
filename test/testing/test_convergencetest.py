#!/usr/bin/env python3
# SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

"""
Test the convergence test driver with quadrature errors of known order.

Usage: test_convergencetest.py <path to dumux_convergencetest.py> <path to test_convergencetest_quadrature>
"""

import argparse
import importlib.util
import math
import subprocess
import sys

DRIVER, EXECUTABLE = sys.argv[1], sys.argv[2]
FAILURES = []


def runDriver(name, levels, levelArgs, checks, options=()):
    """Run the driver on the quadrature executable"""
    name = f"test_convergencetest_{name}"
    call = [DRIVER, "--command", f"{EXECUTABLE} {name}", "--errors", f"{name}-errors.csv"]
    call += ["--levels", *levels, "--level-args", levelArgs, *options]
    for check in checks:
        call += ["--check", check]
    return subprocess.run(call, capture_output=True, text=True, check=False)


def expect(name, condition, result=None):
    """Record a failed expectation together with the driver output"""
    if not condition:
        output = f"\n{result.stdout}\n{result.stderr}" if result is not None else ""
        FAILURES.append(f"{name}{output}")


spec = importlib.util.spec_from_file_location("dumux_convergencetest", DRIVER)
driver = importlib.util.module_from_spec(spec)
spec.loader.exec_module(driver)

rates = driver.convergenceRates([1.0, 0.5, 0.25], [1.0, 0.25, 0.0625])
expect("rates of h^2", all(abs(r - 2.0) < 1e-14 for r in rates))
expect("rate of a zero error is NaN", math.isnan(driver.convergenceRates([1.0, 0.5], [1.0, 0.0])[0]))
expect("rate for an unchanged h is NaN", math.isnan(driver.convergenceRates([0.5, 0.5], [1.0, 0.5])[0]))
expect("check with min and max", driver.parseCheck("e finest 1 2") == {"error": "e", "statistic": "finest", "min": 1.0, "max": 2.0})
expect("check without max", driver.parseCheck("e mean 1")["max"] == math.inf)
for invalid in ["e median 1", "e mean", "e mean 1 2 3"]:
    try:
        driver.parseCheck(invalid)
        expect(f"invalid check '{invalid}' is rejected", False)
    except argparse.ArgumentTypeError:
        pass

result = runDriver("pass", ["4", "8", "16", "32"], "{}", ["midpoint all 1.95 2.05", "leftsum all 0.95 1.05"])
expect("midpoint second and left sum first order", result.returncode == 0, result)
expect("table and checks are printed", "midpoint: all rate" in result.stdout and "leftsum: all rate" in result.stdout, result)

result = runDriver("tuple", ["x 4", "x 8"], "{1}", ["midpoint mean 1.95 2.05"])
expect("level with several values", result.returncode == 0, result)

result = subprocess.run(
    [DRIVER, "--command", EXECUTABLE, "--levels", "4", "8", "--level-args", "test_convergencetest_level{0} {0}",
     "--errors", "test_convergencetest_level{0}-errors.csv", "--check", "midpoint mean 1.95 2.05"],
    capture_output=True, text=True, check=False,
)
expect("error file named after the level", result.returncode == 0 and "midpoint: mean rate" in result.stdout, result)

result = runDriver("template", ["4", "8"], "{} {}", [])
expect("level with too few values fails", result.returncode == 1 and "has no value" in result.stderr, result)

result = runDriver("toohigh", ["4", "8", "16"], "{}", ["midpoint mean 2.5"])
expect("rate below its minimum fails", result.returncode == 1 and "midpoint: mean rate" in result.stdout and "NOT in" in result.stdout, result)

result = runDriver("toolow", ["4", "8", "16"], "{}", ["leftsum finest 0.0 0.5"])
expect("rate above its maximum fails", result.returncode == 1 and "NOT in" in result.stdout, result)

result = runDriver("unknown", ["4", "8"], "{}", ["pressure mean 1.0"])
expect("check of an unknown error fails", result.returncode == 1 and "'pressure' is not in" in result.stderr, result)

result = runDriver("nan", ["4", "8"], "{} nan", ["midpoint mean 1.9"])
expect("non-finite checked error fails", result.returncode != 0 and "Non-finite value of midpoint" in result.stderr, result)

result = runDriver("nanunchecked", ["4", "8"], "{} nan", ["leftsum mean 0.9"])
expect("non-finite unchecked error is allowed", result.returncode == 0, result)

result = runDriver("exact", ["4", "8", "16"], "{}", ["linear all 1.9 2.1"])
expect("rates of an exactly integrated error fail without round-off value", result.returncode == 1 and "NOT in" in result.stdout, result)

result = runDriver("roundoff", ["4", "8", "16"], "{}", ["linear all 1.9 2.1", "midpoint all 1.95 2.05"], ["--round-off", "1e-12"])
expect("exactly integrated error passes below the round-off value",
       result.returncode == 0 and "linear: below the round-off value" in result.stdout and "midpoint: all rate" in result.stdout, result)

result = runDriver("roundoffrates", ["4", "8", "16"], "{}", ["midpoint mean 2.5"], ["--round-off", "1e-12"])
expect("error above the round-off value is still checked", result.returncode == 1 and "NOT in" in result.stdout, result)

result = runDriver("roundoffpartial", ["4", "8", "16"], "{}", ["midpoint all 1.95 2.05"], ["--round-off", "2e-3"])
expect("error below the round-off value on some levels only is still checked",
       result.returncode == 0 and "midpoint: all rate" in result.stdout, result)

result = runDriver("crash", ["4", "0"], "{}", ["leftsum mean 0.9"])
expect("exit code of a failing run is passed on", result.returncode == 2, result)

result = subprocess.run(
    [DRIVER, "--command", EXECUTABLE + "_missing", "--errors", "missing.csv", "--levels", "4", "8", "--level-args", "{}"],
    capture_output=True, text=True, check=False,
)
expect("missing executable fails", result.returncode == 1 and "Could not run" in result.stderr, result)

result = subprocess.run(
    [DRIVER, "--command", sys.executable + " -c pass", "--errors", "missing.csv", "--levels", "4", "8", "--level-args", ""],
    capture_output=True, text=True, check=False,
)
expect("missing error file fails", result.returncode == 1 and "did not write the error file" in result.stderr, result)

if FAILURES:
    sys.stderr.write("Failed:\n" + "\n".join(FAILURES) + "\n")
    sys.exit(1)
print("All convergence test driver checks passed")
