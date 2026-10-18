#!/usr/bin/env python3
# SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later
"""Meshes the two half disks, runs the Hertz contact test and checks the contact pressure
against Hertz's solution for the computed force: the peak pressure within 2 % of p0 and the
L2 error of the pressure profile below a/R, the order to which Hertz's half-space solution
describes disks of finite radius.

usage: hertztest.py <gmsh> <executable> <parameter file> [arguments]
"""
import subprocess
import sys

gmsh, executable, params = sys.argv[1:4]
subprocess.run([gmsh, "-2", "hertz.geo", "-o", "hertz.msh"], check=True, stdout=subprocess.DEVNULL)
output = subprocess.run(["./" + executable, params] + sys.argv[4:], check=True, capture_output=True, text=True).stdout
print(output)

rows = [line.split() for line in output.splitlines() if line[:1].isdigit()]
if not rows:
    sys.exit("no results")
for approach, force, aOverR, p0, peakError, profileError, numPairs in rows:
    aOverR, peakError, profileError = float(aOverR), float(peakError), float(profileError)
    if abs(peakError) > 0.02 or profileError > aOverR:
        sys.exit(f"approach {approach}: a/R {aOverR}, peak pressure error {peakError}, profile error {profileError}")
print("pressure agrees with Hertz's solution")
