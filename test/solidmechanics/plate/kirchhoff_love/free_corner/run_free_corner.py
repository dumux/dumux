#!/usr/bin/env python3
# SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later
"""Check corner-force reduction under mesh refinement at convex free corners."""
import os
import re
import subprocess
import sys
from math import isfinite

CASES = [
    ("square", "../../square.geo", [], "1 1  1 0  0 1", {(1.0, 1.0)}),
    ("L-shape", "../../lshape.geo", [], "1 0.5  0.5 1", {(1.0, 0.5), (0.5, 1.0)}),
    ("wedge", "../../wedge.geo", ["-Problem.ClampX0", "false"], "0 1  1 0  0 0", {(0.0, 1.0)}),
]
MESHES = [0.06, 0.03, 0.018]


def run(exe, geo, extra, corners, clmax):
    mesh = "corner.msh"
    subprocess.check_call(["gmsh", "-2", "-format", "msh2", "-clmax", str(clmax), geo,
                           "-o", mesh], stdout=subprocess.DEVNULL)
    out = subprocess.check_output([exe, "-Grid.File", mesh, "-Problem.Corners", corners] + extra,
                                  text=True)
    os.remove(mesh)
    coordinates = [float(value) for value in corners.split()]
    expected = set(zip(coordinates[::2], coordinates[1::2]))
    found = {}
    for line in out.split("\n"):
        m = re.search(r"Corner \((\S+),(\S+)\) (\S+)\s+R = (\S+)\s+ratio = (\S+)", line)
        if m:
            point = (float(m.group(1)), float(m.group(2)))
            force, ratio = float(m.group(4)), float(m.group(5))
            if point in found or not isfinite(force) or not isfinite(ratio) or ratio < 0:
                raise RuntimeError(f"Invalid corner data at {point} for {geo}, clmax={clmax}")
            found[point] = (m.group(3), ratio)
    if set(found) != expected:
        raise RuntimeError(f"Missing or unexpected corners for {geo}, clmax={clmax}: {set(found)}")
    return found


if __name__ == "__main__":
    exe = sys.argv[1] if len(sys.argv) > 1 else "./test_kirchhoff_love_free_corner"
    exe = os.path.abspath(exe)

    failed = False
    for name, geo, extra, corners, convex_free in CASES:
        series = [run(exe, geo, extra, corners, clmax) for clmax in MESHES]
        print(f"\n{name}")
        for point in series[0]:
            kind = series[0][point][0]
            ratios = [s[point][1] for s in series]
            label = "convex free-free" if point in convex_free else kind
            print(f"  corner {str(point):>12} {label:>17}: "
                  + "  ".join(f"{r:8.4f}" for r in ratios))
            if point in convex_free and any(s[point][0] != "free-free" for s in series):
                raise RuntimeError(f"{name}: expected a free-free corner at {point}")
            if point in convex_free and ratios[-1] > 0.6*ratios[0]:
                sys.stderr.write(f"{name}: the corner force at the free corner {point} "
                                 f"is not decreasing under refinement\n")
                failed = True

    if not failed:
        print("\nCorner forces decrease under refinement at all tested convex free corners.")
    sys.exit(1 if failed else 0)
