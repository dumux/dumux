#!/usr/bin/env python3
# SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later
"""
Runs the complex-valued facet-coupled Helmholtz test on a sequence of grids generated with Gmsh and
checks the convergence rates of the L2 errors in the bulk and in the facet domain.

Usage: convergencetest.py <executable> <name> <minimum rate> [extra arguments for the executable]
The name prefixes the log file and the temporary grid files, so that tests sharing an executable can run concurrently.
"""
import math
import os
import subprocess
import sys


def main():
    exe, name, min_rate, extra = sys.argv[1], sys.argv[2], float(sys.argv[3]), sys.argv[4:]
    if not os.path.isfile(exe):
        sys.exit(77)
    log = name + ".log"
    if os.path.exists(log):
        os.remove(log)
    with open("grids/hybridgrid.geo") as f:
        template = f.read().split("\n", 1)[1]
    for cells in [10, 20, 30, 40, 50, 60, 70, 80]:
        geo = f"grids/{name}.geo"
        with open(geo, "w") as f:
            f.write(f"numElemsPerSide = {cells};\n" + template)
        subprocess.run(["gmsh", "-format", "msh2", "-2", geo], check=True, stdout=subprocess.DEVNULL)
        subprocess.run([f"./{exe}", "params.input", "-Grid.File", f"grids/{name}.msh", "-Grid.NumElemsPerSide", str(cells),
                        "-Problem.OutputFileName", log] + extra, check=True, stdout=subprocess.DEVNULL)
        os.remove(geo)
        os.remove(f"grids/{name}.msh")

    with open(log) as f:
        rows = [[float(v) for v in line.split(",")] for line in f if line.strip()]
    eps = [r[0] for r in rows]
    bulk = [r[1] for r in rows]
    fracture = [r[3] for r in rows]
    for domain, errors in [("matrix", bulk), ("fracture", fracture)]:
        print(f"{domain}:")
        for i in range(len(errors) - 1):
            rate = (math.log(errors[i]) - math.log(errors[i + 1]))/(math.log(eps[i]) - math.log(eps[i + 1]))
            print(f"  a/dx = {eps[i]:.3e}  error = {errors[i]:.4e}  rate = {rate:.4f}")
        print(f"  a/dx = {eps[-1]:.3e}  error = {errors[-1]:.4e}")
    for domain, errors in [("matrix", bulk), ("fracture", fracture)]:
        final = abs((math.log(errors[-2]) - math.log(errors[-1]))/(math.log(eps[-2]) - math.log(eps[-1])))
        if any(math.isnan(e) for e in errors) or final < min_rate:
            sys.exit(f"{domain} convergence rate {final:.3f} below {min_rate}")
    print("check passed")


if __name__ == "__main__":
    main()
