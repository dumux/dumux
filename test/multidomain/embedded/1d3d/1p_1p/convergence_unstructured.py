#!/usr/bin/env python3
# SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later
"""Kernel coupling on unstructured tetrahedral meshes: conservation and error decay.

The meshes are generated here rather than stored: gmsh is called once per resolution and the
.msh is removed afterwards. The executable itself throws if the kernel is not conservative,
so this script adds the part a single run cannot show, namely that the error decays at the
expected rate under refinement.
"""
import math
import os
import subprocess
import sys

exe = sys.argv[1] if len(sys.argv) > 1 else "./test_md_embedded_1d3d_1p1p_boxtpfa_kernel_unstructured"
params = sys.argv[2] if len(sys.argv) > 2 else "params_unstructured.input"
sizes = [0.30, 0.20, 0.15]

errors = []
for clmax in sizes:
    subprocess.run(
        ["gmsh", "-3", "-format", "msh2", "-clmax", str(clmax), "box.geo", "-o", "box.msh"],
        check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
    )
    out = subprocess.check_output([exe, params, "-Tissue.Grid.File", "box.msh"], text=True)
    for line in out.splitlines():
        if "[conservation]" in line:
            print(f"  clmax={clmax}: {line.split('rel. error =')[1].strip()}")
        if " L2_p3d:" in line:
            errors.append(float(line.split()[1]))

if os.path.exists("box.msh"):
    os.remove("box.msh")

if len(errors) < 2:
    sys.stderr.write("No error norms were produced\n")
    sys.exit(1)

rate = math.log(errors[-2] / errors[-1]) / math.log(sizes[-2] / sizes[-1])
print("  L2_p3d: " + " ".join(f"{e:.4e}" for e in errors) + f"   finest-pair rate {rate:.2f}")
if rate < 1.5:
    sys.stderr.write(f"Convergence rate {rate:.2f} is below the expected second order\n")
    sys.exit(1)
print("Kernel coupling is conservative on unstructured meshes and converges as expected.")
