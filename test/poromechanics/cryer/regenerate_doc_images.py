#!/usr/bin/env python3
# SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later
"""Regenerate the Cryer benchmark images for the Doxygen documentation.

Runs the full parameter study (run_benchmark.py: three loads x two permeability
models) plus a high-fidelity single small-load run for the analytical comparison,
then copies all result figures into doc/doxygen/images/ under their documentation
names:
  cryer_analytical.png, cryer_pressure_vs_analytical.png,
  cryer_pressure_kozenycarman.png, cryer_pressure_constant.png,
  cryer_volume_kozenycarman.png, cryer_volume_constant.png

Note: the runs (sphere_eighth.msh, DIRK3) are expensive (~minutes per case). Pass
--csv-only to skip all simulations and just (re)plot/copy from existing CSVs.

Usage:
  python3 regenerate_doc_images.py <build_dir> [--csv-only]
"""

import argparse
import os
import shutil
import subprocess
import sys

parser = argparse.ArgumentParser()
parser.add_argument("build_dir", help="DuMux build directory (e.g. .../build-cmake)")
parser.add_argument("--csv-only", action="store_true",
                    help="skip the simulations, only (re)plot and copy")
args = parser.parse_args()

build_dir = os.path.abspath(args.build_dir)
src_dir = os.path.dirname(os.path.abspath(__file__))
doc_images = os.path.normpath(
    os.path.join(src_dir, "..", "..", "..", "doc", "doxygen", "images"))

run_dir = os.path.join(build_dir, "test", "poromechanics", "cryer")
exe = os.path.join(run_dir, "test_cryer_poroelastic_large_def")

# 1. Full parameter study -> the four pressure/volume figures (+ analytical)
study_cmd = [sys.executable, os.path.join(src_dir, "run_benchmark.py"), build_dir]
if args.csv_only:
    study_cmd.append("--skip-run")
subprocess.check_call(study_cmd)

# 2. High-fidelity single small-load run (DIRK3, fine log-spaced steps) for the
#    analytical comparison figure. Run last so cryer_center_pressure.csv holds it.
if not args.csv_only:
    print("\nHigh-fidelity single run for the analytical comparison...")
    subprocess.check_call(
        [exe,
         "-Problem.PressureLoad", "100.0",
         "-SpatialParams.PermeabilityModel", "Constant",
         "-Grid.File", "sphere_eighth.msh",
         "-TimeLoop.Scheme", "DIRK3",
         "-TimeLoop.Dt", "1.0",
         "-TimeLoop.DtGrowthFactor", "1.15",
         "-TimeLoop.MaxTimeStepSize", "8.0"],
        cwd=run_dir,
    )
subprocess.check_call(
    [sys.executable, os.path.join(src_dir, "plot_cryer.py")], cwd=run_dir)

# 3. Copy all figures into doc/doxygen/images/ under their documentation names
images = [
    ("fig_simulation_vs_analytical.png", "cryer_pressure_vs_analytical.png"),
    ("pressure_kozenycarman.png",       "cryer_pressure_kozenycarman.png"),
    ("pressure_constant.png",           "cryer_pressure_constant.png"),
    ("volume_kozenycarman.png",         "cryer_volume_kozenycarman.png"),
    ("volume_constant.png",             "cryer_volume_constant.png"),
]
for img, dest in images:
    src = os.path.join(run_dir, img)
    dst = os.path.join(doc_images, dest)
    shutil.copy(src, dst)
    print(f"Copied {img} -> {os.path.relpath(dst, src_dir)}")

# 4. Boundary-condition schematic (no simulation needed)
print("\nDrawing the setup schematic...")
setup_svg = os.path.join(doc_images, "cryer_setup.svg")
subprocess.check_call(
    [sys.executable, os.path.join(src_dir, "make_setup_svg.py"), "--out", setup_svg])
print(f"Wrote {os.path.relpath(setup_svg, src_dir)}")

print("\nDone. Cryer images updated in doc/doxygen/images/.")
