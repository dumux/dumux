#!/usr/bin/env python3
# SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later
"""Run the Cryer benchmark and plot the results with plot_cryer.py: the large-deformation model
for the loads p0/K = 0.0001, 0.25, 0.5 with constant and with Kozeny-Carman permeability, and the
small-deformation model (test/poromechanics/poroelastic/cryer), whose normalised results do not
depend on the load. The mesh is generated from sphere_eighth.geo in test/poromechanics with gmsh
if it does not exist. The material parameters are those of params.input in the run directories of
the two models. CSV files and figures are written to the run directory of the large-deformation
model.

Usage:
  python3 run_benchmark.py <build_dir>
"""
import argparse
import os
import shutil
import subprocess
import sys

P0K_VALUES = [0.0001, 0.25, 0.5]
PERM_MODELS = ["KozenyCarman", "Constant"]
# must equal SpatialParams.BulkModulus of both params.input files
BULK_MODULUS = 1.0e6


def run(exe, cwd, args, extra):
    subprocess.check_call(
        [exe,
         "-Grid.File", args.meshPath,
         "-TimeLoop.Scheme", args.scheme,
         "-TimeLoop.Dt", args.dt,
         "-TimeLoop.DtGrowthFactor", args.growth,
         "-TimeLoop.MaxTimeStepSize", args.maxdt,
         "-TimeLoop.TEnd", args.tend] + extra,
        cwd=cwd,
    )


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("build_dir", help="DuMux build directory (e.g. .../build-cmake)")
    ap.add_argument("--mesh", default="sphere_eighth.msh")
    ap.add_argument("--lc", default="0.095", help="maximum element size of a generated mesh [m]")
    ap.add_argument("--scheme", default="DIRK3",
                    help="time scheme: ImplicitEuler | CrankNicolson | DIRK3 (default)")
    ap.add_argument("--dt", default="1.0", help="initial time step [s]")
    ap.add_argument("--growth", default="1.3", help="geometric dt growth factor")
    ap.add_argument("--maxdt", default="20.0", help="max time step [s]")
    ap.add_argument("--tend", default="333.33", help="end time [s] (= R0^2/cv)")
    ap.add_argument("--skip-run", action="store_true",
                    help="skip the simulations, only re-plot from existing CSVs")
    args = ap.parse_args()

    src_dir = os.path.dirname(os.path.abspath(__file__))
    poromechanics_dir = os.path.join(os.path.abspath(args.build_dir), "test", "poromechanics")
    run_dir = os.path.join(poromechanics_dir, "poroelastic_large_deformations", "cryer")
    small_dir = os.path.join(poromechanics_dir, "poroelastic", "cryer")
    large_exe = os.path.join(run_dir, "test_cryer_poroelastic_large_def")
    small_exe = os.path.join(small_dir, "test_cryer_poroelastic_small_def")
    args.meshPath = os.path.join(poromechanics_dir, args.mesh)

    if not args.skip_run:
        for exe in (large_exe, small_exe):
            if not os.path.exists(exe):
                sys.exit(f"Executable not found: {exe}\n"
                         f"Build it first: make {os.path.basename(exe)}")

        if not os.path.exists(args.meshPath):
            subprocess.check_call(
                ["gmsh", "-3", "-format", "msh2", "-setnumber", "lc", args.lc,
                 "sphere_eighth.geo", "-o", args.mesh],
                cwd=poromechanics_dir,
            )

        for perm in PERM_MODELS:
            for p0k in P0K_VALUES:
                p0 = p0k * BULK_MODULUS
                print(f"\n{'='*60}\nlarge deformation, {perm}, p0/K={p0k} (p0={p0:g} Pa)\n{'='*60}")
                run(large_exe, run_dir, args,
                    ["-Problem.PressureLoad", f"{p0:g}", "-SpatialParams.PermeabilityModel", perm])
                dst_csv = os.path.join(run_dir, f"cryer_center_pressure_{perm}_p0k{p0k}.csv")
                shutil.copy(os.path.join(run_dir, "cryer_center_pressure.csv"), dst_csv)
                print(f"  -> {os.path.basename(dst_csv)}")

        print(f"\n{'='*60}\nsmall deformation\n{'='*60}")
        p0 = P0K_VALUES[0] * BULK_MODULUS
        run(small_exe, small_dir, args, ["-Problem.PressureLoad", f"{p0:g}"])
        dst_csv = os.path.join(run_dir, "cryer_center_pressure_small_deformation.csv")
        shutil.copy(os.path.join(small_dir, "cryer_center_pressure.csv"), dst_csv)
        print(f"  -> {os.path.basename(dst_csv)}")

    print(f"\n{'='*60}\nplotting\n{'='*60}")
    subprocess.check_call([sys.executable, os.path.join(src_dir, "plot_cryer.py")],
                          cwd=run_dir)
    print(f"\nDone. Figures written to {run_dir}")


if __name__ == "__main__":
    main()
