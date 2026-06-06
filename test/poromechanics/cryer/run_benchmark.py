#!/usr/bin/env python3
# SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later
"""Driver that runs the full Cryer benchmark parameter study and produces the figures.

Sweeps the three pressure loads p0/K = {0.0001, 0.25, 0.5} combined with the two
permeability models {Kozeny-Carman, Constant}, i.e. six runs. Each run writes its
centre time series to
  cryer_center_pressure_<perm>_p0k<p0k>.csv
and afterwards plot_cryer.py produces:
  pressure_kozenycarman.png / pressure_constant.png   (normalised centre pressure)
  volume_kozenycarman.png   / volume_constant.png     (volume ratio J = V/V0)
  fig_analytical.png             (Cryer 1963 analytical reference)
  fig_simulation_vs_analytical.png (last run vs analytical)

The material parameters (G, K, porosity, ...) come from params.input in the run
directory; this driver only overrides the load and the permeability model (plus the
time-stepping flags). The default DIRK3 + log-spaced schedule on the moderate mesh
reproduces the published figures but takes a while (~minutes per case); use
--scheme ImplicitEuler and/or --mesh sphere_eighth_coarse.msh for a quick check.

Usage:
  python3 run_benchmark.py <build_dir>
  python3 run_benchmark.py <build_dir> --scheme ImplicitEuler
  python3 run_benchmark.py <build_dir> --skip-run        # only re-plot existing CSVs
"""
import argparse
import os
import shutil
import subprocess
import sys

P0K_VALUES = [0.0001, 0.25, 0.5]
PERM_MODELS = ["KozenyCarman", "Constant"]
BULK_MODULUS = 1.0e6   # must match SpatialParams.BulkModulus in params.input (p0 = p0k * K)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("build_dir", help="DuMux build directory (e.g. .../build-cmake)")
    ap.add_argument("--mesh", default="sphere_eighth.msh")
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
    run_dir = os.path.join(os.path.abspath(args.build_dir),
                           "test", "poromechanics", "cryer")
    exe = os.path.join(run_dir, "test_cryer_poroelastic_large_def")

    if not args.skip_run:
        if not os.path.exists(exe):
            sys.exit(f"Executable not found: {exe}\n"
                     f"Build it first: make test_cryer_poroelastic_large_def")

        for perm in PERM_MODELS:
            for p0k in P0K_VALUES:
                p0 = p0k * BULK_MODULUS
                print(f"\n{'='*60}\n{perm}, p0/K={p0k} (p0={p0:g} Pa)\n{'='*60}")
                subprocess.check_call(
                    [exe,
                     "-Problem.PressureLoad", f"{p0:g}",
                     "-SpatialParams.PermeabilityModel", perm,
                     "-Grid.File", args.mesh,
                     "-TimeLoop.Scheme", args.scheme,
                     "-TimeLoop.Dt", args.dt,
                     "-TimeLoop.DtGrowthFactor", args.growth,
                     "-TimeLoop.MaxTimeStepSize", args.maxdt,
                     "-TimeLoop.TEnd", args.tend],
                    cwd=run_dir,
                )
                src_csv = os.path.join(run_dir, "cryer_center_pressure.csv")
                dst_csv = os.path.join(run_dir, f"cryer_center_pressure_{perm}_p0k{p0k}.csv")
                shutil.copy(src_csv, dst_csv)
                print(f"  -> {os.path.basename(dst_csv)}")

    # Produce all figures (Fig 3-6 + analytical + single comparison) from the CSVs
    print(f"\n{'='*60}\nplotting\n{'='*60}")
    subprocess.check_call([sys.executable, os.path.join(src_dir, "plot_cryer.py")],
                          cwd=run_dir)
    print(f"\nDone. Figures written to {run_dir}")


if __name__ == "__main__":
    main()
