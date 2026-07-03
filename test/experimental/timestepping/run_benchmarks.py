#!/usr/bin/env python3
# SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later
"""Run the time-stepping benchmark and reproduce its figures.

Runs the test_timestepmethods, test_timestepmethods_stabilityregions and
test_timestepmethods_symplectic executables and generates the five figures
shown in README.md via the corresponding plot_*.py scripts. With
--copy-to-doc, the figures are additionally copied into doc/doxygen/images/
for the Doxygen documentation.

Usage:
  python3 run_benchmarks.py [build_dir] [--copy-to-doc]

If build_dir is omitted, the current working directory is used as the
directory containing the test_timestepmethods* executables (e.g. run from
within <build_dir>/test/experimental/timestepping).

Example:
  python3 run_benchmarks.py /path/to/dumux/build-cmake
"""

import argparse
import os
import shutil
import subprocess
import sys


def main() -> int:
    parser = argparse.ArgumentParser(description="Run the time-stepping benchmark and reproduce its figures.")
    parser.add_argument("build_dir", nargs="?", default=None,
                         help="path to the DuMux build directory "
                              "(default: current directory is the executable directory)")
    parser.add_argument("--copy-to-doc", action="store_true", default=False,
                         help="also copy the figures into doc/doxygen/images/ (default: False)")
    args = parser.parse_args()

    src_dir = os.path.dirname(os.path.abspath(__file__))
    if args.build_dir is None:
        bin_dir = os.getcwd()
    else:
        bin_dir = os.path.join(os.path.abspath(args.build_dir), "test", "experimental", "timestepping")

    print(f"\n{'='*60}")
    print("Running test_timestepmethods")
    print(f"{'='*60}")
    subprocess.check_call([os.path.join(bin_dir, "test_timestepmethods")], cwd=bin_dir)
    subprocess.check_call([
        sys.executable, os.path.join(src_dir, "plot_timestepmethods_convergence.py"),
        "--input", os.path.join(bin_dir, "test_timestepmethods_data.json"),
        "--output", os.path.join(bin_dir, "test_timestepmethods_convergence.png"),
    ])

    print(f"\n{'='*60}")
    print("Running test_timestepmethods_stabilityregions")
    print(f"{'='*60}")
    subprocess.check_call([os.path.join(bin_dir, "test_timestepmethods_stabilityregions")], cwd=bin_dir)
    subprocess.check_call([
        sys.executable, os.path.join(src_dir, "plot_timestepmethods_stabilityregions.py"),
        "--input", os.path.join(bin_dir, "test_timestepmethods_stabilityregions_data.json"),
        "--output", os.path.join(bin_dir, "test_timestepmethods_stability.png"),
    ])

    print(f"\n{'='*60}")
    print("Running test_timestepmethods_symplectic")
    print(f"{'='*60}")
    subprocess.check_call([os.path.join(bin_dir, "test_timestepmethods_symplectic")], cwd=bin_dir)
    subprocess.check_call([
        sys.executable, os.path.join(src_dir, "plot_timestepmethods_symplectic.py"),
        "--input", os.path.join(bin_dir, "test_timestepmethods_symplectic_data.json"),
        "--output", os.path.join(bin_dir, "test_timestepmethods_symplectic.png"),
    ])

    images = [
        ("test_timestepmethods_convergence.png", "timestepping_convergence.png"),
        ("test_timestepmethods_stability_explicit.png", "timestepping_stability_explicit.png"),
        ("test_timestepmethods_stability_implicit.png", "timestepping_stability_implicit.png"),
        ("test_timestepmethods_stability_comparison.png", "timestepping_stability_comparison.png"),
        ("test_timestepmethods_symplectic.png", "timestepping_symplectic.png"),
    ]

    print(f"\n{'='*60}")
    print("Figures")
    print(f"{'='*60}")
    for src_name, _ in images:
        print(os.path.join(bin_dir, src_name))

    if args.copy_to_doc:
        doc_images = os.path.normpath(os.path.join(src_dir, "..", "..", "..", "doc", "doxygen", "images"))
        for src_name, dst_name in images:
            src = os.path.join(bin_dir, src_name)
            dst = os.path.join(doc_images, dst_name)
            shutil.copy(src, dst)
            print(f"Copied {src_name} -> {os.path.relpath(dst, src_dir)}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
