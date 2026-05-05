#!/usr/bin/env python3
# SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later
"""Regenerate plate benchmark images for the Doxygen documentation.

Runs plot_convergence.py for each plate test, then copies the resulting
convergence.png and deformation_3d.png into doc/doxygen/images/.

Usage:
  python3 regenerate_doc_images.py <build_dir>

Example:
  python3 regenerate_doc_images.py /path/to/dumux/build-cmake
"""

import os
import shutil
import subprocess
import sys

if len(sys.argv) < 2:
    sys.stderr.write("Usage: python3 regenerate_doc_images.py <build_dir>\n")
    sys.exit(1)

build_dir = os.path.abspath(sys.argv[1])
plate_src = os.path.dirname(os.path.abspath(__file__))
doc_images = os.path.join(plate_src, "..", "..", "..", "doc", "doxygen", "images")
doc_images = os.path.normpath(doc_images)

tests = [
    ("kirchhoff_love", "test_kirchhoff_love_plate"),
    ("mindlin_reissner", "test_mindlin_reissner_plate"),
    ("membrane", "test_membrane_plate"),
]

for subdir, testname in tests:
    exe = os.path.join(build_dir, "test", "solidmechanics", "plate", subdir, testname)
    src_dir = os.path.join(plate_src, subdir)

    print(f"\n{'='*60}")
    print(f"Running {testname}")
    print(f"{'='*60}")

    subprocess.check_call(
        [sys.executable, os.path.join(plate_src, "plot_convergence.py"), exe],
        cwd=src_dir,
    )

    prefix = f"plate_{subdir}"
    for img, dest in [
        ("convergence.png", f"{prefix}_convergence.png"),
        ("deformation_3d.png", f"{prefix}_deformation_3d.png"),
    ]:
        src = os.path.join(src_dir, img)
        dst = os.path.join(doc_images, dest)
        shutil.copy(src, dst)
        print(f"Copied {img} -> {os.path.relpath(dst, plate_src)}")
        os.remove(src)

print("\nDone. All images updated in doc/doxygen/images/.")
