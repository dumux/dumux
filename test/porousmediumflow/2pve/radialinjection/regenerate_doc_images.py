#!/usr/bin/env python3
# SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later
"""Regenerate the images of the radial injection benchmark for the Doxygen documentation.

Runs plot_schematic.py and plot_convergence.py in the build directory of the test and copies
domain.svg, profile.png, plume.png and convergence.png into doc/doxygen/images/.

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
script_dir = os.path.dirname(os.path.abspath(__file__))
doc_images = os.path.normpath(os.path.join(script_dir, "..", "..", "..", "..", "doc", "doxygen", "images"))
test_build_dir = os.path.join(build_dir, "test", "porousmediumflow", "2pve", "radialinjection")
testname = "test_2pve_radialinjection_tpfa"
image_prefix = "2pve_radialinjection"

subprocess.check_call(["make", testname], cwd=test_build_dir)
subprocess.check_call([sys.executable, os.path.join(script_dir, "plot_schematic.py")], cwd=test_build_dir)
subprocess.check_call([sys.executable, os.path.join(script_dir, "plot_convergence.py"), testname], cwd=test_build_dir)

for image in ("domain.svg", "profile.png", "plume.png", "convergence.png"):
    src = os.path.join(test_build_dir, image)
    dst = os.path.join(doc_images, f"{image_prefix}_{image}")
    shutil.copy(src, dst)
    os.remove(src)
    print(f"Copied {image} -> {os.path.relpath(dst, script_dir)}")

print("\nDone. All images updated in doc/doxygen/images/.")
