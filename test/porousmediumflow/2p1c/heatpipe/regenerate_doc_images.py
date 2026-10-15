#!/usr/bin/env python3
# SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later
"""Regenerate the heat pipe benchmark images for the Doxygen documentation.

Builds and runs the box and TPFA tests through plot_convergence.py, then copies
profile.png and convergence.png into doc/doxygen/images/.

Usage:
  python3 regenerate_doc_images.py <build_dir>
"""

import os
import shutil
import subprocess
import sys

if len(sys.argv) < 2:
    sys.stderr.write(__doc__)
    sys.exit(1)

build_dir = os.path.abspath(sys.argv[1])
script_dir = os.path.dirname(os.path.abspath(__file__))
doc_images = os.path.normpath(os.path.join(script_dir, "..", "..", "..", "..", "doc", "doxygen", "images"))
test_build_dir = os.path.join(build_dir, "test", "porousmediumflow", "2p1c", "heatpipe")
executables = ["test_2p1cni_heatpipe_box", "test_2p1cni_heatpipe_tpfa"]

subprocess.check_call(["make", *executables], cwd=test_build_dir)
subprocess.check_call([sys.executable, os.path.join(script_dir, "plot_convergence.py"), *executables], cwd=test_build_dir)

for image, destination in [("profile.png", "heatpipe2p1c_profile.png"),
                           ("convergence.png", "heatpipe2p1c_convergence.png")]:
    shutil.copy(os.path.join(test_build_dir, image), os.path.join(doc_images, destination))
    print(f"Copied {image} to {os.path.join(doc_images, destination)}")
