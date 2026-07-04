#!/usr/bin/env python3
# SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later
"""Regenerate the interactive 3D deformation scenes for the plate benchmark docs.

For each plate test this runs the (already built) executable on a single fine
mesh, converts the resulting VTK output into a small interactive-scene data
file via ``bin/postprocessing/vtk_to_scene.py``, and writes it into
``doc/doxygen/scenes/`` as ``scene_plate_<name>.js`` (geometry embedded in a
script so the scene also renders from ``file://``).

The generated files are committed to the repository; Doxygen copies them into
its html output (see the ``HTML_EXTRA_FILES`` entries in ``doc/doxygen/Doxylocal``)
where the plate READMEs embed them, rendered in the browser by
``doc/doxygen/dumux-scene.js`` using vtk.js. Run this by hand to refresh them:

  python3 regenerate_doc_scenes.py <build_dir>

This must be run with a Python interpreter that has ``pyvista`` installed; that
same interpreter is used for the scene conversion.
"""

import os
import subprocess
import sys

# One entry per plate benchmark: build subdirectory, executable name, and the
# scalar field to color by / warp the geometry with (vertical deformation).
TESTS = [
    ("kirchhoff_love", "test_kirchhoff_love_plate", "w"),
    ("mindlin_reissner", "test_mindlin_reissner_plate", "w"),
    ("membrane", "test_membrane_plate", "w"),
]

# Mesh resolution (Gmsh characteristic length) used for the scene. Fine enough
# to look smooth, coarse enough to keep the embedded scene small.
CLMAX = "0.0625"


def generate_scene(build_dir, out_dir, subdir, testname, scalar, converter):
    """Run one plate test on a single mesh and export its deformation scene."""
    build_subdir = os.path.join(build_dir, "test", "solidmechanics", "plate", subdir)
    exe = os.path.join(build_subdir, testname)
    if not os.path.isfile(exe):
        print(f"  skipping {testname}: executable not found (build it first)")
        return False

    # The executable reads 'disk.msh' from the working directory and 'disk.geo'
    # is symlinked into the parent build directory (../disk.geo), mirroring the
    # convergence test setup.
    print(f"\n{'='*60}\nGenerating scene for {testname}\n{'='*60}")
    subprocess.check_call(
        ["gmsh", "-2", "-format", "msh2", "-clmax", CLMAX, "../disk.geo", "-o", "disk.msh"],
        cwd=build_subdir,
    )
    subprocess.check_call([exe], cwd=build_subdir)

    pvd = os.path.join(build_subdir, testname + ".pvd")
    scene = os.path.join(out_dir, f"scene_plate_{subdir}.js")
    subprocess.check_call(
        [
            sys.executable, converter, pvd, scene,
            "--scalar", scalar,
            "--warp-scalar", scalar,
        ]
    )
    print(f"Wrote {os.path.relpath(scene, out_dir)}")

    # Clean up the VTK output and mesh generated in the build directory
    for name in os.listdir(build_subdir):
        if name.startswith(testname) and (name.endswith(".vtu") or name.endswith(".pvd")):
            os.remove(os.path.join(build_subdir, name))
    msh = os.path.join(build_subdir, "disk.msh")
    if os.path.exists(msh):
        os.remove(msh)
    return True


def main(argv):
    if len(argv) < 2:
        sys.stderr.write("Usage: python3 regenerate_doc_scenes.py <build_dir>\n")
        return 1

    build_dir = os.path.abspath(argv[1])

    # Resolve symlinks so this works whether run from the source tree or from a
    # symlinked copy in the build tree.
    plate_src = os.path.dirname(os.path.realpath(__file__))
    out_dir = os.path.normpath(os.path.join(plate_src, "..", "..", "..", "doc", "doxygen", "scenes"))
    os.makedirs(out_dir, exist_ok=True)
    converter = os.path.normpath(
        os.path.join(plate_src, "..", "..", "..", "bin", "postprocessing", "vtk_to_scene.py")
    )

    generated = 0
    for subdir, testname, scalar in TESTS:
        if generate_scene(build_dir, out_dir, subdir, testname, scalar, converter):
            generated += 1

    print(f"\nDone. Generated {generated}/{len(TESTS)} scenes in {out_dir}.")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
