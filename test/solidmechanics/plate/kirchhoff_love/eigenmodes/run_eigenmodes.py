# SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

import sys
import subprocess

try:
    import pyvista as pv
except ImportError:
    print("pyvista is required to run this script. Install it with: pip install pyvista")
    sys.exit(77)

subprocess.check_call(["./test_kirchhoff_love_plate_eigenmodes"])

mesh = pv.read("test_kirchhoff_love_plate_eigenmodes-00000.vtu")
plotter = pv.Plotter(shape=(3, 4), window_size=[1600, 1200])
for i in range(12):
    mode_mesh = mesh.copy()
    mode_mesh.clear_data()
    scalar_name = f"h_mode_{i:0>3d}"
    mode_mesh[scalar_name] = mesh[scalar_name]
    plotter.subplot(i // 4, i % 4)
    plotter.add_text(f"Mode {i}", font_size=10)
    plotter.add_mesh(mode_mesh, scalars=scalar_name, cmap="jet", show_scalar_bar=False)
    plotter.view_xy()
plotter.show()
