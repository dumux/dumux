#!/usr/bin/env python3
# SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later
"""Draw the deflection of the wrinkled sheet from the VTK output of a run.

Usage:
  python3 plot_wrinkles.py <vtu file> [output.png]
"""

import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.tri as tri
import pyvista as pv

mesh = pv.read(sys.argv[1])
output = sys.argv[2] if len(sys.argv) > 2 else "wrinkles.png"

points = mesh.points
cells = mesh.cells.reshape(-1, 4)[:, 1:]
w = mesh.point_data["w"]
triangulation = tri.Triangulation(points[:, 0], points[:, 1], cells)

length = points[:, 0].max() - points[:, 0].min()
width = points[:, 1].max() - points[:, 1].min()
limit = np.percentile(np.abs(w), 99)

fig, ax = plt.subplots(figsize=(7, 7*width/length + 0.8))
image = ax.tripcolor(triangulation, w, shading="gouraud", cmap="RdBu_r", vmin=-limit, vmax=limit)
ax.set_aspect("equal")
ax.set_xlabel("$x/W$")
ax.set_ylabel("$y/W$")
ax.set_title("Sheared sheet: deflection $w/t$")
colorbar = fig.colorbar(image, ax=ax, orientation="vertical", pad=0.03, shrink=0.85)
colorbar.set_label("$w/t$")
fig.tight_layout()
fig.savefig(output, dpi=150)
print(f"Saved {output}")
