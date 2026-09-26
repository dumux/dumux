#!/usr/bin/env python3
# SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later
"""Draw the domain and the boundary conditions of the radial injection benchmark.

Usage:
  python3 plot_schematic.py

Produces domain.svg. The interface between the phases follows the similarity solution
of Nordbotten and Celia (2006) for the mobility ratio of the benchmark. The drawing is not to scale.
"""

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import FancyArrowPatch, Polygon, Rectangle

matplotlib.rcParams["mathtext.fontset"] = "cm"
matplotlib.rcParams["font.family"] = "serif"
matplotlib.rcParams["font.serif"] = ["cmr10"]
matplotlib.rcParams["axes.formatter.use_mathtext"] = True
matplotlib.rcParams["svg.fonttype"] = "path"

DIRICHLET_COLOR = "#002ebb"
NEUMANN_COLOR = "#910028"
PLUME_COLOR = "#f1e3b0"
FONT_SIZE = 13

# drawing coordinates: the well at x = 0, the outer boundary at x = LENGTH, the aquifer between y = 0 and y = HEIGHT
LENGTH, HEIGHT = 10.0, 3.0
AXIS = -0.9
TIP = 7.0
MOBILITY_RATIO = 5.11e-4/6.1e-5


def arrow(ax, start, end, color="black", style="-|>", lw=1.2):
    ax.add_patch(FancyArrowPatch(start, end, arrowstyle=style, mutation_scale=12, color=color, lw=lw, shrinkA=0, shrinkB=0))


def dimension(ax, start, end, label, offset, **text_kwargs):
    arrow(ax, start, end, style="<|-|>", lw=1.0)
    center = 0.5*(np.array(start) + np.array(end))
    ax.text(center[0] + offset[0], center[1] + offset[1], label, fontsize=FONT_SIZE, **text_kwargs)


def hatched_wall(ax, y, below):
    ax.plot([0.0, LENGTH], [y, y], color="black", lw=1.5)
    sign = -1.0 if below else 1.0
    for x in np.arange(0.05, LENGTH, 0.2):
        ax.plot([x, x + 0.12], [y + sign*0.02, y + sign*0.2], color="black", lw=0.8)


fig, ax = plt.subplots(figsize=(10.5, 4.3))

# plume of the injected fluid below the top of the aquifer, the interface height follows eq. (14) of Nordbotten and Celia (2006)
x = np.linspace(1e-3, TIP, 400)
thickness = np.clip((TIP/x - 1.0)/(MOBILITY_RATIO - 1.0), 0.0, 1.0)
interface = HEIGHT*(1.0 - thickness)
ax.add_patch(Polygon(np.column_stack([np.concatenate([[0.0], x, [TIP, 0.0]]),
                                      np.concatenate([[0.0], interface, [HEIGHT, HEIGHT]])]),
                     closed=True, facecolor=PLUME_COLOR, edgecolor="none"))
ax.plot(x, interface, color="black", lw=1.0)
ax.add_patch(Rectangle((0.0, 0.0), LENGTH, HEIGHT, fill=False, edgecolor="none"))
ax.text(0.2, 2.2, r"$\mathrm{CO_2}$", fontsize=FONT_SIZE + 2)
ax.text(6.3, 1.1, r"$\mathrm{brine}$", fontsize=FONT_SIZE + 2)

# plume thickness and plume extent
x_h = 2.0
z_h = HEIGHT*(1.0 - np.clip((TIP/x_h - 1.0)/(MOBILITY_RATIO - 1.0), 0.0, 1.0))
dimension(ax, (x_h, HEIGHT), (x_h, z_h), r"$h(r, t)$", (0.12, 0.0))
ax.plot([TIP, TIP], [HEIGHT, HEIGHT + 0.25], color="black", lw=0.8, ls=":")
ax.text(TIP - 0.1, HEIGHT + 0.3, r"$r_p(t)$", fontsize=FONT_SIZE)

# no-flow boundaries at the top and the bottom
hatched_wall(ax, HEIGHT, below=False)
hatched_wall(ax, 0.0, below=True)
ax.text(8.2, HEIGHT + 0.3, r"$\mathrm{no\ flow}$", fontsize=FONT_SIZE)
ax.text(8.2, -0.5, r"$\mathrm{no\ flow}$", fontsize=FONT_SIZE)

# injection well at r = r_w
ax.plot([0.0, 0.0], [0.0, HEIGHT], color=NEUMANN_COLOR, lw=3.5)
for z in np.linspace(0.3, HEIGHT - 0.3, 6):
    arrow(ax, (-0.55, z), (-0.05, z), color=NEUMANN_COLOR)
ax.text(-4.3, 1.75, r"$q_n = \dfrac{\varrho_n Q}{2 \pi r_w H}$", fontsize=FONT_SIZE + 1, color=NEUMANN_COLOR)
ax.text(-4.3, 1.05, r"$q_w = 0$", fontsize=FONT_SIZE + 1, color=NEUMANN_COLOR)

# axis of rotation
ax.plot([AXIS, AXIS], [-0.9, HEIGHT + 0.6], color="black", lw=0.8, ls="-.")
ax.text(AXIS - 0.95, HEIGHT + 0.75, r"$\mathrm{axis\ of\ rotation}$", fontsize=FONT_SIZE - 1)
dimension(ax, (AXIS, HEIGHT + 0.45), (0.0, HEIGHT + 0.45), r"$r_w$", (-0.15, 0.1))

# fixed pressure at the outer radius
ax.plot([LENGTH, LENGTH], [0.0, HEIGHT], color=DIRICHLET_COLOR, lw=3.5)
ax.text(LENGTH + 0.3, 1.75, r"$P_w = 3.15 \cdot 10^{7}\,\mathrm{Pa}$", fontsize=FONT_SIZE + 1, color=DIRICHLET_COLOR)
ax.text(LENGTH + 0.3, 1.05, r"$\bar S_n = 0$", fontsize=FONT_SIZE + 1, color=DIRICHLET_COLOR)

# dimensions, coordinates and gravity
dimension(ax, (LENGTH - 0.3, 0.0), (LENGTH - 0.3, HEIGHT), r"$H = 15\,\mathrm{m}$", (-1.65, -0.05))
dimension(ax, (AXIS, -0.75), (LENGTH, -0.75), r"$R = 2000\,\mathrm{m}$", (-0.7, -0.45))
arrow(ax, (AXIS - 1.8, -0.75), (AXIS - 0.9, -0.75))
arrow(ax, (AXIS - 1.8, -0.75), (AXIS - 1.8, 0.15))
ax.text(AXIS - 0.95, -1.0, r"$r$", fontsize=FONT_SIZE)
ax.text(AXIS - 2.05, 0.1, r"$z$", fontsize=FONT_SIZE)
arrow(ax, (AXIS - 2.6, 0.15), (AXIS - 2.6, -0.75))
ax.text(AXIS - 2.9, -0.35, r"$g$", fontsize=FONT_SIZE)

ax.text(LENGTH + 0.3, -1.25, r"$\mathrm{not\ to\ scale}$", fontsize=FONT_SIZE - 2)

ax.set_xlim(-4.4, LENGTH + 3.6)
ax.set_ylim(-1.35, HEIGHT + 1.05)
ax.set_aspect("equal")
ax.axis("off")
fig.tight_layout()
fig.savefig("domain.svg", bbox_inches="tight", transparent=False, facecolor="white")
print("Saved domain.svg")
