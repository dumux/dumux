# SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later
"""Common figure style of the shallow water benchmark scripts.

Series are distinguished by color, marker and line style, so that the figures remain readable in
gray scale. The exact solution is a solid black line without markers.
"""

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker

COLORS = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100"]
MARKERS = ["o", "s", "^", "D"]
LINESTYLES = ["--", "-.", ":", (0, (6, 1.5, 1, 1.5, 1, 1.5))]

plt.rcParams.update({
    "font.size": 13,
    "axes.titlesize": 13,
    "axes.labelsize": 13,
    "legend.fontsize": 12,
    "xtick.labelsize": 12,
    "ytick.labelsize": 12,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "axes.grid": True,
    "grid.color": "0.88",
    "grid.linewidth": 0.6,
    "legend.frameon": False,
    "savefig.dpi": 150,
})


def series(i, points=None):
    """Line properties of the i-th series; with many points, only some of them get a marker"""
    style = {"color": COLORS[i], "marker": MARKERS[i], "ls": LINESTYLES[i], "lw": 1.6, "ms": 6.5,
             "mec": "white", "mew": 0.7}
    if points is not None and points > 15:
        step = max(points//10, 1)
        style["markevery"] = (i*step//len(COLORS), step)
    return style


EXACT = {"color": "black", "ls": "-", "lw": 1.2, "zorder": 1}
REFERENCE = {"color": "0.5", "ls": "-", "lw": 0.8, "zorder": 0}


def slope_triangle(ax, x0, x1, y0, order, above=False):
    """Slope triangle of the given order with the hypotenuse from (x0, y0) to x1 on a log-log plot

    Below a decreasing series, the right angle is at the lower left; above it, at the upper right.
    """
    y1 = y0*(x0/x1)**order
    corner = (x1, y0) if above else (x0, y1)
    ax.loglog([x0, x1, corner[0], x0], [y0, y1, corner[1], y0], **REFERENCE)
    ax.annotate(f"order {order}", xy=(corner[0], (y0*y1)**0.5), xytext=(3 if above else -3, 0),
                textcoords="offset points", ha="left" if above else "right", va="center",
                fontsize=12, color=REFERENCE["color"])


def cell_ticks(ax, cells):
    """Label the logarithmic axis with the numbers of cells"""
    ax.set_xticks(cells, [str(c) for c in cells])
    ax.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
