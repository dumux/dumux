#!/usr/bin/env python3
# SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later
"""
Plot the Henry problem output (*.pvd) as an animated GIF or as a PNG of the final time step.

Pass one .pvd for a single plot, or two (Test Case 1 and 2) for two stacked plots.
The test case, and with it the Fahs et al. (2016) reference table, is taken from the
file name ("case2" in the name -> Table D2, otherwise Table D1).

  post_processing.py test_1p2c_henry_fahs_case1_box.pvd test_1p2c_henry_fahs_case2_box.pvd
  post_processing.py test_1p2c_henry_fahs_case1_box.pvd --out henry_case1_final.png
  post_processing.py adaptive_case1.pvd adaptive_case2.pvd --grid --out henry_adaptive_grid.gif
"""

import argparse
import csv
import io
import os
import shutil

import numpy as np
import pyvista as pv
import matplotlib.pyplot as plt
import matplotlib.tri as mtri
from matplotlib.collections import LineCollection
from matplotlib.colors import LinearSegmentedColormap, Normalize
from PIL import Image
from vtkmodules.util.vtkConstants import VTK_TRIANGLE, VTK_QUAD

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))

# what is plotted
FIELD = "X^solute_liq"
SEAWATER_MASS_FRACTION = 0.035
ISOCHLORS = {0.1: (":", "o"), 0.5: ("--", "s"), 0.9: ("-", "^")}  # c: (line style, marker)

TEST_CASES = {
    1: dict(title="Test Case 1 (purely diffusive)",
            reference="fahs2016_table_d1.csv"),
    2: dict(title=r"Test Case 2 ($\alpha_L=0.1$ m, $\alpha_T=0.01$ m)",
            reference="fahs2016_table_d2.csv"),
}

# how it looks
COLORMAP = LinearSegmentedColormap.from_list("salinity", ["#4C72B0", "#9E3D22"])
LINE_COLOR = "#4D4D4D"
PANEL_SIZE = (8.0, 3.0)  # [inch] per test case
DPI = 200
MESH_LINE_WIDTH = 0.15   # [pt]

# the animation
NUM_FRAMES = 30
FRAME_DURATION = 120     # [ms], GIF stores multiples of 10 ms


def readReferenceTable(testCase):
    """Fahs et al. (2016) isochlor positions: {c: (x values, z values)}"""
    path = os.path.join(SCRIPT_DIR, TEST_CASES[testCase]["reference"])
    with open(path, newline="", encoding="utf-8") as f:
        rows = list(csv.DictReader(line for line in f if not line.startswith("#")))
    z = [float(row["Z"]) for row in rows]
    return {c: ([float(row[f"X{round(c * 100)}"]) for row in rows], z) for c in ISOCHLORS}


# Step 1: pick the time steps to plot and read them
def readTimeSteps(pvdFile, times):
    """Read the output closest to each of the given times."""
    reader = pv.get_reader(pvdFile)
    available = np.asarray(reader.time_values)
    indices = np.unique([np.argmin(np.abs(available - t)) for t in times])
    timeSteps = []
    for i in indices:
        reader.set_active_time_point(i)
        data = reader.read()
        mesh = data.combine() if data.n_blocks > 1 else data[0]
        timeSteps.append((available[i], mesh))
    return timeSteps


def splitIntoTriangles(mesh):
    """matplotlib can only plot triangles, so each quad is split into two (for plotting only)."""
    cells = mesh.cells_dict
    triangles = [cells[VTK_TRIANGLE]] if VTK_TRIANGLE in cells else []
    if VTK_QUAD in cells:
        quads = cells[VTK_QUAD]
        triangles += [quads[:, [0, 1, 2]], quads[:, [0, 2, 3]]]
    x, z = mesh.points[:, 0], mesh.points[:, 1]
    return mtri.Triangulation(x, z, np.concatenate(triangles))


def cellEdges(mesh):
    """The actual cell edges (no quad diagonals), as pairs of point indices."""
    cells = mesh.cells_dict
    edges = []
    for cellType in (VTK_TRIANGLE, VTK_QUAD):
        if cellType in cells:
            corners = cells[cellType]
            edges += [corners[:, [i, (i + 1) % corners.shape[1]]] for i in range(corners.shape[1])]
    return np.unique(np.sort(np.concatenate(edges), axis=1), axis=0)


# Step 2: plot each time step with matplotlib
def plotSolution(ax, mesh, c, reference):
    triangles = splitIntoTriangles(mesh)
    # levels slightly past [0, 1], otherwise regions with exactly c = 1 are left blank
    ax.tricontourf(triangles, c, levels=np.linspace(-1e-3, 1 + 1e-3, 101), cmap=COLORMAP)
    for level, (lineStyle, marker) in ISOCHLORS.items():
        ax.tricontour(triangles, c, levels=[level], colors=LINE_COLOR, linewidths=1.5,
                      linestyles=lineStyle)
        x, z = reference[level]
        ax.plot(x, z, linestyle="none", marker=marker, markersize=4, markerfacecolor="none",
                markeredgecolor=LINE_COLOR, label=rf"Fahs et al. (2016), $c={level}$")


def plotGrid(ax, mesh, c):
    edges = cellEdges(mesh)
    points = mesh.points[:, :2]
    lines = LineCollection(points[edges], cmap=COLORMAP, norm=Normalize(0, 1),
                           linewidths=MESH_LINE_WIDTH)
    lines.set_array(c[edges].mean(axis=1))  # color each edge by its mean concentration
    ax.add_collection(lines)


def plotTimeStep(fig, axes, panels, showGrid):
    """panels: one (test case, time, mesh) per test case, plotted below each other"""
    for ax, (testCase, time, mesh) in zip(axes, panels):
        ax.clear()
        c = np.clip(np.ravel(mesh.point_data[FIELD]) / SEAWATER_MASS_FRACTION, 0.0, 1.0)
        if showGrid:
            plotGrid(ax, mesh, c)
        else:
            plotSolution(ax, mesh, c, readReferenceTable(testCase))
        ax.set_xlim(mesh.bounds[0], mesh.bounds[1])
        ax.set_ylim(mesh.bounds[2], mesh.bounds[3])
        ax.set_aspect("equal")
        ax.set_xlabel(r"$x$ [m]")
        ax.set_ylabel(r"$z$ [m]")
        ax.set_title(rf"{TEST_CASES[testCase]['title']} -- $t = {time / 86400:.3f}$ d, {mesh.n_cells} cells")


def createFigure(numPanels):
    fig, axes = plt.subplots(numPanels, 1, figsize=(PANEL_SIZE[0], numPanels * PANEL_SIZE[1]),
                             sharex=True, squeeze=False)
    bottom, top = (0.10, 0.92) if numPanels > 1 else (0.22, 0.85)
    fig.subplots_adjust(left=0.08, right=0.88, bottom=bottom, top=top, hspace=0.30)
    colorbar = fig.colorbar(plt.cm.ScalarMappable(Normalize(0, 1), COLORMAP),
                            cax=fig.add_axes([0.90, bottom, 0.02, top - bottom]),
                            ticks=np.linspace(0, 1, 5))
    colorbar.set_label(r"$c = X/X_\mathrm{seawater}$")
    return fig, axes[:, 0]


def toImage(fig):
    buffer = io.BytesIO()
    fig.savefig(buffer, dpi=DPI, format="png")
    buffer.seek(0)
    return Image.open(buffer).convert("RGB")


# Step 3: put the plots together into a GIF with Pillow
def writeGif(images, fileName):
    # One palette shared by all frames and no dithering. By default, Pillow reduces each
    # frame to its own 256 colors with dithering, which speckles flat regions and makes
    # the speckles flicker from frame to frame.
    width, height = images[0].size
    allImages = Image.new("RGB", (width, height * len(images)))
    for i, image in enumerate(images):
        allImages.paste(image, (0, i * height))
    palette = allImages.reduce(2).quantize(colors=256, method=Image.Quantize.MEDIANCUT)
    frames = [image.quantize(palette=palette, dither=Image.Dither.NONE) for image in images]
    frames[0].save(fileName, save_all=True, append_images=frames[1:],
                   duration=FRAME_DURATION, loop=0)


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("pvd", nargs="+", help="one or two .pvd files (Test Case 1 and/or 2)")
    parser.add_argument("--grid", action="store_true", help="plot the mesh instead of the solution")
    parser.add_argument("--out", default="henry.gif", help="output file, .gif (animated) or .png (final time step)")
    args = parser.parse_args()

    if shutil.which("latex") and shutil.which("dvipng"):
        plt.rcParams.update({"text.usetex": True, "font.family": "serif",
                             "font.serif": ["Computer Modern Roman"], "font.size": 12})
    else:
        plt.rcParams.update({"font.family": "serif", "mathtext.fontset": "cm", "font.size": 12})

    # Step 1: the time steps to plot -- only the last one for a PNG, otherwise evenly
    # spaced up to the end of the shortest run (t = 0 is skipped, it is just seawater)
    testCases = [2 if "case2" in os.path.basename(f) else 1 for f in args.pvd]
    endTime = min(pv.get_reader(f).time_values[-1] for f in args.pvd)
    animate = args.out.endswith(".gif")
    times = np.linspace(endTime / NUM_FRAMES, endTime, NUM_FRAMES) if animate else [endTime]
    timeSteps = [readTimeSteps(f, times) for f in args.pvd]
    numFrames = min(len(t) for t in timeSteps)

    # Step 2: one matplotlib plot per time step
    fig, axes = createFigure(len(args.pvd))
    images = []
    for frame in range(numFrames):
        panels = [(testCase, *steps[frame]) for testCase, steps in zip(testCases, timeSteps)]
        plotTimeStep(fig, axes, panels, args.grid)
        if animate:
            images.append(toImage(fig))

    # Step 3: write the GIF, or the final time step as PNG
    if animate:
        writeGif(images, args.out)
    else:
        if not args.grid:
            fig.legend(*axes[-1].get_legend_handles_labels(), loc="lower center",
                       bbox_to_anchor=(0.5, -0.04), ncol=3, fontsize=8, frameon=False)
        fig.savefig(args.out, dpi=DPI, bbox_inches="tight")
    print(f"wrote {args.out} ({numFrames} time step{'s' if numFrames > 1 else ''})")


if __name__ == "__main__":
    main()
