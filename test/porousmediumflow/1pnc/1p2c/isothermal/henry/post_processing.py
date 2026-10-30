#!/usr/bin/env python3
# SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later
"""
Animation of the Henry problem: concentration and 10/50/90 % isochlors of test case 1 (top) and 2 (bottom)
compared with the tables of Fahs et al. (2016).

Run in the build directory after running both tests:  python3 post_processing.py
Reads test_1p2c_henry_diffusion_box.pvd and test_1p2c_henry_dispersion_box.pvd and writes henry_combined.gif.
Requires numpy, matplotlib, pyvista and pillow.
"""

import csv
import io
import shutil

import numpy as np
import pyvista as pv
import matplotlib.pyplot as plt
import matplotlib.tri as mtri
from matplotlib.colors import LinearSegmentedColormap, Normalize
from PIL import Image

# (title, output of the test, reference table)
TEST_CASES = [
    ("Test Case 1 (purely diffusive)", "test_1p2c_henry_diffusion_box.pvd", "fahs2016_table_d1.csv"),
    (r"Test Case 2 ($\alpha_L=0.1$ m, $\alpha_T=0.01$ m)", "test_1p2c_henry_dispersion_box.pvd", "fahs2016_table_d2.csv"),
]
ISOCHLORS = {0.1: (":", "o"), 0.5: ("--", "s"), 0.9: ("-", "^")}  # c: (line style, marker)

COLORMAP = LinearSegmentedColormap.from_list("salinity", ["#4C72B0", "#9E3D22"])
LINE_COLOR = "#4D4D4D"
NUM_FRAMES = 15
FRAME_DURATION = 120  # [ms]


def readReferenceTable(fileName):
    """Fahs et al. (2016) isochlor positions: {c: (x values, z values)}"""
    with open(fileName, newline="", encoding="utf-8") as f:
        rows = list(csv.DictReader(line for line in f if not line.startswith("#")))
    z = [float(row["Z"]) for row in rows]
    return {c: ([float(row[f"X{round(c * 100)}"]) for row in rows], z) for c in ISOCHLORS}


def readTimeSteps(pvdFile, times):
    """Read the output closest to each of the given times"""
    reader = pv.get_reader(pvdFile)
    available = np.asarray(reader.time_values)
    indices = np.unique([np.argmin(np.abs(available - t)) for t in times])
    timeSteps = []
    for i in indices:
        reader.set_active_time_point(i)
        timeSteps.append((available[i], reader.read()[0]))
    return timeSteps


def plotTimeStep(ax, title, time, mesh, reference):
    ax.clear()
    c = np.clip(np.ravel(mesh.point_data["X^Salt_liq"]), 0.0, 1.0)

    # matplotlib can only plot triangles, so each quad is split into two
    quads = mesh.cells_dict[9]
    triangles = mtri.Triangulation(mesh.points[:, 0], mesh.points[:, 1],
                                   np.concatenate([quads[:, [0, 1, 2]], quads[:, [0, 2, 3]]]))

    # levels slightly past [0, 1], otherwise regions with exactly c = 1 are left blank
    ax.tricontourf(triangles, c, levels=np.linspace(-1e-3, 1 + 1e-3, 101), cmap=COLORMAP)
    for level, (lineStyle, marker) in ISOCHLORS.items():
        ax.tricontour(triangles, c, levels=[level], colors=LINE_COLOR, linewidths=1.5, linestyles=lineStyle)
        x, z = reference[level]
        ax.plot(x, z, linestyle="none", marker=marker, markersize=4, markerfacecolor="none",
                markeredgecolor=LINE_COLOR, label=rf"Fahs et al. (2016), $c={level}$")

    ax.set_xlim(mesh.bounds[0], mesh.bounds[1])
    ax.set_ylim(mesh.bounds[2], mesh.bounds[3])
    ax.set_aspect("equal")
    ax.set_xlabel(r"$x$ [m]")
    ax.set_ylabel(r"$z$ [m]")
    ax.set_title(rf"{title} -- $t = {time / 86400:.3f}$ d")


def main():
    if shutil.which("latex") and shutil.which("dvipng"):
        plt.rcParams.update({"text.usetex": True, "font.family": "serif",
                             "font.serif": ["Computer Modern Roman"], "font.size": 12})
    else:
        plt.rcParams.update({"font.family": "serif", "mathtext.fontset": "cm", "font.size": 12})

    # the time steps to plot: evenly spaced up to the end of the shortest run (t = 0 is just seawater)
    endTime = min(pv.get_reader(pvd).time_values[-1] for _, pvd, _ in TEST_CASES)
    times = np.linspace(endTime / NUM_FRAMES, endTime, NUM_FRAMES)
    timeSteps = [readTimeSteps(pvd, times) for _, pvd, _ in TEST_CASES]
    references = [readReferenceTable(table) for _, _, table in TEST_CASES]

    # one plot of both test cases below each other per time step
    fig, axes = plt.subplots(2, 1, figsize=(8.0, 6.0), sharex=True)
    fig.subplots_adjust(left=0.08, right=0.88, bottom=0.10, top=0.92, hspace=0.30)
    colorbar = fig.colorbar(plt.cm.ScalarMappable(Normalize(0, 1), COLORMAP),
                            cax=fig.add_axes([0.90, 0.10, 0.02, 0.82]), ticks=np.linspace(0, 1, 5))
    colorbar.set_label(r"$c = X/X_\mathrm{seawater}$")

    images = []
    for frame in range(min(len(steps) for steps in timeSteps)):
        for ax, (title, _, _), steps, reference in zip(axes, TEST_CASES, timeSteps, references):
            plotTimeStep(ax, title, *steps[frame], reference)
        buffer = io.BytesIO()
        fig.savefig(buffer, dpi=100, format="png")
        buffer.seek(0)
        images.append(Image.open(buffer).convert("RGB"))

    # One palette shared by all frames and no dithering. By default, Pillow reduces each
    # frame to its own 256 colors with dithering, which speckles flat regions and makes
    # the speckles flicker from frame to frame.
    width, height = images[0].size
    allImages = Image.new("RGB", (width, height * len(images)))
    for i, image in enumerate(images):
        allImages.paste(image, (0, i * height))
    palette = allImages.reduce(2).quantize(colors=256, method=Image.Quantize.MEDIANCUT)
    frames = [image.quantize(palette=palette, dither=Image.Dither.NONE) for image in images]
    frames[0].save("henry_combined.gif", save_all=True, append_images=frames[1:], duration=FRAME_DURATION, loop=0)
    print(f"wrote henry_combined.gif ({len(images)} frames)")


if __name__ == "__main__":
    main()
