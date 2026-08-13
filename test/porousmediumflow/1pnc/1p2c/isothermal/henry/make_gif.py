#!/usr/bin/env python3
# SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later
"""
Runs a Henry problem test case (see README.md) and renders the transient approach to
steady state as an animated GIF, together with a final static comparison against the
digitized semianalytical isochlor positions of Fahs et al. (2016, WRR,
doi:10.1002/2016WR019288).

pyvista is used purely to read the VTU/PVD time series and extract the salt mass
fraction field (no off-screen 3D rendering is used, since the domain is a flat 2D
slab -- matplotlib, with real LaTeX text rendering, does all of the actual plotting).

Run from the build directory (the executables and symlinked *.input files must be in
the current working directory), e.g.:

    cd <build-dir>/test/porousmediumflow/1pnc/1p2c/isothermal/henry
    python3 <source-dir>/test/porousmediumflow/1pnc/1p2c/isothermal/henry/make_gif.py --combined

--combined runs both cases and renders them stacked into a single GIF/PNG (this is
what the README embeds): a single image element animates in lockstep by
construction, whereas two independent GIF files drift out of sync in a browser
regardless of matching frame timing, since each one's play head starts on its own
load/decode schedule. --case 1/--case 2 remain available for quick single-case
iteration (their own separate GIF/PNG).
"""

import argparse
import glob
import os
import shutil
import shlex
import subprocess
import sys

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib.colors import LinearSegmentedColormap

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from validate_fahs2016 import readReferenceTable

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))

CASES = {
    1: dict(
        executable="test_1p2c_henry_fahs_box",
        params="params.input",
        reference=os.path.join(SCRIPT_DIR, "fahs2016_table_d1.csv"),
        title="Henry Problem, Test Case 1 (purely diffusive)",
        out_prefix="henry_case1",
    ),
    2: dict(
        executable="test_1p2c_henry_fahs_case2_box",
        params="params_case2.input",
        reference=os.path.join(SCRIPT_DIR, "fahs2016_table_d2.csv"),
        title=r"Henry Problem, Test Case 2 ($\alpha_L=0.1$ m, $\alpha_T=0.01$ m)",
        out_prefix="henry_case2",
    ),
}

FIELD = "X^solute_liq"
SEAWATER_FRACTION = 0.035
LEVELS = (10, 50, 90)

# Both cases are already converged well before this (see README); this equals
# TimeLoop.TEnd exactly, so the animation covers the full run for both cases. Both
# cases animate over this same [0, ANIMATION_DURATION] window with the same number of
# frames, so the two GIFs play in lockstep in real simulated time.
ANIMATION_DURATION = 86400.0  # [s], 1 d

# freshwater (blue) -> brine (red-brown), interpolated directly in RGB so the
# midpoint is a muted plum rather than washing out through white
SALINITY_CMAP = LinearSegmentedColormap.from_list(
    "henry_salinity", ["#4C72B0", "#9E3D22"]
)

# all three isochlors are compared against the same reference color; they are
# told apart by linestyle/marker shape instead
COMPARISON_COLOR = "#4D4D4D"
LEVEL_STYLE = {
    10: dict(linestyle=":", marker="o"),
    50: dict(linestyle="--", marker="s"),
    90: dict(linestyle="-", marker="^"),
}


def setupLatexFonts():
    """Use real LaTeX text rendering if a LaTeX toolchain is available, otherwise
    fall back to matplotlib's built-in (LaTeX-like) Computer Modern mathtext."""
    if shutil.which("latex") and shutil.which("dvipng"):
        plt.rcParams.update(
            {
                "text.usetex": True,
                "font.family": "serif",
                "font.serif": ["Computer Modern Roman"],
                "font.size": 12,
            }
        )
    else:
        print("No LaTeX toolchain found (latex/dvipng); falling back to mathtext.")
        plt.rcParams.update(
            {
                "text.usetex": False,
                "font.family": "serif",
                "mathtext.fontset": "cm",
                "font.size": 12,
            }
        )


def runSimulation(executable, params):
    if not os.path.exists(executable):
        sys.stderr.write(
            f"'{executable}' not found in {os.getcwd()}. Build it first, "
            f"e.g. `cmake --build . --target {executable}`.\n"
        )
        sys.exit(1)

    # remove any stale output from a previous run so the frame list below is unambiguous
    for f in glob.glob(f"{executable}-*.vtu") + [f"{executable}.pvd"]:
        os.remove(f)

    result = subprocess.call(shlex.split(f"./{executable} {params}"))
    if result:
        sys.exit(result)


def loadFrames(executable, numFrames, duration):
    """Returns a list of (time, points[:,:2], field values) tuples, sampled at
    numFrames equally spaced times in [0, duration] (rather than equally spaced
    output-frame indices), so that two calls with the same numFrames/duration --
    e.g. for different test cases with different output time steps -- animate at
    the same real-time speed and stop at the same simulated time."""
    import pyvista as pv

    pvd = f"{executable}.pvd"
    reader = pv.get_reader(pvd)
    timeValues = np.asarray(reader.time_values)

    targetTimes = np.linspace(0.0, duration, numFrames)
    indices = np.unique([np.argmin(np.abs(timeValues - t)) for t in targetTimes])

    frames = []
    for i in indices:
        reader.set_active_time_value(timeValues[i])
        mesh = reader.read()[0]
        points = mesh.points[:, :2]
        values = np.ravel(mesh.point_data[FIELD])
        frames.append((timeValues[i], points, values))
    return frames


def toStructuredGrid(points, values):
    """Reshape an unstructured (but grid-aligned) point cloud from a structured
    YaspGrid into 2D arrays suitable for contourf/contour, regardless of the VTU's
    native node ordering."""
    order = np.lexsort((points[:, 0], points[:, 1]))
    xs, ys, vs = points[order, 0], points[order, 1], values[order]
    nx, ny = len(np.unique(xs)), len(np.unique(ys))
    return xs.reshape(ny, nx), ys.reshape(ny, nx), vs.reshape(ny, nx)


def drawFrame(ax, X, Y, c, time, title, reference):
    ax.clear()
    # clip (physically the field is bounded by [0, 1] anyway) rather than using
    # contourf's extend="both", so the colorbar is a plain rectangle without
    # triangular over/under arrows
    c = np.clip(c, 0.0, 1.0)
    contourf = ax.contourf(X, Y, c, levels=np.linspace(0, 1, 21), cmap=SALINITY_CMAP)

    contourSet = None
    for lvl in LEVELS:
        contourSet = ax.contour(
            X, Y, c, levels=[lvl / 100.0],
            colors=COMPARISON_COLOR, linewidths=1.5,
            linestyles=LEVEL_STYLE[lvl]["linestyle"],
        )

    for lvl in LEVELS:
        zRef = [row["Z"] for row in reference]
        xRef = [row[lvl] for row in reference]
        ax.plot(
            xRef, zRef, linestyle="none", marker=LEVEL_STYLE[lvl]["marker"], markersize=4,
            markerfacecolor="none", markeredgecolor=COMPARISON_COLOR,
            label=rf"Fahs et al. (2016), $c={lvl/100:.1f}$",
        )

    ax.set_xlim(X.min(), X.max())
    ax.set_ylim(Y.min(), Y.max())
    ax.set_aspect("equal")
    ax.set_xlabel(r"$x$ [m]")
    ax.set_ylabel(r"$z$ [m]")
    days = time / 86400.0
    ax.set_title(title + rf" -- $t = {days:.2f}$ d")

    return contourf, contourSet


def makeFigures(case):
    cfg = CASES[case]
    reference = readReferenceTable(cfg["reference"])

    frames = loadFrames(cfg["executable"], numFrames=60, duration=ANIMATION_DURATION)

    fig, ax = plt.subplots(figsize=(8, 4.0))
    fig.subplots_adjust(left=0.08, right=0.88, bottom=0.28, top=0.88)
    cbar_ax = fig.add_axes([0.90, 0.28, 0.02, 0.60])
    cbar = None

    def update(frameIdx):
        nonlocal cbar
        time, points, values = frames[frameIdx]
        X, Y, C = toStructuredGrid(points, values / SEAWATER_FRACTION)
        contourf, _ = drawFrame(ax, X, Y, C, time, cfg["title"], reference)
        if cbar is None:
            cbar = fig.colorbar(contourf, cax=cbar_ax, ticks=np.linspace(0, 1, 5))
            cbar.set_label(r"$c = X/X_\mathrm{seawater}$")
        return ()

    anim = animation.FuncAnimation(fig, update, frames=len(frames), blit=False)
    gifPath = f"{cfg['out_prefix']}.gif"
    anim.save(gifPath, writer=animation.PillowWriter(fps=8))
    print(f"Saved {gifPath}")

    # final frame as a static comparison figure
    update(len(frames) - 1)
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.22), ncol=3, fontsize=9, frameon=False)
    pngPath = f"{cfg['out_prefix']}_final.png"
    fig.savefig(pngPath, dpi=200)
    print(f"Saved {pngPath}")


def makeCombinedFigure():
    """Renders both cases stacked into a single GIF/PNG, sharing frame times and a
    colorbar, so the two panels are structurally in sync -- unlike two separate GIF
    files, which a browser plays independently (and thus out of sync) regardless of
    matching frame timing."""
    references = {case: readReferenceTable(cfg["reference"]) for case, cfg in CASES.items()}
    frames = {
        case: loadFrames(cfg["executable"], numFrames=60, duration=ANIMATION_DURATION)
        for case, cfg in CASES.items()
    }
    numFrames = min(len(frames[1]), len(frames[2]))

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(8, 6.4), sharex=True)
    fig.subplots_adjust(left=0.08, right=0.88, bottom=0.16, top=0.92, hspace=0.35)
    cbar_ax = fig.add_axes([0.90, 0.16, 0.02, 0.76])
    cbar = None

    def update(frameIdx):
        nonlocal cbar
        for ax, case in ((ax1, 1), (ax2, 2)):
            time, points, values = frames[case][frameIdx]
            X, Y, C = toStructuredGrid(points, values / SEAWATER_FRACTION)
            contourf, _ = drawFrame(ax, X, Y, C, time, CASES[case]["title"], references[case])
        if cbar is None:
            cbar = fig.colorbar(contourf, cax=cbar_ax, ticks=np.linspace(0, 1, 5))
            cbar.set_label(r"$c = X/X_\mathrm{seawater}$")
        return ()

    anim = animation.FuncAnimation(fig, update, frames=numFrames, blit=False)
    gifPath = "henry_combined.gif"
    anim.save(gifPath, writer=animation.PillowWriter(fps=8))
    print(f"Saved {gifPath}")

    # final frame as a static comparison figure
    update(numFrames - 1)
    handles, labels = ax2.get_legend_handles_labels()
    fig.legend(handles, labels, loc="lower center", bbox_to_anchor=(0.5, 0.0), ncol=3, fontsize=9, frameon=False)
    pngPath = "henry_combined_final.png"
    fig.savefig(pngPath, dpi=200)
    print(f"Saved {pngPath}")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--case", type=int, choices=sorted(CASES), help="Render only this case's own GIF/PNG (quick single-case iteration)")
    group.add_argument("--combined", action="store_true", help="Run both cases and render them stacked into a single, structurally in-sync GIF/PNG (what the README embeds)")
    parser.add_argument("--skip-run", action="store_true", help="Reuse existing VTU/PVD output instead of re-running the simulation(s)")
    args = parser.parse_args()

    try:
        import pyvista  # noqa: F401
    except ImportError:
        sys.stderr.write("pyvista is required for this script: pip install pyvista\n")
        sys.exit(1)

    setupLatexFonts()

    if args.combined:
        if not args.skip_run:
            for cfg in CASES.values():
                runSimulation(cfg["executable"], cfg["params"])
        makeCombinedFigure()
    else:
        cfg = CASES[args.case]
        if not args.skip_run:
            runSimulation(cfg["executable"], cfg["params"])
        makeFigures(args.case)


if __name__ == "__main__":
    main()
