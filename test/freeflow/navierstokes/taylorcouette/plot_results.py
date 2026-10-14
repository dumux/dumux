#!/usr/bin/env python3
# SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later
"""Run the Taylor-Couette benchmark on four grids and compare with the analytical solution.

Usage:
  python3 plot_results.py <executable> [output directory]

Runs the executable with params.input of the current directory on the grid of params.input and on
three global refinements of it, prints the discrete L2 errors with the convergence rates and writes
taylorcouette_comparison.png (tangential velocity and pressure over the radius) to the output
directory (default: current directory). Requires numpy, matplotlib and pyvista.
"""

import sys
from pathlib import Path

from convergencetest import run, rates

LEVELS = [0, 1, 2, 3]
PLOT_LEVELS = [0, 1]


def radialProfile(vtu):
    """Radius, numerical and analytical tangential velocity and pressure in one column of cells."""
    import numpy as np
    import pyvista as pv

    mesh = pv.read(vtu)
    centers = mesh.cell_centers().points
    angle = np.arctan2(centers[:, 1], centers[:, 0])
    column = np.flatnonzero(np.isclose(angle, angle[angle > 0].min(), rtol=0.0, atol=1e-5))
    column = column[np.argsort(np.hypot(centers[column, 0], centers[column, 1]))]

    x, y = centers[column, 0], centers[column, 1]
    r = np.hypot(x, y)
    tangential = lambda v: (-y*v[column, 0] + x*v[column, 1])/r
    data = lambda name: np.asarray(mesh.cell_data[name])
    return {
        "cells": (len(column), mesh.n_cells // len(column)),
        "r": r,
        "u": tangential(data("velocity_liq (m/s)")),
        "uExact": tangential(data("velocityExact")),
        "p": data("p")[column],
        "pExact": data("pressureExact")[column],
    }


def main():
    if len(sys.argv) < 2:
        sys.stderr.write(__doc__)
        sys.exit(1)

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    executable = str(Path(sys.argv[1]).resolve())
    outdir = Path(sys.argv[2]).resolve() if len(sys.argv) > 2 else Path.cwd()

    results = [run(executable, level) for level in LEVELS]
    profiles = [radialProfile(f"{name}-00000.vtu") for name, _, _ in results]

    pErrors = [r[1] for r in results]
    uErrors = [r[2] for r in results]
    pRates = ["-"] + [f"{rate:.2f}" for rate in rates(pErrors)]
    uRates = ["-"] + [f"{rate:.2f}" for rate in rates(uErrors)]
    print("| Grid (radial × angular cells) | L2 error pressure (Pa) | rate | L2 error velocity (m/s) | rate |")
    print("|:--|--:|--:|--:|--:|")
    for profile, ep, rp, eu, ru in zip(profiles, pErrors, pRates, uErrors, uRates):
        nr, nt = profile["cells"]
        print(f"| {nr} × {nt} | {ep:.3e} | {rp} | {eu:.3e} | {ru} |")

    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))
    finest = profiles[-1]
    for ax, field, label in ((axes[0], "u", r"$u_\theta$ [m/s]"), (axes[1], "p", "$p$ [Pa]")):
        ax.plot(finest["r"], finest[field + "Exact"], color="black", label="analytical")
        for level in PLOT_LEVELS:
            profile = profiles[level]
            nr, nt = profile["cells"]
            ax.plot(profile["r"], profile[field], "o", markerfacecolor="none", markersize=4,
                    label=f"DuMux, {nr} × {nt} cells")
        ax.set_xlabel("$r$ [m]")
        ax.set_ylabel(label)
        ax.grid(alpha=0.3)
        ax.legend(fontsize=9)
    fig.tight_layout()
    fig.savefig(outdir / "taylorcouette_comparison.png", dpi=150)

    print(f"\nWrote {outdir / 'taylorcouette_comparison.png'}")


if __name__ == "__main__":
    main()
