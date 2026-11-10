#!/usr/bin/env python3
# SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later
"""Run the unsteady 2D Navier-Stokes test of the FVCA8 benchmark and compare with the analytical solution.

Usage:
  python3 plot_results.py <executable> [output directory]

Runs the executable with params.input of the current directory on the Cartesian grids of the
benchmark for both viscosities with central and upwind weighting, prints the relative space-time L2
errors with the convergence orders and writes angeli_comparison.png (profiles of velocity and
pressure and the errors) to the output directory (default: current directory). The runs take about
40 minutes on eight cores. Requires numpy, matplotlib and pyvista.
"""

import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

from convergencetest import run, rates, velocityUnknowns, pressureUnknowns, CENTRAL

UPWIND = 1.0
VISCOSITIES = [0.1, 0.01]
# mesh_cart_1 to mesh_cart_5 of the benchmark
CELLS = [4, 8, 16, 32, 64]
POINTWISE_CELLS = [4, 8, 16, 32]
PROFILE_CELLS = 16
PROFILE_VISCOSITY = 0.1
PARALLEL_RUNS = 4


def exactProfile(x, y, viscosity, t):
    """The analytical velocity component in x-direction and pressure on the line through x."""
    import numpy as np

    decay = np.exp(-5.0*viscosity*np.pi**2*t)
    u = -2.0*np.pi*decay*np.cos(np.pi*x)*np.sin(2.0*np.pi*y)
    p = -0.25*decay**2*np.pi**2*(4.0*np.cos(2.0*np.pi*x) + np.cos(4.0*np.pi*y))
    return u, p


def profile(name):
    """Numerical velocity component in x-direction and mean-free pressure at the final time in the column of cells left of x = 0.25."""
    import numpy as np
    import pyvista as pv

    mesh = pv.read(sorted(Path.cwd().glob(f"{name}-*.vtu"))[-1])
    centers = mesh.cell_centers().points
    x = centers[centers[:, 0] < 0.25, 0].max()
    column = np.flatnonzero(np.isclose(centers[:, 0], x))
    column = column[np.argsort(centers[column, 1])]
    p = np.asarray(mesh.cell_data["p"])
    return x, centers[column, 1], np.asarray(mesh.cell_data["velocity_liq (m/s)"])[column, 0], (p - p.mean())[column]


def printTable(title, cells, results):
    pErrors = [r[1] for r in results]
    uErrors = [r[2] for r in results]
    guErrors = [r[4] for r in results]
    uRates = ["-"] + [f"{rate:.2f}" for rate in rates(uErrors, velocityUnknowns(cells))]
    pRates = ["-"] + [f"{rate:.2f}" for rate in rates(pErrors, pressureUnknowns(cells))]
    guRates = ["-"] + [f"{rate:.2f}" for rate in rates(guErrors, velocityUnknowns(cells))]
    print(f"\n{title}")
    print("| Mesh | Cells | errgu | ordgu | erru | ordu | errp | ordp | errdivu | Newton iterations |")
    print("|--:|:--|--:|--:|--:|--:|--:|--:|--:|--:|")
    for i, n in enumerate(cells):
        print(f"| {i + 1} | {n} × {n} | {guErrors[i]:.3e} | {guRates[i]} | {uErrors[i]:.3e} | {uRates[i]} "
              f"| {pErrors[i]:.3e} | {pRates[i]} | {results[i][5]:.1e} | {results[i][3]:.1f} |")


def printComplexityTable(cells, results):
    print("\nComplexity")
    print("| Mesh | Cells | nuu | npu | nnzu | nnzp | nnzup |")
    print("|--:|:--|--:|--:|--:|--:|--:|")
    for i, n in enumerate(cells):
        nonZeros = results[i][6]
        print(f"| {i + 1} | {n} × {n} | {velocityUnknowns(cells)[i]} | {pressureUnknowns(cells)[i]} "
              f"| {nonZeros['nnzu']} | {nonZeros['nnzp']} | {nonZeros['nnzup']} |")


def main():
    if len(sys.argv) < 2:
        sys.stderr.write(__doc__)
        sys.exit(1)

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import numpy as np

    executable = str(Path(sys.argv[1]).resolve())
    outdir = Path(sys.argv[2]).resolve() if len(sys.argv) > 2 else Path.cwd()

    with ThreadPoolExecutor(PARALLEL_RUNS) as pool:
        futures = {(nu, weight): [pool.submit(run, executable, n, nu, weight) for n in CELLS]
                   for nu in VISCOSITIES for weight in (CENTRAL, UPWIND)}
        pointwiseFutures = [pool.submit(run, executable, n, PROFILE_VISCOSITY, CENTRAL, False) for n in POINTWISE_CELLS]
        results = {key: [f.result() for f in values] for key, values in futures.items()}
        pointwise = [f.result() for f in pointwiseFutures]

    for (nu, weight), values in results.items():
        printTable(f"nu = {nu}, {'central' if weight == CENTRAL else 'upwind'}", CELLS, values)
    printTable(f"nu = {PROFILE_VISCOSITY}, central, pointwise initial and boundary velocity", POINTWISE_CELLS, pointwise)
    printComplexityTable(CELLS, results[(PROFILE_VISCOSITY, CENTRAL)])

    fig, axes = plt.subplots(2, 2, figsize=(11, 8))
    colors = {CENTRAL: "#1f6fb4", UPWIND: "#c8501e"}
    names = {CENTRAL: "central", UPWIND: "upwind"}

    # the profiles at the final time on one of the grids with both weightings
    profileIdx = CELLS.index(PROFILE_CELLS)
    tEnd = 1.0/(10.0*PROFILE_VISCOSITY)
    x, y, uCentral, pCentral = profile(results[(PROFILE_VISCOSITY, CENTRAL)][profileIdx][0])
    _, _, uUpwind, pUpwind = profile(results[(PROFILE_VISCOSITY, UPWIND)][profileIdx][0])
    yFine = np.linspace(0.0, 1.0, 400)
    uExact, pExact = exactProfile(x, yFine, PROFILE_VISCOSITY, tEnd)
    for ax, exact, numCentral, numUpwind, label in (
        (axes[0, 0], uExact, uCentral, uUpwind, "$u_x$ [m/s]"),
        (axes[0, 1], pExact, pCentral, pUpwind, "$p$ [Pa]"),
    ):
        ax.plot(yFine, exact, color="black", linewidth=1.5, label="analytical")
        ax.plot(y, numCentral, "o", color=colors[CENTRAL], markerfacecolor="none", markersize=6, label="DuMux, central")
        ax.plot(y, numUpwind, "s", color=colors[UPWIND], markerfacecolor="none", markersize=6, label="DuMux, upwind")
        ax.set_xlabel("$y$ [m]")
        ax.set_ylabel(label)
        ax.set_title(rf"$\nu = {PROFILE_VISCOSITY}$ m²/s, $t = {tEnd:g}$ s, $x = {x:g}$ m, {PROFILE_CELLS} × {PROFILE_CELLS} cells", fontsize=10)
        ax.legend(fontsize=9)

    # the errors for both viscosities on the same scale
    h = np.array([1.0/n for n in CELLS])
    for ax, nu in zip(axes[1], VISCOSITIES):
        for weight in (CENTRAL, UPWIND):
            values = results[(nu, weight)]
            ax.loglog(h, [r[2] for r in values], "-o", color=colors[weight], markersize=5, linewidth=1.5, label=f"erru, {names[weight]}")
            ax.loglog(h, [r[1] for r in values], "--s", color=colors[weight], markersize=5, linewidth=1.5, label=f"errp, {names[weight]}")
        ref = results[(nu, CENTRAL)][-1][2]
        ax.loglog(h, ref*(h/h[-1])**2, color="black", linewidth=0.8, label=r"$\propto h^2$")
        ax.loglog(h, results[(nu, UPWIND)][-1][2]*(h/h[-1]), color="gray", linewidth=0.8, label=r"$\propto h$")
        ax.set_xlabel("$h$ [m]")
        ax.set_ylabel("relative space-time $L^2$ error")
        ax.set_title(rf"$\nu = {nu}$ m²/s, $T = {1.0/(10.0*nu):g}$ s", fontsize=10)
        ax.legend(fontsize=8)
    axes[1, 1].sharey(axes[1, 0])

    for ax in axes.flat:
        ax.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(outdir / "angeli_comparison.png", dpi=150)

    print(f"\nWrote {outdir / 'angeli_comparison.png'}")


if __name__ == "__main__":
    main()
