#!/usr/bin/env python3
# SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later
"""Run the unsteady 2D Navier-Stokes test of the FVCA8 benchmark and compare with the solution.

Usage:
  python3 plot_results.py <executable> [output directory]

Runs the executable with params.input of the current directory on the Cartesian grids of the
benchmark for both viscosities with central and upwind weighting, prints the relative space-time L2
errors with the convergence orders and writes angeli_comparison.png (profiles of velocity and
pressure and the errors) to the output directory (default: current directory). The runs take about
40 minutes on eight cores. Requires numpy, matplotlib and pyvista.
"""

import sys
from concurrent import futures
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pyvista as pv

from convergencetest import run, rates, velocityUnknowns, pressureUnknowns, CENTRAL

UPWIND = 1.0
VISCOSITIES = [0.1, 0.01]
# mesh_cart_1 to mesh_cart_5 of the benchmark
CELLS = [4, 8, 16, 32, 64]
POINTWISE_CELLS = [4, 8, 16, 32]
PROFILE_CELLS = 16
PROFILE_VISCOSITY = 0.1
PARALLEL_RUNS = 4
COLORS = {CENTRAL: "#1f6fb4", UPWIND: "#c8501e"}
NAMES = {CENTRAL: "central", UPWIND: "upwind"}


def exactProfile(x, y, viscosity, t):
    """The analytical velocity component in x-direction and pressure on the line through x."""
    decay = np.exp(-5.0 * viscosity * np.pi**2 * t)
    u = -2.0 * np.pi * decay * np.cos(np.pi * x) * np.sin(2.0 * np.pi * y)
    p = -0.25 * decay**2 * np.pi**2 * (4.0 * np.cos(2.0 * np.pi * x) + np.cos(4.0 * np.pi * y))
    return u, p


def profile(name):
    """Numerical velocity in x-direction and mean-free pressure at the final time.

    The values are taken in the column of cells left of x = 0.25.
    """
    mesh = pv.read(sorted(Path.cwd().glob(f"{name}-*.vtu"))[-1])
    centers = mesh.cell_centers().points
    x = centers[centers[:, 0] < 0.25, 0].max()
    column = np.flatnonzero(np.isclose(centers[:, 0], x))
    column = column[np.argsort(centers[column, 1])]
    p = np.asarray(mesh.cell_data["p"])
    return (
        x,
        centers[column, 1],
        np.asarray(mesh.cell_data["velocity_liq (m/s)"])[column, 0],
        (p - p.mean())[column],
    )


def formatRates(errors, unknowns):
    """The orders as strings for a table, starting with a dash for the coarsest grid."""
    return ["-"] + [f"{rate:.2f}" for rate in rates(errors, unknowns)]


def printTable(title, cells, results):
    """Print the accuracy table of the benchmark for one case."""
    pErrors = [r.pressureError for r in results]
    uErrors = [r.velocityError for r in results]
    guErrors = [r.gradientError for r in results]
    pRates = formatRates(pErrors, pressureUnknowns(cells))
    uRates = formatRates(uErrors, velocityUnknowns(cells))
    guRates = formatRates(guErrors, velocityUnknowns(cells))
    print(f"\n{title}")
    print(
        "| Mesh | Cells | errgu | ordgu | erru | ordu | errp | ordp | errdivu | Newton iterations |"
    )
    print("|--:|:--|--:|--:|--:|--:|--:|--:|--:|--:|")
    for i, n in enumerate(cells):
        print(
            f"| {i + 1} | {n} × {n} | {guErrors[i]:.3e} | {guRates[i]} "
            f"| {uErrors[i]:.3e} | {uRates[i]} | {pErrors[i]:.3e} | {pRates[i]} "
            f"| {results[i].divergenceError:.1e} | {results[i].newtonIterations:.1f} |"
        )


def printComplexityTable(cells, results):
    """Print the complexity table of the benchmark."""
    print("\nComplexity")
    print("| Mesh | Cells | nuu | npu | nnzu | nnzp | nnzup |")
    print("|--:|:--|--:|--:|--:|--:|--:|")
    for i, n in enumerate(cells):
        nonZeros = results[i].nonZeros
        print(
            f"| {i + 1} | {n} × {n} | {velocityUnknowns(cells)[i]} | {pressureUnknowns(cells)[i]} "
            f"| {nonZeros['nnzu']} | {nonZeros['nnzp']} | {nonZeros['nnzup']} |"
        )


def runStudy(executable):
    """Run all cases in parallel and return the results by viscosity and weighting.

    Also returns the results with pointwise initial and boundary velocities.
    """
    with futures.ThreadPoolExecutor(PARALLEL_RUNS) as pool:
        submitted = {
            (nu, weight): [pool.submit(run, executable, n, nu, weight) for n in CELLS]
            for nu in VISCOSITIES
            for weight in (CENTRAL, UPWIND)
        }
        pointwiseSubmitted = [
            pool.submit(run, executable, n, PROFILE_VISCOSITY, CENTRAL, False)
            for n in POINTWISE_CELLS
        ]
        results = {key: [f.result() for f in values] for key, values in submitted.items()}
        pointwise = [f.result() for f in pointwiseSubmitted]
    return results, pointwise


def plotProfiles(axes, results):
    """Plot the profiles at the final time on one of the grids with both weightings."""
    tEnd = 1.0 / (10.0 * PROFILE_VISCOSITY)
    # position, coordinates, velocity and pressure for both weightings
    profiles = {
        weight: profile(results[(PROFILE_VISCOSITY, weight)][CELLS.index(PROFILE_CELLS)].name)
        for weight in (CENTRAL, UPWIND)
    }
    x, y = profiles[CENTRAL][0], profiles[CENTRAL][1]
    yFine = np.linspace(0.0, 1.0, 400)
    exact = exactProfile(x, yFine, PROFILE_VISCOSITY, tEnd)
    title = (
        rf"$\nu = {PROFILE_VISCOSITY}$ m²/s, $t = {tEnd:g}$ s, $x = {x:g}$ m, "
        f"{PROFILE_CELLS} × {PROFILE_CELLS} cells"
    )
    for ax, component, label in ((axes[0], 0, "$u_x$ [m/s]"), (axes[1], 1, "$p$ [Pa]")):
        ax.plot(yFine, exact[component], color="black", linewidth=1.5, label="analytical")
        for weight, marker in ((CENTRAL, "o"), (UPWIND, "s")):
            ax.plot(
                y,
                profiles[weight][2 + component],
                marker,
                color=COLORS[weight],
                markerfacecolor="none",
                markersize=6,
                label=f"DuMux, {NAMES[weight]}",
            )
        ax.set_xlabel("$y$ [m]")
        ax.set_ylabel(label)
        ax.set_title(title, fontsize=10)
        ax.legend(fontsize=9)


def plotErrors(axes, results):
    """Plot the errors for both viscosities on the same scale."""
    h = np.array([1.0 / n for n in CELLS])
    for ax, nu in zip(axes, VISCOSITIES):
        for weight in (CENTRAL, UPWIND):
            for style, error, label in (
                ("-o", [r.velocityError for r in results[(nu, weight)]], "erru"),
                ("--s", [r.pressureError for r in results[(nu, weight)]], "errp"),
            ):
                ax.loglog(
                    h,
                    error,
                    style,
                    color=COLORS[weight],
                    markersize=5,
                    linewidth=1.5,
                    label=f"{label}, {NAMES[weight]}",
                )
        secondOrder = results[(nu, CENTRAL)][-1].velocityError * (h / h[-1]) ** 2
        firstOrder = results[(nu, UPWIND)][-1].velocityError * (h / h[-1])
        ax.loglog(h, secondOrder, color="black", linewidth=0.8, label=r"$\propto h^2$")
        ax.loglog(h, firstOrder, color="gray", linewidth=0.8, label=r"$\propto h$")
        ax.set_xlabel("$h$ [m]")
        ax.set_ylabel("relative space-time $L^2$ error")
        tEnd = 1.0 / (10.0 * nu)
        ax.set_title(rf"$\nu = {nu}$ m²/s, $T = {tEnd:g}$ s", fontsize=10)
        ax.legend(fontsize=8)
    axes[1].sharey(axes[0])


def main():
    """Run the study, print the tables and write the figure."""
    if len(sys.argv) < 2:
        sys.stderr.write(__doc__)
        sys.exit(1)

    executable = str(Path(sys.argv[1]).resolve())
    outdir = Path(sys.argv[2]).resolve() if len(sys.argv) > 2 else Path.cwd()

    results, pointwise = runStudy(executable)
    for (nu, weight), values in results.items():
        printTable(f"nu = {nu}, {NAMES[weight]}", CELLS, values)
    printTable(
        f"nu = {PROFILE_VISCOSITY}, central, pointwise initial and boundary velocity",
        POINTWISE_CELLS,
        pointwise,
    )
    printComplexityTable(CELLS, results[(PROFILE_VISCOSITY, CENTRAL)])

    fig, axes = plt.subplots(2, 2, figsize=(11, 8))
    plotProfiles(axes[0], results)
    plotErrors(axes[1], results)
    for ax in axes.flat:
        ax.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(outdir / "angeli_comparison.png", dpi=150)

    print(f"\nWrote {outdir / 'angeli_comparison.png'}")


if __name__ == "__main__":
    main()
