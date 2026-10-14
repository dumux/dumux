#!/usr/bin/env python3
# SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later
"""Run the 2D hyperelastic Cook's membrane on a sequence of refined meshes and plot the results.

Usage:
  python3 plot_results.py <mini executable> <taylorhood executable> [output directory]

Runs both executables with params_mixed.input of the current directory on refinement levels 1-5
(MINI) and 1-4 (Taylor-Hood) and writes to the output directory (default: current directory)
- cooks_membrane_hyperelastic_convergence.png: tip displacement u_y(A) versus the number of
  degrees of freedom for Poisson's ratios 0.35 and 0.4999,
- cooks_membrane_hyperelastic_deformation.png: deformed membrane on the finest Taylor-Hood mesh.
The tip displacements of all runs are printed as a table.
Requires numpy, matplotlib and pyvista.
"""

import configparser
import subprocess
import sys
from pathlib import Path

LEVELS = {"MINI": [1, 2, 3, 4, 5], "Taylor-Hood": [1, 2, 3, 4]}
TIP = (48.0, 60.0)
# Table 3 of Schröder et al. (2021), plane strain, Poisson's ratio 0.35
REFERENCES = {"p-FEM": 10.5953, "T2": 10.5919}
COLORS = {"MINI": "#2a78d6", "Taylor-Hood": "#eb6834"}
MARKERS = {"MINI": "s", "Taylor-Hood": "o"}


def readParameters(fileName):
    parser = configparser.ConfigParser(inline_comment_prefixes=("#",))
    parser.optionxform = str
    parser.read(fileName)
    lambdas = [float(v) for v in parser["SpatialParams"]["LambdaValues"].split()]
    return lambdas, float(parser["SpatialParams"]["Mu"])


def outputFile(executable, lam):
    return f"{executable.name}_lambda{int(lam)}-00000.vtu"


def numDofs(mesh, element):
    """Displacement and pressure degrees of freedom of the mixed scheme on a triangle mesh."""
    import numpy as np
    triangles = mesh.cells.reshape(-1, 4)[:, 1:]
    numVertices, numElements = mesh.n_points, triangles.shape[0]
    if element == "MINI":
        return 2*(numVertices + numElements) + numVertices
    edges = np.sort(np.concatenate([triangles[:, [0, 1]], triangles[:, [1, 2]], triangles[:, [2, 0]]]), axis=1)
    numEdges = len(np.unique(edges, axis=0))
    return 2*(numVertices + numEdges) + numVertices


def tipDisplacement(mesh):
    import numpy as np
    distance = np.linalg.norm(mesh.points[:, :2] - np.array(TIP), axis=1)
    vertex = np.argmin(distance)
    if distance[vertex] > 1e-8:
        raise RuntimeError("The tip A is not a vertex of the mesh")
    return float(np.asarray(mesh.point_data["d"])[vertex, 1])


def run(executable, level, lambdas):
    """Run one refinement level and return the output mesh per lambda."""
    import pyvista as pv
    subprocess.run([str(executable), "params_mixed.input", "-Grid.Refinement", str(level)],
                   check=True, stdout=subprocess.DEVNULL)
    return {lam: pv.read(outputFile(executable, lam)) for lam in lambdas}


def plotConvergence(results, lambdas, mu, outdir):
    import matplotlib.pyplot as plt
    from matplotlib.ticker import FixedLocator, FuncFormatter, NullFormatter

    shown = [lambdas[0], lambdas[-1]]
    fig, axes = plt.subplots(1, 2, figsize=(10, 4.2))
    for ax, lam in zip(axes, shown):
        for element, levels in results.items():
            dofs = [r["dofs"] for r in levels]
            uy = [r["uy"][lam] for r in levels]
            ax.semilogx(dofs, uy, marker=MARKERS[element], color=COLORS[element], linewidth=2,
                        markersize=8, markeredgecolor="white", markeredgewidth=1.5, label=element)
        ax.set_title(f"$\\nu = {lam/(2.0*(lam + mu)):.4g}$, $\\lambda = {lam:g}$ MPa", fontsize=11)
        ax.set_xlabel("degrees of freedom (displacement and pressure)")
        ax.set_ylabel("$u_y(A)$ [mm]")
        ax.xaxis.set_major_locator(FixedLocator([2e3, 5e3, 1e4, 2e4, 4e4]))
        ax.xaxis.set_major_formatter(FuncFormatter(lambda value, _: f"{value:,.0f}"))
        ax.xaxis.set_minor_formatter(NullFormatter())
        ax.grid(alpha=0.3)

    for (name, value), style in zip(REFERENCES.items(), [":", "-."]):
        axes[0].axhline(value, color="#52514e", linestyle=style, linewidth=1.2,
                        label=f"{name}, Schröder et al. (2021)")
    axes[0].set_ylim(top=max(REFERENCES.values()) + 0.01)
    axes[0].legend(fontsize=9)
    fig.tight_layout()
    fig.savefig(outdir/"cooks_membrane_hyperelastic_convergence.png", dpi=150)


def plotDeformation(meshes, lambdas, mu, outdir):
    import numpy as np
    import matplotlib.pyplot as plt
    import matplotlib.tri as tri
    from matplotlib.colors import ListedColormap

    shown = [lambdas[0], lambdas[-1]]
    blues = plt.get_cmap("Blues")
    cmap = ListedColormap(blues(np.linspace(0.15, 1.0, 256)))
    vmax = max(np.linalg.norm(np.asarray(meshes[lam].point_data["d"]), axis=1).max() for lam in shown)
    outline = np.array([[0.0, 0.0], [48.0, 44.0], [48.0, 60.0], [0.0, 44.0], [0.0, 0.0]])

    fig, axes = plt.subplots(1, 2, figsize=(9, 4.6), sharey=True)
    for ax, lam in zip(axes, shown):
        mesh = meshes[lam]
        d = np.asarray(mesh.point_data["d"])[:, :2]
        x = mesh.points[:, :2] + d
        triangles = mesh.cells.reshape(-1, 4)[:, 1:]
        triangulation = tri.Triangulation(x[:, 0], x[:, 1], triangles)
        image = ax.tripcolor(triangulation, np.linalg.norm(d, axis=1), shading="gouraud",
                             cmap=cmap, vmin=0.0, vmax=vmax)
        ax.triplot(triangulation, color="white", linewidth=0.15, alpha=0.5)
        ax.plot(outline[:, 0], outline[:, 1], color="#52514e", linestyle="--", linewidth=1.0)
        tip = np.argmin(np.linalg.norm(mesh.points[:, :2] - np.array(TIP), axis=1))
        ax.plot(*x[tip], "o", color="black", markersize=5)
        ax.annotate("$A$", x[tip], xytext=(6, 2), textcoords="offset points", fontsize=12)
        ax.set_title(f"$\\nu = {lam/(2.0*(lam + mu)):.4g}$", fontsize=11)
        ax.set_aspect("equal")
        ax.set_xlabel("$x$ [mm]")
    axes[0].set_ylabel("$y$ [mm]")
    fig.colorbar(image, ax=axes, label="$|\\mathbf{u}|$ [mm]", shrink=0.85)
    fig.savefig(outdir/"cooks_membrane_hyperelastic_deformation.png", dpi=120, bbox_inches="tight")


def main():
    if len(sys.argv) < 3:
        sys.stderr.write(__doc__)
        sys.exit(1)

    import matplotlib
    matplotlib.use("Agg")

    executables = {"MINI": Path(sys.argv[1]).resolve(), "Taylor-Hood": Path(sys.argv[2]).resolve()}
    outdir = Path(sys.argv[3]).resolve() if len(sys.argv) > 3 else Path.cwd()
    lambdas, mu = readParameters("params_mixed.input")

    results, finestMeshes = {}, {}
    for element, executable in executables.items():
        results[element] = []
        for level in LEVELS[element]:
            meshes = run(executable, level, lambdas)
            results[element].append({
                "level": level,
                "dofs": numDofs(meshes[lambdas[0]], element),
                "uy": {lam: tipDisplacement(mesh) for lam, mesh in meshes.items()},
            })
            print(f"{element} level {level}: {results[element][-1]['dofs']} dofs", flush=True)
        finestMeshes[element] = meshes

    plotConvergence(results, lambdas, mu, outdir)
    plotDeformation(finestMeshes["Taylor-Hood"], lambdas, mu, outdir)

    print("\n| element | refinement | dofs | " + " | ".join(f"lambda = {lam:g}" for lam in lambdas) + " |")
    for element, levels in results.items():
        for r in levels:
            print(f"| {element} | {r['level']} | {r['dofs']} | " + " | ".join(f"{r['uy'][lam]:.3f}" for lam in lambdas) + " |")
    print(f"\nWrote {outdir/'cooks_membrane_hyperelastic_convergence.png'} and {outdir/'cooks_membrane_hyperelastic_deformation.png'}")


if __name__ == "__main__":
    main()
