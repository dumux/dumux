#!/usr/bin/env python3
# SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later
"""Run the 2D hyperelastic Cook's membrane on a sequence of refined meshes and plot the results.

Usage:
  python3 plot_results.py <mini executable> <taylorhood executable> [output directory]

Runs both executables with params_mixed.input of the current directory on refinement levels 1-8
(MINI) and 1-7 (Taylor-Hood) and writes to the output directory (default: current directory)
- cooks_membrane_hyperelastic_convergence.png: relative error of the tip displacement u_y(A) at
  Poisson's ratio 0.35 with respect to the p-FEM value of Schröder et al. (2021) versus the number
  of degrees of freedom,
- cooks_membrane_hyperelastic_deformation.png: deformed membrane on Taylor-Hood refinement level 4.
The tip displacements of all runs and their extrapolation from every second refinement level are
printed.
Requires numpy, matplotlib and pyvista.
"""

import configparser
import subprocess
import sys
from pathlib import Path

LEVELS = {"MINI": [1, 2, 3, 4, 5, 6, 7, 8], "Taylor-Hood": [1, 2, 3, 4, 5, 6, 7]}
DEFORMATION_LEVEL = 4
TIP = (48.0, 60.0)
# Table 3 of Schröder et al. (2021), plane strain, Poisson's ratio 0.35
P_FEM = 10.5953
T2 = {"dofs": 115680, "uy": 10.5919}
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


def extrapolated(levels, lam):
    """Richardson extrapolation of u_y(A) from the last three even and the last three odd levels.

    Every refinement level bisects each triangle once, so the mesh size halves every second level.
    Returns the extrapolated values and the ratios of successive changes over two levels.
    """
    estimates, ratios = [], []
    for parity in (0, 1):
        uy = [r["uy"][lam] for r in levels if r["level"] % 2 == parity][-3:]
        ratios.append((uy[1] - uy[0])/(uy[2] - uy[1]))
        estimates.append(uy[2] + (uy[2] - uy[1])/(ratios[-1] - 1.0))
    return estimates, ratios


def plotConvergence(results, lambdas, outdir):
    import numpy as np
    import matplotlib.pyplot as plt

    lam = lambdas[0]
    limit = np.mean([value for levels in results.values() for value in extrapolated(levels, lam)[0]])
    fig, ax = plt.subplots(figsize=(6.4, 4.4))
    for element, levels in results.items():
        dofs = [r["dofs"] for r in levels]
        error = [abs(r["uy"][lam] - P_FEM)/P_FEM for r in levels]
        ax.loglog(dofs, error, marker=MARKERS[element], color=COLORS[element], linewidth=2,
                  markersize=8, markeredgecolor="white", markeredgewidth=1.5, label=element)
    ax.loglog(T2["dofs"], abs(T2["uy"] - P_FEM)/P_FEM, "D", color="#52514e", markersize=8,
              label="T2, Schröder et al. (2021)")
    ax.axhline(abs(limit - P_FEM)/P_FEM, color="#52514e", linestyle=":", linewidth=1.2,
               label="extrapolated, MINI and Taylor-Hood")
    guide = np.array([3e3, 8e4])
    ax.loglog(guide, 3e-4*(guide/1e4)**-0.75, color="#52514e", linestyle="--", linewidth=1.0)
    ax.annotate("$\\propto N^{-3/4}$", (guide[1], 3e-4*(guide[1]/1e4)**-0.75), xytext=(4, 0),
                textcoords="offset points", fontsize=10, va="center")
    ax.set_xlabel("degrees of freedom $N$ (displacement and pressure)")
    ax.set_ylabel("relative error of $u_y(A)$")
    ax.set_title("$\\nu = 0.35$, reference: p-FEM, Schröder et al. (2021)", fontsize=11)
    ax.grid(alpha=0.3)
    ax.legend(fontsize=9)
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

    results = {}
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
            if element == "Taylor-Hood" and level == DEFORMATION_LEVEL:
                deformationMeshes = meshes

    plotConvergence(results, lambdas, outdir)
    plotDeformation(deformationMeshes, lambdas, mu, outdir)

    print("\n| element | refinement | dofs | " + " | ".join(f"lambda = {lam:g}" for lam in lambdas) + " |")
    for element, levels in results.items():
        for r in levels:
            print(f"| {element} | {r['level']} | {r['dofs']} | " + " | ".join(f"{r['uy'][lam]:.4f}" for lam in lambdas) + " |")
    print()
    for element, levels in results.items():
        for lam in lambdas:
            (even, odd), (evenRatio, oddRatio) = extrapolated(levels, lam)
            print(f"{element}, lambda = {lam:g}: extrapolated {even:.4f} (even levels), {odd:.4f} (odd levels), "
                  f"ratio of changes over two levels {evenRatio:.2f} (even), {oddRatio:.2f} (odd)")
    print(f"\nWrote {outdir/'cooks_membrane_hyperelastic_convergence.png'} and {outdir/'cooks_membrane_hyperelastic_deformation.png'}")


if __name__ == "__main__":
    main()
