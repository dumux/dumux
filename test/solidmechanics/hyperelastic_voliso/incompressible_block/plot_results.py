#!/usr/bin/env python3
# SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later
"""Run the incompressible block on a finer hexahedral mesh and render the deformed block.

Usage:
  python3 plot_results.py <executable> [output directory]

<executable> is test_hyperelastic_incompressible_block_q2_q1. It runs with params.input of the
current directory on 8 x 8 x 8 hexahedra and writes incompressible_block_deformation.png to the
output directory (default: current directory): the deformed quarter that is computed and the full
block obtained by mirroring it at the symmetry planes, coloured by the displacement magnitude.
Requires numpy, matplotlib and pyvista.
"""

import subprocess
import sys
from pathlib import Path

CELLS = "8 8 8"
NAME = "incompressible_block_plot"
P = (50.0, 50.0, 50.0)


def run(executable):
    subprocess.run([str(executable), "params.input", "-Problem.Name", NAME, "-Grid.Cells", CELLS],
                   check=True, stdout=subprocess.DEVNULL)


def readFinalState():
    import numpy as np
    import pyvista as pv
    reader = pv.get_reader(f"{NAME}.pvd")
    reader.set_active_time_value(reader.time_values[-1])
    mesh = reader.read()[0]
    mesh.point_data["|u|"] = np.linalg.norm(np.asarray(mesh.point_data["d"]), axis=1)
    return mesh


def mirrored(quarter):
    """The full block from the deformed quarter, mirrored at the planes x = 50 mm and y = 50 mm."""
    copyX = quarter.reflect((1.0, 0.0, 0.0), point=(P[0], 0.0, 0.0))
    copyY = quarter.reflect((0.0, 1.0, 0.0), point=(0.0, P[1], 0.0))
    copyXY = copyX.reflect((0.0, 1.0, 0.0), point=(0.0, P[1], 0.0))
    return quarter.merge([copyX, copyY, copyXY])


def render(meshes, cameraPosition, label=None):
    """Off-screen render of the given (mesh, style) pairs, cropped to the drawn content."""
    import numpy as np
    import pyvista as pv
    plotter = pv.Plotter(off_screen=True, window_size=(900, 760))
    plotter.set_background("white")
    for mesh, style in meshes:
        plotter.add_mesh(mesh, **style)
    if label is not None:
        plotter.add_point_labels([label[0]], [label[1]], font_size=26, text_color="black", shape=None,
                                 show_points=False, always_visible=True, italic=True)
    plotter.camera_position = cameraPosition
    image = plotter.screenshot(return_img=True)
    plotter.close()
    rows, cols = np.where(image.min(axis=2) < 250)
    margin = 10
    return image[max(rows.min() - margin, 0):rows.max() + margin, max(cols.min() - margin, 0):cols.max() + margin]


def main():
    if len(sys.argv) < 2:
        sys.stderr.write(__doc__)
        sys.exit(1)

    import numpy as np
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.colors import ListedColormap, Normalize
    import pyvista as pv

    executable = Path(sys.argv[1]).resolve()
    outdir = Path(sys.argv[2]).resolve() if len(sys.argv) > 2 else Path.cwd()

    run(executable)
    mesh = readFinalState()
    vertexP = np.argmin(np.linalg.norm(mesh.points - np.array(P), axis=1))
    uzP = float(np.asarray(mesh.point_data["d"])[vertexP, 2])

    quarter = mesh.warp_by_vector("d")
    block = mirrored(quarter)
    deformedP = quarter.points[vertexP]

    cmap = ListedColormap(plt.get_cmap("Blues")(np.linspace(0.15, 1.0, 256)))
    norm = Normalize(0.0, float(mesh.point_data["|u|"].max()))
    meshStyle = dict(scalars="|u|", cmap=cmap, clim=[norm.vmin, norm.vmax], show_scalar_bar=False,
                     show_edges=True, edge_color="white", line_width=0.6,
                     ambient=0.35, diffuse=0.7, specular=0.0)
    outlineStyle = dict(color="#52514e", line_width=1.5)

    quarterImage = render(
        [(quarter, meshStyle),
         (pv.Box(bounds=(0.0, 50.0, 0.0, 50.0, 0.0, 50.0)).outline(), outlineStyle),
         (pv.Sphere(radius=0.8, center=deformedP), dict(color="black"))],
        [(185.0, 175.0, 150.0), (25.0, 25.0, 22.0), (0.0, 0.0, 1.0)],
        label=(deformedP + np.array([0.0, 0.0, 7.0]), "P"))
    blockImage = render(
        [(block, meshStyle),
         (pv.Box(bounds=(0.0, 100.0, 0.0, 100.0, 0.0, 50.0)).outline(), outlineStyle)],
        [(255.0, 225.0, 190.0), (50.0, 50.0, 18.0), (0.0, 0.0, 1.0)])

    fig, axes = plt.subplots(1, 2, figsize=(10, 4.2),
                             gridspec_kw=dict(width_ratios=[quarterImage.shape[1]/quarterImage.shape[0],
                                                            blockImage.shape[1]/blockImage.shape[0]]))
    for ax, image, title in zip(axes, [quarterImage, blockImage], ["modelled quarter", "full block, mirrored"]):
        ax.imshow(image)
        ax.set_title(title, fontsize=11)
        ax.axis("off")
    fig.colorbar(plt.cm.ScalarMappable(norm=norm, cmap=cmap), ax=axes, shrink=0.8,
                 label="$|\\mathbf{u}|$ [mm]")
    fig.savefig(outdir/"incompressible_block_deformation.png", dpi=120, bbox_inches="tight")

    print(f"|u_z(P)| = {abs(uzP):.2f} mm on {CELLS.replace(' ', ' x ')} hexahedra")
    print(f"Wrote {outdir/'incompressible_block_deformation.png'}")


if __name__ == "__main__":
    main()
