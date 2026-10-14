#!/usr/bin/env python3
# SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later
"""
Render the animation of the decaying vortex, see README.md, with the kinetic energy and the
enstrophy over time (bottom) and on top

  * 2D (run of 'benchmark.py --study animation'): the velocity magnitude with arrows, the pressure
    and the vorticity magnitude, plotted with matplotlib,
  * 3D (run of the energy study 'benchmark.py --study energy --dim 3'): isosurfaces of the velocity
    magnitude and the vorticity magnitude, rendered with pyvista.

The fields are read from the higher-order ParaView output with pyvista and sampled on a regular grid
(which evaluates the polynomials of the higher-order cells). Color ranges and isovalues are fixed
over all time steps, such that the decay of the vortex is visible. The video is written with ffmpeg
to results/taylorgreen_animation_<dim>d_<scheme>.mp4:

    python3 animate.py [--dim 2|3] [--manifest ...] [--scheme ...] [--output-dir results]
"""

import argparse
import json
import sys
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
import pyvista as pv
from matplotlib.animation import FFMpegWriter

from visualize import ENERGY_LABEL, ENSTROPHY_LABEL, Context, readRows


def sample(mesh, resolution):
    """
    Sample the fields of the (higher-order) mesh on a regular grid with the given number of
    points per direction, which evaluates the polynomials of the higher-order cells
    """
    lower, upper = np.array(mesh.bounds[::2]), np.array(mesh.bounds[1::2])
    dimensions = [resolution if high > low else 1 for low, high in zip(lower, upper)]
    spacing = [(high - low) / max(n - 1, 1) for low, high, n in zip(lower, upper, dimensions)]
    return pv.ImageData(dimensions=dimensions, spacing=spacing, origin=lower).sample(mesh)


def vorticityMagnitude(fields):
    """The magnitude of the vorticity, a scalar in 2D and a vector in 3D"""
    vorticity = np.asarray(fields.point_data["vorticity"])
    return np.abs(vorticity) if vorticity.ndim == 1 else np.linalg.norm(vorticity, axis=1)


# the plotted fields in 2D: key, color map, colorbar label and title
FIELDS = (
    ("speed", "viridis", r"$\|\mathbf{u}_h\|$ [m/s]", "velocity magnitude"),
    ("pressure", "coolwarm", r"$p_h$ [Pa]", "pressure"),
    ("vorticity", "inferno", r"$|\omega_h|$ [1/s]", "vorticity magnitude"),
)


@dataclass
class Frames:
    """The sampled fields of all time steps on a regular grid"""

    times: list
    bounds: tuple = ()
    velocity: list = field(default_factory=list)
    fields: dict = field(default_factory=lambda: {key: [] for key, *_ in FIELDS})

    def colorRange(self, key):
        """Color range of a field fixed over all time steps (symmetric for the pressure)"""
        bound = np.abs(self.fields[key]).max()
        return (-bound, bound) if key == "pressure" else (0.0, bound)


def sampleFrames(pvd, resolution):
    """Sample velocity, pressure and vorticity magnitude of all time steps of a .pvd series"""
    reader = pv.get_reader(pvd)
    frames = Frames(times=list(reader.time_values))
    n = resolution
    for time in frames.times:
        reader.set_active_time_value(time)
        sampled = sample(reader.read()[0], n)
        u = np.asarray(sampled.point_data["velocity"])[:, :2]
        # the higher-order output writes scalar fields with three components
        p = np.asarray(sampled.point_data["pressure"]).reshape(n * n, -1)[:, 0]
        frames.velocity.append(u.reshape(n, n, 2))
        frames.fields["speed"].append(np.linalg.norm(u, axis=1).reshape(n, n))
        frames.fields["pressure"].append(p.reshape(n, n))
        frames.fields["vorticity"].append(vorticityMagnitude(sampled).reshape(n, n))
        frames.bounds = sampled.bounds[:4]
    return frames


def drawFields(fig, axes, frames):
    """Draw the fields of the first time step, returns the images"""
    images = {}
    for ax, (key, cmap, label, title) in zip(axes, FIELDS):
        vmin, vmax = frames.colorRange(key)
        images[key] = ax.imshow(
            frames.fields[key][0],
            origin="lower",
            extent=frames.bounds,
            cmap=cmap,
            vmin=vmin,
            vmax=vmax,
        )
        fig.colorbar(images[key], ax=ax, label=label, shrink=0.8)
        ax.set_xlabel(r"$x$ [m]")
        ax.set_ylabel(r"$y$ [m]")
        ax.set_title(title)
    return images


def drawArrows(ax, frames, numArrows):
    """
    Arrows at a regular subset of the sampling points, scaled with the maximum velocity over all
    time steps such that the decay is visible; returns the arrows and the subset
    """
    n = len(frames.velocity[0])
    step = max(n // numArrows, 1)
    subset = (slice(step // 2, None, step), slice(step // 2, None, step))
    x = np.linspace(frames.bounds[0], frames.bounds[1], n)[subset[0]]
    y = np.linspace(frames.bounds[2], frames.bounds[3], n)[subset[1]]
    arrowLength = 0.9 * (frames.bounds[1] - frames.bounds[0]) / numArrows
    velocity = frames.velocity[0][subset]
    arrows = ax.quiver(
        *np.meshgrid(x, y),
        velocity[..., 0],
        velocity[..., 1],
        angles="xy",
        scale_units="xy",
        scale=frames.colorRange("speed")[1] / arrowLength,
        color="black",
    )
    return arrows, subset


def drawIndicators(context, axes, rows, label):
    """
    Kinetic energy and enstrophy (normalized with the analytical values at t = 0) with the
    analytical curves; returns the (initially empty) numerical curves and their values
    """
    t = [row["t"] for row in rows]
    curves = []
    for ax, key, title, symbol in (
        (axes[0], "kineticEnergy", f"kinetic energy {ENERGY_LABEL}", "E"),
        (axes[1], "enstrophy", f"enstrophy {ENSTROPHY_LABEL}", r"\boldsymbol{\varepsilon}"),
    ):
        initial = rows[0][f"{key}Exact"]
        numerical = [row[key] / initial for row in rows]
        ax.plot(t, [row[f"{key}Exact"] / initial for row in rows], "k--", label="exact")
        (line,) = ax.plot([], [], color="C0", linewidth=2.0, label=label)
        (marker,) = ax.plot([], [], "o", color="C0")
        curves.append((line, marker, numerical))
        ax.set_xlim(t[0], t[-1])
        ax.set_ylim(0.0, 1.05 * max(numerical))
        ax.set_xlabel(r"$t$ [s]")
        ax.set_ylabel(context.tex(rf"${symbol}_h(t) / {symbol}(0)$"))
        ax.set_title(context.tex(title), fontsize=12)
        ax.legend(loc="upper right")
        ax.grid(True, alpha=0.3)
    return t, curves


def parseArguments():
    """Parse the command line arguments"""
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("--dim", type=int, choices=(2, 3), default=2)
    parser.add_argument(
        "--manifest",
        help="manifest written by benchmark.py (default: runs/taylorgreen_animation_2d.json "
        "in 2D, runs/taylorgreen_energy_3d.json in 3D)",
    )
    parser.add_argument("--scheme", help="momentum scheme of the run (default: the first one)")
    parser.add_argument("--output-dir", default="results", help="directory of the video")
    parser.add_argument("--fps", type=int, help="frames per second (default: 10 in 2D, 4 in 3D)")
    parser.add_argument(
        "--resolution",
        type=int,
        help="sampling points per direction (default: 200 in 2D, 41 in 3D)",
    )
    parser.add_argument(
        "--arrows", type=int, default=16, help="number of arrows per direction (2D)"
    )
    args = parser.parse_args()
    study = "animation" if args.dim == 2 else "energy"
    args.manifest = args.manifest or f"runs/taylorgreen_{study}_{args.dim}d.json"
    args.fps = args.fps or (10 if args.dim == 2 else 4)
    args.resolution = args.resolution or (200 if args.dim == 2 else 41)
    return args


def drawFigure(context, frames, rows, label, numArrows):
    """
    The figure with the fields (top) and the indicators (bottom) of the first time step and a
    function that updates it to time step i
    """
    fig = context.plt.figure(figsize=(16, 9.6), constrained_layout=True)
    layout = fig.add_gridspec(2, 1, height_ratios=(1.25, 1))
    fieldAxes = layout[0].subgridspec(1, 3).subplots()
    images = drawFields(fig, fieldAxes, frames)
    arrows, subset = drawArrows(fieldAxes[0], frames, numArrows)
    t, curves = drawIndicators(context, layout[1].subgridspec(1, 2).subplots(), rows, label)

    def update(i):
        for key, image in images.items():
            image.set_data(frames.fields[key][i])
        velocity = frames.velocity[i][subset]
        arrows.set_UVC(velocity[..., 0], velocity[..., 1])
        for line, marker, values in curves:
            line.set_data(t[: i + 1], values[: i + 1])
            marker.set_data([t[i]], [values[i]])
        fig.suptitle(f"2D Taylor-Green vortex, {label}, $t = {frames.times[i]:.2f}$ s")

    return fig, update


# the plotted fields in 3D: key, color map, colorbar label and title
ISOSURFACES = (
    ("speed", "viridis", r"$\|\mathbf{u}_h\|$ [m/s]", "velocity magnitude"),
    ("vorticity", "plasma", r"$\|\boldsymbol{\omega}_h\|$ [1/s]", "vorticity magnitude"),
)
# isovalues as fractions of the maximum of the field at t = 0
ISO_FRACTIONS = (0.25, 0.5, 0.75)


@dataclass
class Frames3d:
    """The sampled velocity and vorticity magnitude of all time steps on a regular grid"""

    times: list
    grids: list = field(default_factory=list)

    def isovalues(self, key):
        """Isovalues fixed over all time steps (fractions of the maximum at t = 0)"""
        return [f * self.grids[0].point_data[key].max() for f in ISO_FRACTIONS]


def sampleFrames3d(pvd, resolution):
    """Sample velocity and vorticity magnitude of all time steps of a .pvd series"""
    reader = pv.get_reader(pvd)
    frames = Frames3d(times=list(reader.time_values))
    for time in frames.times:
        reader.set_active_time_value(time)
        sampled = sample(reader.read()[0], resolution)
        grid = pv.ImageData(
            dimensions=sampled.dimensions, spacing=sampled.spacing, origin=sampled.origin
        )
        grid.point_data["speed"] = np.linalg.norm(sampled.point_data["velocity"], axis=1)
        grid.point_data["vorticity"] = vorticityMagnitude(sampled)
        frames.grids.append(grid)
    return frames


class IsosurfaceRenderer:
    """Renders the isosurfaces of a field of a 3D regular grid with pyvista (off screen)"""

    def __init__(self, key, cmap, isovalues):
        self.key, self.cmap, self.isovalues = key, cmap, isovalues
        self.plotter = pv.Plotter(off_screen=True, window_size=(1000, 900))
        self.plotter.set_background("white")

    def render(self, grid):
        """The image of the isosurfaces of the given grid"""
        self.plotter.clear_actors()
        isosurfaces = grid.contour(isosurfaces=self.isovalues, scalars=self.key)
        if isosurfaces.n_points > 0:
            self.plotter.add_mesh(
                isosurfaces,
                scalars=self.key,
                cmap=self.cmap,
                clim=(self.isovalues[0], self.isovalues[-1]),
                opacity=0.6,
                smooth_shading=True,
                show_scalar_bar=False,
            )
        self.plotter.add_mesh(grid.outline(), color="black")
        lower, upper = np.array(grid.bounds[::2]), np.array(grid.bounds[1::2])
        center, length = 0.5 * (lower + upper), max(upper - lower)
        self.plotter.camera_position = [
            center + length * np.array([-2.2, -2.6, 1.9]),
            center,
            (0.0, 0.0, 1.0),
        ]
        return self.plotter.screenshot(return_img=True)


def drawIsosurfaces(context, fig, axes, frames):
    """Image of the isosurfaces of the first time step with a colorbar of the isovalues per field"""
    images = []
    for ax, (key, cmap, colorbarLabel, title) in zip(axes, ISOSURFACES):
        isovalues = frames.isovalues(key)
        renderer = IsosurfaceRenderer(key, cmap, isovalues)
        images.append((ax.imshow(renderer.render(frames.grids[0])), renderer))
        ax.set_axis_off()
        ax.set_title(f"{title}, isovalues {', '.join(map(str, ISO_FRACTIONS))} of max at $t = 0$")
        colorbar = fig.colorbar(
            context.plt.cm.ScalarMappable(context.plt.Normalize(isovalues[0], isovalues[-1]), cmap),
            ax=ax,
            label=context.tex(colorbarLabel),
            shrink=0.7,
        )
        colorbar.set_ticks(isovalues)
    return images


def drawFigure3d(context, frames, rows, label):
    """
    The figure with the isosurfaces (top) and the indicators (bottom) of the first time step and
    a function that updates it to time step i
    """
    fig = context.plt.figure(figsize=(16, 11), constrained_layout=True)
    layout = fig.add_gridspec(2, 1, height_ratios=(1.6, 1))
    images = drawIsosurfaces(context, fig, layout[0].subgridspec(1, 2).subplots(), frames)
    t, curves = drawIndicators(context, layout[1].subgridspec(1, 2).subplots(), rows, label)

    def update(i):
        for image, renderer in images:
            image.set_data(renderer.render(frames.grids[i]))
        for line, marker, values in curves:
            line.set_data(t[: i + 1], values[: i + 1])
            marker.set_data([t[i]], [values[i]])
        fig.suptitle(f"3D Beltrami (ABC) flow, {label}, $t = {frames.times[i]:.2f}$ s")

    return fig, update


def selectGroup(manifest, scheme):
    """The group of runs of the given momentum scheme (default: the first one)"""
    groups = [g for g in manifest["groups"] if scheme is None or g["scheme"] == scheme]
    if not groups:
        sys.exit(f"No run of the scheme {scheme} in the manifest")
    return groups[0]


def main():
    """Sample the fields of all time steps and write one video frame per time step"""
    args = parseArguments()
    with open(args.manifest) as manifestFile:
        group = selectGroup(json.load(manifestFile), args.scheme)
    run = group["runs"][0]
    if not Path(run["pvd"]).is_file():
        sys.exit(f"{run['pvd']} not found (higher-order output requires GridFormat)")
    rows = readRows(run["errors"])

    context = Context(args.output_dir)
    if args.dim == 2:
        frames = sampleFrames(run["pvd"], args.resolution)
        fig, update = drawFigure(context, frames, rows, group["label"], args.arrows)
    else:
        frames = sampleFrames3d(run["pvd"], args.resolution)
        fig, update = drawFigure3d(context, frames, rows, group["label"])
    if len(frames.times) != len(rows):
        sys.exit(f"{len(frames.times)} time steps in {run['pvd']}, but {len(rows)} in the errors")

    context.outputDir.mkdir(parents=True, exist_ok=True)
    video = context.outputDir / f"taylorgreen_animation_{args.dim}d_{group['scheme']}.mp4"
    writer = FFMpegWriter(fps=args.fps, codec="libx264", extra_args=["-pix_fmt", "yuv420p"])
    with writer.saving(fig, video, dpi=100):
        for i in range(len(frames.times)):
            update(i)
            writer.grab_frame()
    print(f"Wrote {video}")


if __name__ == "__main__":
    main()
