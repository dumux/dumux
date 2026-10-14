#!/usr/bin/env python3
# SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later
"""
Render the animation of the decaying 2D Taylor-Green vortex simulated by
'benchmark.py --study animation', see README.md: the velocity magnitude with arrows, the pressure
and the vorticity magnitude (top) and the kinetic energy and the enstrophy over time (bottom).

The fields are read from the higher-order ParaView output with pyvista, sampled on a regular grid
(which evaluates the polynomials of the higher-order cells) and plotted with matplotlib with color
ranges fixed over all time steps, such that the decay of the vortex is visible. The video is written
with ffmpeg to results/taylorgreen_animation_2d_<scheme>.mp4:

    python3 animate.py [--manifest runs/taylorgreen_animation_2d.json] [--output-dir results]
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


# the plotted fields: key, color map, colorbar label and title
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
    parser.add_argument(
        "--manifest",
        default="runs/taylorgreen_animation_2d.json",
        help="manifest of the animation run written by benchmark.py",
    )
    parser.add_argument("--output-dir", default="results", help="directory of the video")
    parser.add_argument("--fps", type=int, default=10, help="frames per second")
    parser.add_argument("--resolution", type=int, default=200, help="sampling points per direction")
    parser.add_argument("--arrows", type=int, default=16, help="number of arrows per direction")
    return parser.parse_args()


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


def main():
    """Sample the fields of all time steps and write one video frame per time step"""
    args = parseArguments()
    with open(args.manifest) as manifestFile:
        group = json.load(manifestFile)["groups"][0]
    run = group["runs"][0]
    if not Path(run["pvd"]).is_file():
        sys.exit(f"{run['pvd']} not found (higher-order output requires GridFormat)")
    rows = readRows(run["errors"])
    frames = sampleFrames(run["pvd"], args.resolution)
    if len(frames.times) != len(rows):
        sys.exit(f"{len(frames.times)} time steps in {run['pvd']}, but {len(rows)} in the errors")

    context = Context(args.output_dir)
    fig, update = drawFigure(context, frames, rows, group["label"], args.arrows)
    context.outputDir.mkdir(parents=True, exist_ok=True)
    video = context.outputDir / f"{run['name']}.mp4"
    writer = FFMpegWriter(fps=args.fps, codec="libx264", extra_args=["-pix_fmt", "yuv420p"])
    with writer.saving(fig, video, dpi=100):
        for i in range(len(frames.times)):
            update(i)
            writer.grab_frame()
    print(f"Wrote {video}")


if __name__ == "__main__":
    main()
