#!/usr/bin/env python3
# SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

"""
DFG flow around a cylinder benchmark (Schäfer & Turek, 1996): video of a time-dependent run

Reads the VTK output of a run of benchmark.py (runs/2d<case>/<run>/vtk/<run>.pvd) and the time
series of the benchmark indicators (runs/2d<case>/<run>/<run>.json) and renders a video with

* the velocity magnitude |u|,
* the vorticity omega = dv/dx - du/dy, which shows the Karman vortex street,
* the pressure p,
* the drag and lift coefficients over time with a marker at the current time.

The fields are sampled on a uniform image of the channel (the velocity of the VTK output is the
element average, it is interpolated to the vertices first). The color ranges are fixed for the
whole video: they are the 1st and 99th percentiles over a subset of the frames (symmetric for the
vorticity, starting at zero for the velocity magnitude), so single extreme values near the
cylinder do not wash out the vortex street.

Default time window: case 2D-2 the last three shedding periods (about 1 s), case 2D-3 the whole run.
The video is written to results/video_2d<case>_<run>.mp4 (--format gif for a GIF).

Usage (in the build directory test/freeflow/navierstokes/flowaroundcylinder):
    python3 video.py --case 2 --run pq1bubble_box_r0_dirk2
    python3 video.py --case 3 --run pq1bubble_box_r0_dirk2 --stride 2 --format gif
"""

import argparse
import os
import sys
import xml.etree.ElementTree as ET

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FFMpegWriter, PillowWriter
import pyvista as pv

from benchmark import RESULTS_DIR, runsDir
from visualize import readRun

plt.rcParams.update({"mathtext.fontset": "cm", "font.family": "serif", "font.size": 11})

# channel [0, 2.2] x [0, 0.41] m and cylinder (center, radius) of the benchmark
CHANNEL = (2.2, 0.41)
CYLINDER = ((0.2, 0.2), 0.05)
# mean inflow velocity [m/s] and cylinder diameter [m], shedding period from Strouhal number 0.3
MEAN_VELOCITY, DIAMETER = 1.0, 0.1
PERIOD = DIAMETER / (0.3 * MEAN_VELOCITY)
# number of frames used to determine the color ranges
NUM_RANGE_FRAMES = 40

VELOCITY = "velocity_liq (m/s)"
PRESSURE = "p"
CASE_TITLES = {2: "2D-2 (periodic)", 3: "2D-3 (time-dependent inflow)"}
# fields of the video (key, color bar label, color map)
FIELD_PANELS = [
    ("speed", r"$|\mathbf{u}|$ [m/s]", "viridis"),
    ("vorticity", r"$\omega = \partial_x v - \partial_y u$ [1/s]", "RdBu_r"),
    ("pressure", r"$p$ [Pa]", "cividis"),
]


def readCollection(pvdFile):
    """The (time, file name) pairs of a ParaView collection file"""
    folder = os.path.dirname(pvdFile)
    root = ET.parse(pvdFile).getroot()
    return [
        (float(dataSet.get("timestep")), os.path.join(folder, dataSet.get("file")))
        for dataSet in root.iter("DataSet")
    ]


def samplingImage(resolution):
    """Uniform image of the channel with the given number of pixels in x direction"""
    spacing = CHANNEL[0] / (resolution - 1)
    numY = int(round(CHANNEL[1] / spacing)) + 1
    return pv.ImageData(dimensions=(resolution, numY, 1), spacing=(spacing, spacing, 1.0))


def sampleFields(fileName, image):
    """Velocity magnitude, vorticity and pressure of a VTK file on the image as 2D arrays (y, x)"""
    mesh = pv.read(fileName)
    if VELOCITY in mesh.cell_data:
        mesh = mesh.cell_data_to_point_data(pass_cell_data=False)
    sampled = image.sample(mesh)
    shape = (image.dimensions[1], image.dimensions[0])
    velocity = [sampled.point_data[VELOCITY][:, i].reshape(shape) for i in range(2)]
    spacing = image.spacing[0]
    dvdx = np.gradient(velocity[1], spacing, axis=1)
    dudy = np.gradient(velocity[0], spacing, axis=0)
    fields = {
        "speed": np.hypot(*velocity),
        "vorticity": dvdx - dudy,
        "pressure": sampled.point_data[PRESSURE].reshape(shape),
    }
    valid = sampled.point_data["vtkValidPointMask"].reshape(shape).astype(bool)
    return {key: np.where(valid & fluidMask(image), field, np.nan) for key, field in fields.items()}


def fluidMask(image):
    """Pixels of the image outside of the cylinder"""
    x, y = np.meshgrid(*(np.arange(n) * image.spacing[0] for n in image.dimensions[:2]))
    (centerX, centerY), radius = CYLINDER
    return (x - centerX) ** 2 + (y - centerY) ** 2 > radius**2


def colorRanges(frames, image):
    """Fixed color ranges from the 1st and 99th percentiles over a subset of the frames"""
    subset = frames[:: max(1, len(frames) // NUM_RANGE_FRAMES)]
    values = {"speed": [], "vorticity": [], "pressure": []}
    for _, fileName in subset:
        for key, field in sampleFields(fileName, image).items():
            values[key].append(field[np.isfinite(field)])
    values = {key: np.concatenate(arrays) for key, arrays in values.items()}
    omega = np.percentile(np.abs(values["vorticity"]), 99)
    return {
        "speed": (0.0, np.percentile(values["speed"], 99)),
        "vorticity": (-omega, omega),
        "pressure": tuple(np.percentile(values["pressure"], [1, 99])),
    }


def defaultWindow(case, frames):
    """Time window of the video: 2D-2 the last three shedding periods, 2D-3 the whole run"""
    tEnd = frames[-1][0]
    return (max(frames[0][0], tEnd - 3 * PERIOD), tEnd) if case == 2 else (frames[0][0], tEnd)


def parseArguments():
    """The command line arguments"""
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("--case", type=int, choices=[2, 3], required=True, help="benchmark case")
    parser.add_argument("--run", default="pq1bubble_box_r0_dirk2", help="name of the run")
    parser.add_argument("--tstart", type=float, default=None, help="start time of the video")
    parser.add_argument("--tend", type=float, default=None, help="end time of the video")
    parser.add_argument("--stride", type=int, default=1, help="use every n-th VTK file")
    parser.add_argument("--fps", type=int, default=25, help="frames per second")
    parser.add_argument("--resolution", type=int, default=1100, help="sampling points in x")
    parser.add_argument("--dpi", type=int, default=100, help="resolution of the frames")
    parser.add_argument("--format", choices=["mp4", "gif"], default="mp4", help="video format")
    parser.add_argument("--output", default=None, help="output file (default in results/)")
    return parser.parse_args()


def selectFrames(case, pvdFile, args):
    """The (time, file name) pairs of the video and its time window"""
    frames = readCollection(pvdFile)
    window = defaultWindow(case, frames)
    tStart = window[0] if args.tstart is None else args.tstart
    tEnd = window[1] if args.tend is None else args.tend
    frames = [frame for frame in frames if tStart - 1e-12 <= frame[0] <= tEnd + 1e-12]
    frames = frames[:: args.stride]
    if not frames:
        sys.exit(f"No VTK output in [{tStart}, {tEnd}]")
    return frames, (tStart, tEnd)


def createFigure(fields, ranges, indicators, window):
    """Figure with the field panels and the drag and lift coefficients over time"""
    panelHeight = 10.0 * CHANNEL[1] / CHANNEL[0]
    figure, axes = plt.subplots(
        4,
        1,
        figsize=(10.0, 3 * panelHeight + 2.6),
        gridspec_kw={"height_ratios": [panelHeight] * 3 + [2.0]},
    )
    images = {}
    for axis, (key, label, colormap) in zip(axes, FIELD_PANELS):
        images[key] = axis.imshow(
            fields[key],
            origin="lower",
            extent=(0.0, CHANNEL[0], 0.0, CHANNEL[1]),
            cmap=colormap,
            vmin=ranges[key][0],
            vmax=ranges[key][1],
            interpolation="bilinear",
        )
        axis.set_aspect("equal")
        axis.set_xticks([])
        axis.set_yticks([])
        figure.colorbar(images[key], ax=axis, label=label, fraction=0.03, pad=0.01)

    marker = plotIndicators(axes[3], indicators, window)
    title = figure.suptitle("")
    figure.tight_layout()
    return figure, images, marker, title


def plotIndicators(dragAxis, indicators, window):
    """Drag and lift coefficients over time, returns the marker of the current time"""
    time = np.asarray(indicators["time"])
    inWindow = (time >= window[0]) & (time <= window[1])
    dragAxis.plot(time[inWindow], np.asarray(indicators["cDrag"])[inWindow], color="C0")
    dragAxis.set_ylabel(r"$C_D$", color="C0")
    dragAxis.set_xlabel(r"$t$ [s]")
    dragAxis.set_xlim(*window)
    liftAxis = dragAxis.twinx()
    liftAxis.plot(time[inWindow], np.asarray(indicators["cLift"])[inWindow], color="C3")
    liftAxis.set_ylabel(r"$C_L$", color="C3")
    return dragAxis.axvline(window[0], color="k", linewidth=1.0)


def videoWriter(videoFormat, fps):
    """Writer of the video frames (even frame sizes for H.264)"""
    if videoFormat == "gif":
        return PillowWriter(fps=fps)
    return FFMpegWriter(
        fps=fps,
        codec="libx264",
        extra_args=["-pix_fmt", "yuv420p", "-vf", "scale=trunc(iw/2)*2:trunc(ih/2)*2"],
    )


def main():
    """Render the video of a time-dependent run"""
    args = parseArguments()
    runFolder = os.path.join(runsDir(args.case), args.run)
    pvdFile = os.path.join(runFolder, "vtk", f"{args.run}.pvd")
    if not os.path.exists(pvdFile):
        sys.exit(f"{pvdFile} not found")
    frames, window = selectFrames(args.case, pvdFile, args)
    print(f"{len(frames)} frames in [{frames[0][0]:g}, {frames[-1][0]:g}] s")
    indicators, note = readRun(runFolder, args.run)
    if indicators is None:
        sys.exit(f"No indicators found for {args.run}")
    if note is not None:
        print(note)

    image = samplingImage(args.resolution)
    ranges = colorRanges(frames, image)
    for key, (low, high) in ranges.items():
        print(f"color range {key}: [{low:.3g}, {high:.3g}]")
    figureParts = createFigure(sampleFields(frames[0][1], image), ranges, indicators, window)

    output = args.output or os.path.join(
        RESULTS_DIR, f"video_2d{args.case}_{args.run}.{args.format}"
    )
    os.makedirs(os.path.dirname(os.path.abspath(output)), exist_ok=True)
    renderVideo(frames, image, figureParts, output, args)
    print(f"\nWrote {output}")


def renderVideo(frames, image, figureParts, output, args):
    """Sample each VTK file, update the figure and write the frame"""
    figure, images, marker, suptitle = figureParts
    writer = videoWriter(args.format, args.fps)
    with writer.saving(figure, output, dpi=args.dpi):
        for index, (frameTime, fileName) in enumerate(frames):
            for key, field in sampleFields(fileName, image).items():
                images[key].set_data(field)
            marker.set_xdata([frameTime, frameTime])
            suptitle.set_text(
                f"DFG benchmark {CASE_TITLES[args.case]}, {args.run}, $t = {frameTime:.3f}$ s"
            )
            writer.grab_frame()
            print(f"\rframe {index + 1}/{len(frames)}", end="", flush=True)


if __name__ == "__main__":
    main()
