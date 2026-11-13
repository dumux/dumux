#!/usr/bin/env python3
# SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later
"""
Create the overview GIF of the 1pni heat conduction benchmark: numerical temperature profile over twice the benchmark end time.

Usage (run from the build directory of the test, which contains params.input):
  python3 make_gif.py <executable> [--output conduction.gif]
"""

import argparse
import glob
import os
import subprocess
import xml.etree.ElementTree as ET

from plot_results import (BASE_NAME, INPUT_FILE, TEMPERATURE_FIELD, TEMPERATURE_LIMITS,
                          X_MAX, format_time, read_params, require_plot_modules, sorted_profile)

NUM_FRAMES = 30
FRAME_MS = 100
NUMERICAL_COLOR = "#1676CA"
TIME_FACTOR = 2  # simulate TIME_FACTOR times the end time of the input file


def run_outputs(executable):
    for file_name in glob.glob(f"{BASE_NAME}_gif*"):
        os.remove(file_name)
    name = f"{BASE_NAME}_gif"
    exe = executable if os.path.isabs(executable) else "./" + executable
    end_time = TIME_FACTOR*float(read_params()["TimeLoop.TEnd"].split()[0])
    subprocess.run([exe, INPUT_FILE, "-TimeLoop.TEnd", str(end_time), "-Problem.Name", name, "-Problem.OutputInterval", "5",
                    "-Vtk.Precision", "Float64"], check=True, stdout=subprocess.DEVNULL)
    collection = ET.parse(f"{name}.pvd").getroot()
    return sorted((float(d.attrib["timestep"]), d.attrib["file"]) for d in collection.iter("DataSet"))


def frame_times(results, count):
    """Equidistant in sqrt(t): the diffusion front then advances at a steady pace."""
    transient = [r for r in results if r[0] > 0.0]
    first, last = transient[0][0]**0.5, transient[-1][0]**0.5
    targets = [(first + (last - first)*i/(count - 1))**2 for i in range(count)]
    picked = sorted({min(range(len(transient)), key=lambda k: abs(transient[k][0] - t)) for t in targets})
    return [transient[k] for k in picked]


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("executable")
    parser.add_argument("--output", default="conduction.gif")
    args = parser.parse_args()

    pv, plt = require_plot_modules()
    from PIL import Image
    import io

    frames = []
    for time, vtu in frame_times(run_outputs(args.executable), NUM_FRAMES):
        mesh = pv.read(vtu)
        x, temperature = sorted_profile(mesh, TEMPERATURE_FIELD)

        fig, ax = plt.subplots(figsize=(4.8, 3.0), dpi=100, constrained_layout=True)
        ax.plot(x, temperature, color=NUMERICAL_COLOR, linewidth=2.2)
        ax.set_xlim(0, X_MAX)
        ax.set_ylim(*TEMPERATURE_LIMITS)
        ax.set_xlabel("x (m)")
        ax.set_ylabel("Temperature $T$ (K)")
        ax.grid(True, color="#D6D6D6", linewidth=0.6)
        ax.text(0.97, 0.93, f"$t$ = {format_time(time)}", transform=ax.transAxes, ha="right", va="top")
        buffer = io.BytesIO()
        fig.savefig(buffer, format="png")
        plt.close(fig)
        buffer.seek(0)
        frames.append(Image.open(buffer).convert("RGB").quantize(16, dither=0))

    frames[0].save(args.output, save_all=True, append_images=frames[1:] + [frames[-1]]*8,
                   duration=FRAME_MS, loop=0, optimize=True)
    print(f"Wrote {args.output} ({len(frames)} frames, {os.path.getsize(args.output)/1024:.0f} kB)")


if __name__ == "__main__":
    main()
