#!/usr/bin/env python3
# SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

"""
Collapse of a thin ring under external pressure: contact force and onset of line contact for
the meshes below against the inextensible elastica (elastica.py). Prints the table of the benchmark
page and, if an output directory is given, writes the figures. Runs the meshes in parallel,
about twelve minutes in total.

usage: run_benchmark.py [output directory of the figures]
"""

import os
import subprocess
import sys
import xml.etree.ElementTree as ET

import numpy as np

from elastica import reference

EXECUTABLE = "./test_contact_ring"
MESHES = [(640, 4), (1280, 8), (2560, 16)]
PRESSURES = [5.5 + 0.5 * k for k in range(14)]
SHAPE_PRESSURES = [6.0, 9.0, 12.0]
RADIUS, THICKNESS = 1.0, 0.04


def start(cells, layers):
    return subprocess.Popen(
        [EXECUTABLE, "params.input", "-Grid.Cells", str(cells), "-Grid.Layers", str(layers),
         "-Problem.Name", f"ring_{cells}x{layers}", "-Problem.Pressures",
         " ".join(map(str, PRESSURES))],
        stdout=subprocess.PIPE,
        text=True,
        env=dict(os.environ, DUMUX_NUM_THREADS="1"),
    )


def parse(output):
    rows = [[float(v) for v in line.split()] for line in output.splitlines()
            if line[:1].isdigit() and len(line.split()) == 7]
    return np.array(rows)


def fit(p, force):
    """Slope and zero of the linear fit of the contact force over 6 <= p <= 9"""
    selected = (p >= 6.0) & (p <= 9.0)
    slope, intercept = np.polyfit(p[selected], force[selected], 1)
    return slope, -intercept / slope


def read_vtu(filename):
    piece = ET.parse(filename).getroot().find("UnstructuredGrid/Piece")

    def to_array(array):
        return np.array(array.text.split(), dtype=float).reshape(
            -1, int(array.get("NumberOfComponents", "1")))
    points = to_array(piece.find("Points/DataArray"))
    pointData = {array.get("Name"): to_array(array) for array in piece.find("PointData")}
    return points, pointData


def plot(results, elastica, onset, outputDir):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    plt.rcParams.update({"font.size": 13, "xtick.labelsize": 12, "ytick.labelsize": 12,
                         "legend.fontsize": 12})
    styles = [("o", "--", "#1f77b4"), ("s", ":", "#d62728"), ("^", "-.", "#2ca02c")]

    fig, axes = plt.subplots(1, 2, figsize=(12, 4.6), layout="constrained")
    order = np.argsort(elastica[:, 0])
    pRef, qRef = elastica[order, 0], elastica[order, 1]
    axes[0].plot(pRef, qRef, "k-", lw=1.8, label="elastica")
    axes[0].plot([onset[0]], [onset[1]], "kD", ms=7, label="elastica: line contact starts")
    axes[1].axhline(0.0, color="k", lw=1.8)
    for (marker, line, color), (cells, layers) in zip(styles, MESHES):
        rows = results[(cells, layers)]
        axes[0].plot(rows[:, 0], rows[:, 4], marker=marker, ls=line, color=color, mfc="none",
                     label=f"{cells} × {layers} cells")
        axes[1].plot(rows[:, 0], rows[:, 4] - np.interp(rows[:, 0], pRef, qRef), marker=marker,
                     ls=line, color=color, mfc="none")
    axes[0].set_ylabel(r"$Q R^2/D$")
    axes[0].set_ylim(-0.3, 11.5)
    axes[1].set_ylabel(r"$(Q - Q_\mathrm{elastica}) R^2/D$")
    for ax in axes:
        ax.set_xlabel(r"$p R^3/D$")
        ax.set_xlim(5.0, 12.2)
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="outside lower center", ncol=5)
    fig.savefig(os.path.join(outputDir, "ring_contactforce.png"), dpi=150)

    cells, layers = MESHES[-1]
    fig, axes = plt.subplots(1, len(SHAPE_PRESSURES), figsize=(12, 3.2), layout="constrained",
                             sharey=True)
    for ax, pressure in zip(axes, SHAPE_PRESSURES):
        index = PRESSURES.index(pressure) + 1
        points, data = read_vtu(f"ring_{cells}x{layers}-{index:05d}.vtu")
        reference = points[:, :2]
        current = reference + data["d"][:, :2]
        radius = np.linalg.norm(reference, axis=1)
        for surface in [radius < RADIUS - 0.49 * THICKNESS, radius > RADIUS + 0.49 * THICKNESS]:
            nodes = np.flatnonzero(surface)
            nodes = nodes[np.argsort(np.arctan2(reference[nodes, 1], reference[nodes, 0]))]
            ax.plot(current[nodes, 0], current[nodes, 1], "k-", lw=1.2)
        inContact = np.ravel(data["contactPressure"]) > 0.0
        ax.plot(current[inContact, 0], current[inContact, 1], "o", color="#d62728", ms=3,
                label="nodes in contact")
        ax.set_aspect("equal")
        ax.set_xlim(-1.45, 1.45)
        ax.set_ylim(-0.4, 0.4)
        ax.set_xlabel("$x/R$")
        ax.set_title(f"$p R^3/D = {pressure:g}$")
    axes[0].set_ylabel("$y/R$")
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="outside lower center")
    fig.savefig(os.path.join(outputDir, "ring_shapes.png"), dpi=150)


def main():
    halfThickness = 0.5 * THICKNESS / RADIUS
    processes = {mesh: start(*mesh) for mesh in MESHES}
    elastica, first, onset = reference(halfThickness)
    contact = elastica[elastica[:, 1] > 0.0]
    slopeRef, zeroRef = fit(contact[:, 0], contact[:, 1])

    results, failed = {}, False
    for mesh, process in processes.items():
        output, _ = process.communicate()
        if process.returncode != 0:
            sys.exit(f"run with {mesh[0]} x {mesh[1]} cells failed")
        results[mesh] = parse(output)
        failed = failed or not np.all(np.isfinite(results[mesh]))

    print("| Mesh | $dQ/dp$ | deviation | zero of the fit | deviation | line contact starts |"
          " smallest semi-minor axis |")
    print("| :--- | ---: | ---: | ---: | ---: | :--- | ---: |")
    print(f"| elastica | {slopeRef:.4f} | | {zeroRef:.4f} | | {onset[0]:.2f} | {halfThickness} |")
    for (cells, layers), rows in results.items():
        p, semiMinor, force, peak = rows[:, 0], rows[:, 1], rows[:, 4], rows[:, 5]
        slope, zero = fit(p, force)
        offAxis = np.flatnonzero(peak > 1e-8)
        lineContact = (f"{p[offAxis[0] - 1]:g} to {p[offAxis[0]]:g}"
                       if len(offAxis) and offAxis[0] > 0 else "not reached")
        failed = failed or semiMinor.min() < halfThickness
        print(f"| {cells} × {layers} | {slope:.4f} | {100 * (slope / slopeRef - 1):+.2f} % | "
              f"{zero:.4f} | {100 * (zero / zeroRef - 1):+.2f} % | {lineContact} | "
              f"{semiMinor.min():.5f} |")
    print(f"\nelastica: first contact at p = {first[0]:.4f}, line contact at p = {onset[0]:.4f}")

    if len(sys.argv) > 1:
        plot(results, elastica, onset, sys.argv[1])
    if failed:
        sys.exit("non-finite results or crossing walls")


if __name__ == "__main__":
    main()
