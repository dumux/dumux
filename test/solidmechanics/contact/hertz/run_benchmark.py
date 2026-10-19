#!/usr/bin/env python3
# SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

"""
Hertz contact of two cylinders: contact pressure against Hertz's solution for three meshes and
four approaches of the platens. Prints the tables of the benchmark page and, if an output
directory is given, writes the figures. Runs the meshes in parallel, about an hour in total.

usage: run_benchmark.py <gmsh> [output directory of the figures]
"""

import math
import os
import subprocess
import sys

EXECUTABLE = "./test_contact_hertz"
MESH_SIZES = [0.005, 0.0025, 0.00125]
APPROACHES = [0.00125, 0.0025, 0.005, 0.01]
POISSON_RATIO = 0.3


def start(gmsh, meshSize):
    mesh = f"hertz_hc{meshSize}.msh"
    subprocess.run(
        [gmsh, "-2", "hertz.geo", "-setnumber", "hc", str(meshSize), "-o", mesh],
        check=True,
        stdout=subprocess.DEVNULL,
    )
    return subprocess.Popen(
        [EXECUTABLE, "params.input", "-Grid.File", mesh, "-Problem.Name", f"hertz_hc{meshSize}",
         "-Problem.Approaches", " ".join(map(str, APPROACHES)), "-Problem.PrintProfiles", "true",
         "-Vtk.Disable", "true"],
        stdout=subprocess.PIPE,
        text=True,
        env=dict(os.environ, DUMUX_NUM_THREADS="1"),
    )


def parse(output):
    rows, profiles = {}, {}
    for line in output.splitlines():
        values = line.split()
        if line[:1].isdigit() and len(values) == 7:
            approach, force, aOverR, p0, peak, profile, pairs = map(float, values)
            rows[approach] = dict(force=force, aOverR=aOverR, peak=peak, profile=profile)
        elif values and values[0] == "profile":
            profiles.setdefault(float(values[1]), []).append(tuple(map(float, values[2:6])))
    for approach, profile in profiles.items():
        rows[approach]["peakPenalty"] = max(v[3] for v in profile) - 1.0
    return rows, profiles


def plot(results, outputDir):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import numpy as np

    plt.rcParams.update({"font.size": 13, "xtick.labelsize": 12, "ytick.labelsize": 12,
                         "legend.fontsize": 12})
    styles = [("o", "-", "#1f77b4"), ("s", "--", "#d62728"), ("^", ":", "#2ca02c")]

    finest = results[MESH_SIZES[-1]][1]
    fig, ax = plt.subplots(figsize=(7, 4.2), layout="constrained")
    x = np.linspace(-1.5, 1.5, 601)
    ax.plot(x, np.sqrt(np.clip(1.0 - x**2, 0.0, None)), "k-", lw=1.8, label="Hertz")
    for (marker, _, color), approach in zip(styles, [APPROACHES[1], APPROACHES[-1]]):
        aOverR = results[MESH_SIZES[-1]][0][approach]["aOverR"]
        xs, ps = np.array(sorted(v[:2] for v in finest[approach])).T
        ax.plot(xs, ps, marker, color=color, ms=4.5, mfc="none", ls="none",
                label=f"$a/R = {aOverR:.3f}$")
    ax.set_xlabel("$x/a$")
    ax.set_ylabel("$p/p_0$")
    ax.set_xlim(-1.5, 1.5)
    ax.set_ylim(-0.05, 1.12)
    ax.legend(loc="upper right")
    fig.savefig(os.path.join(outputDir, "hertz_profiles.png"), dpi=150)

    fig, axes = plt.subplots(1, 2, figsize=(11, 4.4), layout="constrained")
    aOverR = np.array([results[MESH_SIZES[-1]][0][d]["aOverR"] for d in APPROACHES])
    for (marker, line, color), meshSize in zip(styles, MESH_SIZES):
        rows = results[meshSize][0]
        for ax, key in zip(axes, ["peak", "profile"]):
            ax.plot([rows[d]["aOverR"] for d in APPROACHES],
                    [100 * rows[d][key] for d in APPROACHES],
                    marker=marker, ls=line, color=color, mfc="none",
                    label=f"$h_c = {meshSize}\\,R$")
    rows = results[MESH_SIZES[-1]][0]
    axes[0].plot(aOverR, [100 * rows[d]["peakPenalty"] for d in APPROACHES], marker="x",
                 ls="-.", color="k", label=f"penalty pressure, $h_c = {MESH_SIZES[-1]}\\,R$")
    for ax in axes:
        ax.set_xlabel("$a/R$")
        ax.set_xticks(np.round(aOverR, 3))
        ax.axhline(0.0, color="0.6", lw=0.8)
    axes[0].set_ylabel(r"$p_\mathrm{max}/p_0 - 1$ [%]")
    axes[1].set_ylabel("relative $L^2$ error [%]")
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="outside lower center", ncol=4)
    fig.savefig(os.path.join(outputDir, "hertz_errors.png"), dpi=150)


def main():
    if len(sys.argv) < 2:
        sys.exit(__doc__)
    gmsh = sys.argv[1]
    processes = {h: start(gmsh, h) for h in MESH_SIZES}
    results = {}
    for meshSize, process in processes.items():
        output, _ = process.communicate()
        if process.returncode != 0:
            sys.exit(f"run with h_c = {meshSize} failed")
        results[meshSize] = parse(output)

    failed = False
    titles = {"peak": "peak pressure p_max/p0 - 1", "profile": "relative L2 error",
              "peakPenalty": "peak penalty pressure (force per reference area) / p0 - 1"}
    for key, title in titles.items():
        print(f"\n{title}")
        print("| $a/R$ | $p_0/E$ | " + " | ".join(f"$h_c = {h}\\,R$" for h in MESH_SIZES) + " |")
        print("| ---: | ---: |" + " ---: |" * len(MESH_SIZES))
        for approach in APPROACHES:
            aOverR = results[MESH_SIZES[-1]][0][approach]["aOverR"]
            values = [results[h][0][approach][key] for h in MESH_SIZES]
            failed = failed or not all(math.isfinite(v) for v in values)
            strain = aOverR / (2.0 * (1.0 - POISSON_RATIO**2))
            cells = [f"{100 * v:.2f} %" if key == "profile" else f"{100 * v:+.2f} %"
                     for v in values]
            print(f"| {aOverR:.4f} | {strain:.4f} | " + " | ".join(cells) + " |")

    print("\nforce per unit length F / (E R)")
    print("| approach $\\Delta/R$ | " + " | ".join(f"$h_c = {h}\\,R$" for h in MESH_SIZES) + " |")
    print("| ---: |" + " ---: |" * len(MESH_SIZES))
    for approach in APPROACHES:
        print(f"| {approach} | " + " | ".join(f"{results[h][0][approach]['force']:.6f}"
              for h in MESH_SIZES) + " |")

    if len(sys.argv) > 2:
        plot(results, sys.argv[2])
    if failed:
        sys.exit("non-finite results")


if __name__ == "__main__":
    main()
