#!/usr/bin/env python3
# SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

"""
Vessel with a closed tip: grid convergence of the 1d-3d embedded coupling methods towards
the exact solution of each method.

Usage: python3 run_benchmark.py [output directory of the figures]

Runs the executables test_md_embedded_1d3d_1p1p_tip_* in the current directory, prints the
error tables and fails if an observed convergence rate on the two finest grids is below the
expected one. The figures are only produced if an output directory is given.
"""

import math
import subprocess
import sys
import xml.etree.ElementTree as ET

import numpy as np

# executable suffix: label and minimum rates (p_t, p_v, q) between the two finest grids
CASES = {
    "tpfatpfa_average": ("line source, tpfa", (1.35, 1.35, 1.35)),
    "tpfatpfa_surface": ("surface source, tpfa", (1.3, 1.75, 1.75)),
    "tpfatpfa_kernel": ("kernel, tpfa", (1.9, 1.9, 1.9)),
    "boxtpfa_kernel": ("kernel, box", (1.9, 1.9, 1.9)),
}
# cells per unit length; the vessel surface and the tip keep their position relative to the grid
CELLS = [10, 20, 40, 80]
LENGTH = 0.6
RADIUS = 0.1
QUANTITIES = ["p_t", "p_v", "q"]
# color, marker and line style differ per case, so the figures also read in black and white
STYLES = {
    "tpfatpfa_average": dict(color="#2a78d6", marker="o", linestyle="-"),
    "tpfatpfa_surface": dict(color="#eb6834", marker="^", linestyle="--"),
    "tpfatpfa_kernel": dict(color="#1baf7a", marker="s", linestyle="-."),
    "boxtpfa_kernel": dict(color="#eda100", marker="D", linestyle=":"),
}
FONT = {"font.size": 13, "axes.titlesize": 13, "axes.labelsize": 13,
        "xtick.labelsize": 12, "ytick.labelsize": 12, "legend.fontsize": 12}


def run(case, cells, vtk=False):
    """Run one case and return the cell size and the relative errors"""
    vessel_cells = LENGTH * cells / 2
    if abs(vessel_cells - round(vessel_cells)) > 1e-12:
        sys.exit(f"The tip at z = {LENGTH} is not on a cell face for {cells} cells")
    out = subprocess.run(
        [f"./test_md_embedded_1d3d_1p1p_tip_{case}", "params.input",
         "-Tissue.Grid.Cells", f"{cells} {cells} {cells // 2}",
         "-Vessel.Grid.UpperRight", f"0 0 {LENGTH}",
         "-Vessel.Grid.Cells", str(round(vessel_cells)),
         "-SpatialParams.Radius", str(RADIUS),
         "-Vtk.EnableVtkOutput", str(vtk).lower(), "-Vtk.Precision", "Float64",
         "-Vtk.OutputName", f"tip_{case}_{cells}"],
        capture_output=True, text=True, check=True,
    ).stdout
    words = next(line for line in out.splitlines() if line.startswith("Errors:")).split()
    values = dict(zip(words[1::2], map(float, words[2::2])))
    errors = [values[name] for name in QUANTITIES]
    if not all(math.isfinite(e) for e in errors):
        sys.exit(f"Non-finite error norms for {case} with {cells} cells")
    return 2.0 / cells, errors


def read_vtk(filename):
    """Return points, cell connectivity and cell data of an ASCII VTU or VTP file"""
    root = ET.parse(filename).getroot()
    piece = root.find("*/Piece")
    to_array = lambda array: np.array(array.text.split(), dtype=float)
    points = to_array(piece.find("Points/DataArray")).reshape(-1, 3)
    cells = piece.find("Cells") if piece.find("Cells") is not None else piece.find("Lines")
    arrays = {a.get("Name"): to_array(a) for a in cells}
    offsets = arrays["offsets"].astype(int)
    connectivity = np.split(arrays["connectivity"].astype(int), offsets[:-1])
    cell_data = {a.get("Name"): to_array(a) for a in piece.find("CellData")}
    return points, connectivity, cell_data


def slope_triangle(ax, x_fine, y_low, order, width=1.5):
    """Small labeled slope triangle below the data (h decreasing to the right)"""
    x_coarse = width * x_fine
    y_high = y_low * width ** order
    ax.plot([x_coarse, x_fine, x_coarse, x_coarse], [y_high, y_low, y_low, y_high], color="0.2", linewidth=1.2)
    ax.annotate(f"{order:g}", (x_coarse, np.sqrt(y_low * y_high)), textcoords="offset points", xytext=(-5, 0),
                ha="right", va="center", color="0.2")
    ax.annotate("1", (np.sqrt(x_coarse * x_fine), y_low), textcoords="offset points", xytext=(0, -3),
                ha="center", va="top", color="0.2")


def plot_convergence(results, filename):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update(FONT)

    titles = [r"tissue pressure $p_t$", r"vessel pressure $p_v$", r"exchange rate $q$"]
    fig, axes = plt.subplots(1, 3, figsize=(12, 4.6), layout="constrained")
    for k, ax in enumerate(axes):
        for case, rows in results.items():
            h_over_r = [h / RADIUS for h, _ in rows]
            ax.loglog(h_over_r, [errors[k] for _, errors in rows], markersize=7, linewidth=1.8,
                      label=CASES[case][0], **STYLES[case])
        h_over_r = np.array([h / RADIUS for h, _ in next(iter(results.values()))])
        lows = []
        for level, order in [(-1, 2.0), (-2, 1.5)]:
            lows.append(0.3 * min(rows[level][1][k] for rows in results.values()))
            slope_triangle(ax, h_over_r[level], lows[-1], order)
        ax.set_ylim(bottom=0.5 * min(lows))
        ax.invert_xaxis()
        ax.set_xticks(h_over_r, [f"{v:g}" for v in h_over_r])
        ax.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
        ax.set_xlabel(r"$h / R$")
        ax.set_title(titles[k])
        ax.grid(True, which="major", color="0.9", linewidth=0.6)
    axes[0].set_ylabel("relative $L^2$ error")
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="outside lower center", ncol=len(labels), frameon=False)
    fig.savefig(filename, dpi=150)


def plot_solution(cells, filename):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.colors import LinearSegmentedColormap
    plt.rcParams.update(FONT)

    blues = LinearSegmentedColormap.from_list(
        "blues", ["#cde2fb", "#9ec5f4", "#6da7ec", "#3987e5", "#256abf", "#184f95", "#0d366b"])

    fig, axes = plt.subplots(1, 3, figsize=(12.5, 4.8), layout="constrained",
                             gridspec_kw={"width_ratios": [1.3, 1, 1]})

    # tissue pressure in the plane through the vessel axis (cells adjacent to y = 0)
    _, _, data = read_vtk(f"tip_tpfatpfa_kernel_{cells}_3d-00000.vtu")
    shape = (cells // 2, cells, cells)
    numeric = data["p"].reshape(shape)[:, cells // 2, :]
    exact = data["exact pressure"].reshape(shape)[:, cells // 2, :]
    h = 2.0 / cells
    x_faces = np.linspace(-1.0, 1.0, cells + 1)
    z_faces = np.linspace(0.0, 1.0, cells // 2 + 1)
    x_centers, z_centers = 0.5 * (x_faces[1:] + x_faces[:-1]), 0.5 * (z_faces[1:] + z_faces[:-1])
    ax = axes[0]
    image = ax.pcolormesh(x_faces, z_faces, numeric, cmap=blues, shading="flat")
    levels = np.linspace(0.05, numeric.max(), 8)
    ax.contour(x_centers, z_centers, exact, levels=levels, colors="0.1", linewidths=1.0, linestyles="--")
    ax.plot([-RADIUS, -RADIUS, RADIUS, RADIUS], [0.0, LENGTH, LENGTH, 0.0], color="white", linewidth=1.5)
    ax.set_aspect("equal")
    ax.set_xticks([-1.0, -0.5, 0.0, 0.5, 1.0])
    ax.set_xlabel("$x$")
    ax.set_ylabel("$z$")
    ax.set_title(f"$p_t$ at $y = h/2$, kernel, tpfa, $h/R = {h / RADIUS:g}$")
    fig.colorbar(image, ax=ax, shrink=0.8)

    # vessel pressure error and exchange rate along the vessel
    for case in CASES:
        points, connectivity, data = read_vtk(f"tip_{case}_{cells}_1d-00000.vtp")
        s = np.array([points[c, 2].mean() for c in connectivity])
        order = np.argsort(s)
        style = dict(markersize=6, linewidth=1.8, label=CASES[case][0], **STYLES[case])
        axes[1].plot(s[order], (data["p"] - data["exact pressure"])[order], **style)
        axes[2].plot(s[order], data["exchange rate"][order], **style)
    axes[1].axhline(0.0, color="0.1", linewidth=1.0)
    axes[2].axhline(1.0, color="0.1", linewidth=1.0, label="exact solution")
    for ax, title in [(axes[1], r"$p_v - p_{v,e}$"), (axes[2], r"mean exchange rate per cell")]:
        ax.set_xlabel("$s$")
        ax.axvline(LENGTH, color="0.6", linewidth=0.8, linestyle=":")
        ax.set_title(title)
        ax.grid(True, color="0.9", linewidth=0.6)
    handles, labels = axes[2].get_legend_handles_labels()
    fig.legend(handles, labels, loc="outside lower center", ncol=len(labels), frameon=False)
    fig.savefig(filename, dpi=150)


def main():
    image_dir = sys.argv[1] if len(sys.argv) > 1 else None

    results = {}
    failed = []
    for case, (label, minimum_rates) in CASES.items():
        rows = [run(case, cells) for cells in CELLS]
        results[case] = rows
        print(f"\n{label}")
        print("| $h/R$ | $e_{p_t}$ | rate | $e_{p_v}$ | rate | $e_q$ | rate |")
        print("|------|------|------|------|------|------|------|")
        for i, (h, errors) in enumerate(rows):
            cols = [f"{h / RADIUS:g}"]
            for k, error in enumerate(errors):
                rate = "-" if i == 0 else f"{math.log(rows[i - 1][1][k] / error) / math.log(rows[i - 1][0] / h):.2f}"
                cols += [f"{error:.2e}", rate]
            print("| " + " | ".join(cols) + " |")
        for k, name in enumerate(QUANTITIES):
            (h0, e0), (h1, e1) = rows[-2], rows[-1]
            rate = math.log(e0[k] / e1[k]) / math.log(h0 / h1)
            if rate < minimum_rates[k]:
                failed.append(f"{label}: rate {rate:.2f} of {name} below {minimum_rates[k]}")

    if image_dir is not None:
        plot_convergence(results, f"{image_dir}/embedded1d3d_tip_convergence.png")
        cells = 40
        for case in CASES:
            run(case, cells, vtk=True)
        plot_solution(cells, f"{image_dir}/embedded1d3d_tip_solution.png")

    if failed:
        sys.exit("\n".join(failed))
    print("\nAll observed convergence rates reach the expected rates.")


if __name__ == "__main__":
    main()
