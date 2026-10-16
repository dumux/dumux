#!/usr/bin/env python3
# SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

"""
Henry problem of Fahs et al. (2016): steady concentration field with the isochlor positions of Tables D1 and D2,
and the concentration along the midline over time.

Usage: python3 plot_profiles.py test_1p2c_henry_diffusion_box.pvd test_1p2c_henry_dispersion_box.pvd [figure]
"""

import sys
import xml.etree.ElementTree as ET
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.axes import Axes
from matplotlib.cm import ScalarMappable
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.patheffects import withStroke

TABLE_DIR = Path(__file__).resolve().parent
CASES = ("Molecular diffusion", "Velocity-dependent dispersion")
HOURS = (1, 2, 3, 4, 6, 24)
LEVELS = (0.1, 0.5, 0.9)
MIDLINE = 0.5
FIELD_COLORS = LinearSegmentedColormap.from_list("field", ["#ffffff", "#cde9ef", "#5fb4c6", "#00849e", "#0b4f63"])
TIME_COLORS = LinearSegmentedColormap.from_list("time", ["#bfe3ea", "#6cbccc", "#00849e", "#0b4f63", "#08303d"])
LINE_COLOR = "#08303d"
REFERENCE = dict(marker="o", ls="none", mfc="none", mec="#c2185b")
HALO = [withStroke(linewidth=3.0, foreground="white")]

Snapshot = tuple[np.ndarray, np.ndarray]


def read_vtu(path: Path) -> Snapshot:
    """Points and solute mass fraction of an ASCII VTU file"""
    piece = ET.parse(path).getroot().find("UnstructuredGrid/Piece")
    array = lambda name: np.array(piece.find(f".//DataArray[@Name='{name}']").text.split(), dtype=float)
    return array("Coordinates").reshape(-1, 3), array("X^Salt_liq")


def read_output(pvd: Path) -> dict[float, Snapshot]:
    """Output closest to the given hours, keyed by its time in hours"""
    files = {float(d.get("timestep"))/3600: pvd.parent/d.get("file") for d in ET.parse(pvd).getroot().iter("DataSet")}
    times = np.array(list(files))
    return {t: read_vtu(files[t]) for t in (times[np.abs(times - h).argmin()] for h in HOURS)}


def read_table(case: int) -> np.ndarray:
    """Rows (Z, X10, X50, X90) of Table D1 or D2"""
    lines = (TABLE_DIR/f"fahs2016_table_d{case}.csv").read_text().splitlines()
    return np.loadtxt([line for line in lines if line[:1].isdigit()], delimiter=",")


def plot_field(ax: Axes, snapshot: Snapshot, table: np.ndarray) -> ScalarMappable:
    """Concentration field with the computed isochlors, the reference positions and the midline"""
    points, c = snapshot
    x, z = points[:, 0], points[:, 1]
    field = ax.tricontourf(x, z, c, levels=np.linspace(0.0, 1.0, 21), cmap=FIELD_COLORS)
    ax.tricontour(x, z, c, levels=LEVELS, colors=LINE_COLOR, linewidths=1.2).set(path_effects=HALO)
    ax.plot(table[:, 1:].ravel(), np.repeat(table[:, 0], 3), ms=4, mew=1.2, **REFERENCE)
    ax.axhline(MIDLINE, color=LINE_COLOR, lw=1.0, ls=(0, (4, 3)), path_effects=HALO)
    ax.set(xlim=(0, 3), ylim=(0, 1), yticks=[0, 0.5, 1], xlabel="x (m)")
    return field


def plot_midline(ax: Axes, output: dict[float, Snapshot], table: np.ndarray) -> None:
    """Concentration along the midline over time and the reference positions"""
    colors = TIME_COLORS(np.linspace(0.0, 1.0, len(output)))
    widths = np.linspace(1.0, 2.5, len(output))
    for color, width, (t, (points, c)) in zip(colors, widths, output.items()):
        on_line = np.isclose(points[:, 1], MIDLINE)
        order = np.argsort(points[on_line, 0])
        ax.plot(points[on_line, 0][order], c[on_line][order], color=color, lw=width, label=f"t = {t:g} h")
    row = table[np.isclose(table[:, 0], MIDLINE)][0]
    ax.plot(row[1:], LEVELS, ms=6, mew=1.4, label="Fahs et al. (2016)", **REFERENCE)
    ax.text(0.02, 0.92, f"z = {MIDLINE:g} m", transform=ax.transAxes, va="top")
    ax.set(xlim=(0, 3), ylim=(-0.03, 1.03), yticks=[0, 0.5, 1], xlabel="x (m)")
    ax.grid(color="#e4e3df", lw=0.5)


def plot(pvd_files: list[str], filename: str) -> None:
    plt.rcParams.update({"font.size": 11, "axes.titlesize": 12, "xtick.labelsize": 10, "ytick.labelsize": 10,
                         "legend.fontsize": 10})
    # the height of the first row shows the 3:1 domain to scale
    fig, axes = plt.subplots(2, 2, figsize=(9.0, 4.3), layout="constrained", gridspec_kw={"height_ratios": [0.8, 1.0]})
    for col, (pvd, case) in enumerate(zip(pvd_files, CASES)):
        output = read_output(Path(pvd))
        table = read_table(col + 1)
        field = plot_field(axes[0, col], output[max(output)], table)
        axes[0, col].set_title(f"{case}, t = {max(output):g} h")
        plot_midline(axes[1, col], output, table)
    axes[0, 0].set_ylabel("z (m)")
    axes[1, 0].set_ylabel("c (-)")
    axes[1, 1].tick_params(labelleft=False)
    fig.colorbar(field, ax=axes[0, 1], ticks=[0, 0.5, 1], fraction=0.04, pad=0.02, label="c (-)")
    fig.legend(*axes[1, 1].get_legend_handles_labels(), loc="outside lower center", frameon=False, ncol=len(HOURS) + 1)
    fig.savefig(filename, dpi=150)


if __name__ == "__main__":
    if len(sys.argv) < 3:
        sys.exit(__doc__)
    plot(sys.argv[1:3], sys.argv[3] if len(sys.argv) > 3 else "henry_profiles.png")
