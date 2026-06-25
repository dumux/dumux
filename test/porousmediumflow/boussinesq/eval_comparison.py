#!/usr/bin/env python3
"""
Compare vorticity (ψ-C) and pressure (p-C) Boussinesq results.

Produces two figures:
  comparison_concentration.pdf  – side-by-side concentration fields at t≈T_PLOT
                                   (log scale 1e-3..1, grey colourmap)
  comparison_sherwood.pdf       – F_mass and F_grad vs time for both methods
                                   with a reference line at F = 0.017

Usage
-----
  python3 eval_comparison.py <boussinesq_build_dir> [--time T]

  <boussinesq_build_dir>  CMake binary dir for the boussinesq tests,
                          i.e. <build_root>/test/porousmediumflow/boussinesq.
  --time T      Target plot time (default: 8.5).
  --out DIR     Output directory for figures (default: current dir).
"""

import argparse
import os
import sys
import xml.etree.ElementTree as ET

import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.ticker as mticker
import numpy as np
import pandas as pd
import pyvista as pv
import seaborn as sns
from scipy.interpolate import griddata


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def find_vtu(pvd_path, target_time):
    """Return (actual_time, vtu_path) closest to *target_time* in a PVD file."""
    pvd_dir = os.path.dirname(pvd_path)
    tree = ET.parse(pvd_path)
    root = tree.getroot()
    collection = root.find("Collection")
    if collection is None:
        sys.exit(f"No <Collection> found in {pvd_path}")
    entries = [
        (float(ds.get("timestep", 0)), ds.get("file", ""))
        for ds in collection.findall("DataSet")
    ]
    if not entries:
        sys.exit(f"Empty collection in {pvd_path}")
    t_actual, rel_path = min(entries, key=lambda x: abs(x[0] - target_time))
    vtu_path = os.path.join(pvd_dir, rel_path)
    return t_actual, vtu_path


def read_concentration(vtu_path, field_name):
    """Read nodal concentration from a VTU file via pyvista."""
    mesh = pv.read(vtu_path)
    if field_name not in mesh.point_data:
        available = list(mesh.point_data.keys())
        sys.exit(
            f"Field '{field_name}' not found in {vtu_path}.\n"
            f"Available: {available}"
        )
    pts = mesh.points[:, :2]          # (n, 2) – drop z
    conc = np.asarray(mesh.point_data[field_name], dtype=float)
    return pts, conc


def interp_to_grid(pts, values, nx=400, ny=400,
                   x_range=(0, 1), y_range=(0, 1)):
    """Interpolate scattered nodal data to a regular grid."""
    xi = np.linspace(*x_range, nx)
    yi = np.linspace(*y_range, ny)
    Xi, Yi = np.meshgrid(xi, yi)
    Ci = griddata(pts, values, (Xi, Yi), method="linear")
    return Xi, Yi, Ci


# ---------------------------------------------------------------------------
# Figure 1 – concentration field comparison
# ---------------------------------------------------------------------------

def plot_concentration(vort_pts, vort_conc, pres_pts, pres_conc,
                       t_vort, t_pres, outfile):
    cmap = "Greys"          # white = low C, black = high C
    norm = mcolors.LogNorm(vmin=1e-3, vmax=1.0)

    fig, axes = plt.subplots(1, 2, figsize=(9, 4.5),
                             constrained_layout=True)

    panels = [
        (axes[0], vort_pts, vort_conc, t_vort,
         r"Vorticity ($\psi$-C)"),
        (axes[1], pres_pts, pres_conc, t_pres,
         r"Pressure ($p$-C)"),
    ]

    im = None
    for ax, pts, conc, t, title in panels:
        Xi, Yi, Ci = interp_to_grid(pts, conc)
        Ci = np.clip(Ci, 1e-3, 1.0)
        im = ax.pcolormesh(Xi, Yi, Ci, norm=norm, cmap=cmap,
                           shading="auto", rasterized=True)
        ax.set_aspect("equal")
        ax.set_xlabel("$x$", fontsize=11)
        ax.set_ylabel("$z$", fontsize=11)
        ax.set_title(f"{title}\n$t = {t:.3f}$", fontsize=11)
        ax.tick_params(labelsize=9)

    cbar = fig.colorbar(im, ax=axes, orientation="vertical",
                        fraction=0.03, pad=0.02, aspect=30)
    cbar.set_label("Concentration $C$", fontsize=11)
    cbar.ax.yaxis.set_major_formatter(
        mticker.LogFormatterSciNotation(base=10, labelOnlyBase=False))
    cbar.set_ticks([1e-3, 1e-2, 1e-1, 1e0])

    fig.suptitle(
        r"CO$_2$ concentration — Ra = 8000, 200×200, log scale",
        fontsize=12, y=1.01
    )
    fig.savefig(outfile, dpi=200, bbox_inches="tight")
    print(f"Saved: {outfile}")
    plt.close(fig)


# ---------------------------------------------------------------------------
# Figure 2 – Sherwood number comparison (seaborn)
# ---------------------------------------------------------------------------

def plot_sherwood(vort_csv, pres_csv, outfile, f_ref=0.017):
    sns.set_theme(style="whitegrid", context="paper", font_scale=1.2)
    palette = sns.color_palette("colorblind", 4)

    vort = pd.read_csv(vort_csv)
    pres = pd.read_csv(pres_csv)

    fig, ax = plt.subplots(figsize=(7, 4.5))

    # --- vorticity lines ---
    ax.plot(vort["time"], vort["F_grad"], color=palette[0], lw=1.5,
            ls="-", label=r"Vorticity $F_\mathrm{grad}$")
    ax.plot(vort["time"], vort["F_mass"], color=palette[0], lw=1.5,
            ls="--", label=r"Vorticity $F_\mathrm{mass}$")

    # --- pressure lines ---
    ax.plot(pres["time"], pres["F_grad"], color=palette[1], lw=1.5,
            ls="-", label=r"Pressure $F_\mathrm{grad}$")
    ax.plot(pres["time"], pres["F_mass"], color=palette[1], lw=1.5,
            ls="--", label=r"Pressure $F_\mathrm{mass}$")

    # --- reference line ---
    ax.axhline(f_ref, color="0.3", lw=1.2, ls=":",
               label=f"$F = {f_ref}$")

    ax.set_xlabel("Dimensionless time $t$", fontsize=12)
    ax.set_ylabel(r"Dissolution flux $F$", fontsize=12)
    ax.set_title(
        r"Dissolution flux — Ra = 8000, 200×200", fontsize=12
    )

    # legend outside plot area to avoid clutter
    ax.legend(loc="upper right", framealpha=0.9, fontsize=9,
              ncol=2)
    ax.set_xlim(left=0)
    ax.set_ylim(bottom=0)

    sns.despine()
    fig.tight_layout()
    fig.savefig(outfile, dpi=200, bbox_inches="tight")
    print(f"Saved: {outfile}")
    plt.close(fig)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Evaluate vorticity vs. pressure Boussinesq comparison")
    parser.add_argument("build_dir",
                        help="CMake binary dir for boussinesq tests "
                             "(<build_root>/test/porousmediumflow/boussinesq)")
    parser.add_argument("--time", type=float, default=8.5,
                        help="Target plot time (default: 8.5)")
    parser.add_argument("--out", default=".",
                        help="Output directory for figures (default: .)")
    args = parser.parse_args()

    vort_dir = os.path.join(args.build_dir, "vorticity")
    pres_dir = os.path.join(args.build_dir, "pressure")

    os.makedirs(args.out, exist_ok=True)

    # ---- locate PVD and CSV files -----------------------------------------
    vort_pvd = os.path.join(vort_dir, "comparison_vorticity.pvd")
    pres_pvd = os.path.join(pres_dir, "comparison_pressure.pvd")
    vort_csv = os.path.join(vort_dir, "comparison_vorticity_sherwood.csv")
    pres_csv = os.path.join(pres_dir, "comparison_pressure_sherwood.csv")

    for p in [vort_pvd, pres_pvd, vort_csv, pres_csv]:
        if not os.path.isfile(p):
            sys.exit(f"File not found: {p}\nRun run_comparison.py first.")

    # ---- find VTU files at t ≈ args.time -----------------------------------
    t_vort, vtu_vort = find_vtu(vort_pvd, args.time)
    t_pres, vtu_pres = find_vtu(pres_pvd, args.time)
    print(f"Vorticity snapshot: t = {t_vort:.4f}  ({vtu_vort})")
    print(f"Pressure  snapshot: t = {t_pres:.4f}  ({vtu_pres})")

    # ---- read concentration fields -----------------------------------------
    # vorticity model exports "concentration"
    # pressure model exports "X^Solute_liquid"  (OnePNC mass fraction)
    vort_pts, vort_conc = read_concentration(vtu_vort, "concentration")
    pres_pts, pres_conc = read_concentration(vtu_pres, "X^Solute_liquid")

    # ---- Figure 1: concentration fields ------------------------------------
    conc_out = os.path.join(args.out, "comparison_concentration.pdf")
    plot_concentration(vort_pts, vort_conc,
                       pres_pts, pres_conc,
                       t_vort, t_pres,
                       conc_out)

    # ---- Figure 2: Sherwood comparison (seaborn) ---------------------------
    sh_out = os.path.join(args.out, "comparison_sherwood.pdf")
    plot_sherwood(vort_csv, pres_csv, sh_out)


if __name__ == "__main__":
    main()
