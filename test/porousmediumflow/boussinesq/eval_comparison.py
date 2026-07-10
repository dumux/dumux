#!/usr/bin/env python3
"""
Compare vorticity (ψ-C) and pressure (p-C) Boussinesq results.

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
    sns.set_theme(style="white", context="paper", font_scale=1.2)
    cmap = "mako_r"
    norm = mcolors.LogNorm(vmin=1e-3, vmax=1.0)

    fig, axes = plt.subplots(1, 2, figsize=(9, 4.5),
                             constrained_layout=True)

    panels = [
        (axes[0], vort_pts, vort_conc, t_vort,
         r"Vorticity ($\psi$-$C$)"),
        (axes[1], pres_pts, pres_conc, t_pres,
         r"Pressure ($p$-$C$)"),
    ]

    im = None
    for ax, pts, conc, t, title in panels:
        Xi, Yi, Ci = interp_to_grid(pts, conc)
        Ci = np.clip(Ci, 1e-3, 1.0)
        im = ax.pcolormesh(Xi, Yi, Ci, norm=norm, cmap=cmap,
                           shading="auto", rasterized=True)
        ax.set_aspect("equal")
        ax.set_xlabel("$x$")
        ax.set_ylabel("$z$")
        ax.set_title(f"{title}\n$t = {t:.2f}$")
        sns.despine(ax=ax, left=False, bottom=False)

    cbar = fig.colorbar(im, ax=axes, orientation="vertical",
                        fraction=0.03, pad=0.02, aspect=30)
    cbar.set_label("Concentration $C$")
    cbar.set_ticks([1e-3, 1e-2, 1e-1, 1e0])
    cbar.ax.yaxis.set_major_formatter(mticker.LogFormatterMathtext())

    fig.suptitle(r"CO$_2$ concentration — Ra = 8000",
                 fontsize=12)

    fig.savefig(outfile, dpi=200, bbox_inches="tight")
    plt.close(fig)


# ---------------------------------------------------------------------------
# Figure 2 – Sherwood number comparison (seaborn)
# ---------------------------------------------------------------------------

def plot_sherwood(vort_csv, pres_csv, outfile, f_ref=0.017):
    sns.set_theme(style="whitegrid", context="paper", font_scale=1.2)
    c_vort, c_pres = sns.color_palette("deep", 2)

    vort = pd.read_csv(vort_csv)
    pres = pd.read_csv(pres_csv)

    fig, ax = plt.subplots(figsize=(7, 4.5))

    # --- vorticity lines ---
    ax.plot(vort["time"], vort["F_grad"], color=c_vort, lw=1.8,
            ls="-",  label=r"Vorticity $F_\mathrm{grad}$")
    ax.plot(vort["time"], vort["F_mass"], color=c_vort, lw=1.8,
            ls="--", label=r"Vorticity $F_\mathrm{mass}$")

    # --- pressure lines ---
    ax.plot(pres["time"], pres["F_grad"], color=c_pres, lw=1.8,
            ls="-",  label=r"Pressure $F_\mathrm{grad}$")
    ax.plot(pres["time"], pres["F_mass"], color=c_pres, lw=1.8,
            ls="--", label=r"Pressure $F_\mathrm{mass}$")

    # --- reference line ---
    ax.axhline(f_ref, color="0.45", lw=1.2, ls=":",
               label=f"$F_{{\\mathrm{{ref}}}} = {f_ref}$")

    ax.set_yscale("log")
    ax.set_xscale("log")
    ax.set_xlabel("Dimensionless time $t$")
    ax.set_ylabel(r"Dissolution flux $F$")
    ax.set_title(r"Dissolution flux — Ra = 8000")
    ax.set_xlim(left=0)
    ax.legend(ncol=2, framealpha=0.9, fontsize=9)

    sns.despine()
    fig.tight_layout()

    fig.savefig(outfile, dpi=200, bbox_inches="tight")
    plt.close(fig)


# ---------------------------------------------------------------------------
# Figure 3 – y-velocity field comparison (pyvista off-screen)
# ---------------------------------------------------------------------------

def plot_velocity_pyvista(vtu_vort, vtu_pres, t_vort, t_pres, outfile):
    mesh_vort = pv.read(vtu_vort)
    mesh_pres = pv.read(vtu_pres)

    def extract_vy(mesh, label):
        for name, arr in mesh.point_data.items():
            arr = np.asarray(arr)
            if arr.ndim == 2 and arr.shape[1] >= 2:
                return arr[:, 1]
        sys.exit(f"No vector field found in {label}. "
                 f"Available: {list(mesh.point_data.keys())}")

    vy_vort = extract_vy(mesh_vort, "vorticity VTU")
    vy_pres = extract_vy(mesh_pres, "pressure VTU")

    v_max = max(np.abs(vy_vort).max(), np.abs(vy_pres).max())
    clim  = [-v_max, v_max]

    mesh_vort["v_y"] = vy_vort
    mesh_pres["v_y"] = vy_pres

    sargs = dict(title="$v_y$", title_font_size=18, label_font_size=14,
                 vertical=True, width=0.05, height=0.55,
                 position_x=0.92, position_y=0.22,
                 color="black", fmt="%.3f")

    pl = pv.Plotter(shape=(1, 2), off_screen=True,
                    window_size=[1500, 650], border=False)

    for i, (mesh, t, title) in enumerate([
        (mesh_vort, t_vort, f"Vorticity  (ψ-C),  t = {t_vort:.2f}"),
        (mesh_pres, t_pres, f"Pressure   (p-C),  t = {t_pres:.2f}"),
    ]):
        pl.subplot(0, i)
        pl.background_color = "white"
        pl.add_mesh(mesh, scalars="v_y", cmap="RdBu_r", clim=clim,
                    show_scalar_bar=(i == 1),
                    scalar_bar_args=sargs if i == 1 else {})
        pl.add_title(title, font_size=11, color="black")
        pl.view_xy()
        pl.enable_parallel_projection()
        pl.reset_camera()

    pl.screenshot(outfile, transparent_background=False)
    pl.close()
    print(f"Saved: {outfile}")


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

    # both formulations now build in the same application dir, distinguished
    # by filename prefix (-Problem.Name comparison_vorticity/comparison_pressure)
    vort_dir = os.path.join(args.build_dir, "white-noise-perturbations")
    pres_dir = os.path.join(args.build_dir, "white-noise-perturbations")

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
    conc_out = os.path.join(args.out, "comparison_concentration.png")
    plot_concentration(vort_pts, vort_conc,
                       pres_pts, pres_conc,
                       t_vort, t_pres,
                       conc_out)

    # ---- Figure 2: Sherwood comparison (seaborn) ---------------------------
    sh_out = os.path.join(args.out, "comparison_sherwood.png")
    plot_sherwood(vort_csv, pres_csv, sh_out)

    # ---- Figure 3: y-velocity field (pyvista) ------------------------------
    vy_out = os.path.join(args.out, "comparison_velocity_y.png")
    plot_velocity_pyvista(vtu_vort, vtu_pres, t_vort, t_pres, vy_out)


if __name__ == "__main__":
    main()
