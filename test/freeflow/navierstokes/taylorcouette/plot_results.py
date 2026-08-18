#!/usr/bin/env pvpython
# SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later
"""
plot_results.py
================

Post-processing for a single Taylor-Couette benchmark run.
Always generates an analytical comparison plot using the
pressureExact and velocityExact fields already written to the
VTK output; optionally also renders 2D field maps.

Features
--------
1. Analytical comparison plot (always generated)
2. Optional 2D field rendering (pressure, velocity)

Examples
--------
# Point at the directory containing test_ff_taylorcouette.pvd
# (e.g. your build output directory)
pvpython plot_results.py /path/to/build-cmake/test/freeflow/navierstokes/taylorcouette

# Also generate 2D field maps
pvpython plot_results.py /path/to/build-cmake/test/freeflow/navierstokes/taylorcouette --fields
"""

from __future__ import annotations

import argparse
import re
import shutil
import tempfile
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

# ParaView
from paraview import servermanager
from paraview.simple import OpenDataFile, Delete

# =============================================================================
# Configuration
# =============================================================================

ROOT_DIR = Path(__file__).resolve().parent
DEFAULT_RESULTS_DIR = ROOT_DIR / "results"

DPI = 150
VIEW_SIZE = [2000, 1600]

# =============================================================================
# Case Discovery
# =============================================================================


def find_case(results_dir: Path):
    """
    Look for the DuMux VTK output in the given directory.
    """

    pvd_file = results_dir / "test_ff_taylorcouette.pvd"

    if not pvd_file.exists():
        raise RuntimeError(
            f"No test_ff_taylorcouette.pvd found in: {results_dir}"
        )

    return {"path": results_dir, "pvd": pvd_file}


# =============================================================================
# Extraction Utilities
# =============================================================================


def prepare_pvd_for_paraview(case_dir: Path, tmpdir: Path):
    """
    Copy the PVD and its referenced VTU files into a working directory
    and clean the paths inside the PVD, so ParaView can open it
    regardless of where the original run wrote its output.
    """

    pvd_src = case_dir / "test_ff_taylorcouette.pvd"

    if not pvd_src.exists():
        raise RuntimeError(f"No PVD file found in: {case_dir}")

    for vtu_src in case_dir.glob("test_ff_taylorcouette*.vtu"):
        shutil.copy(vtu_src, tmpdir / vtu_src.name)

    pvd_path = tmpdir / pvd_src.name
    shutil.copy(pvd_src, pvd_path)

    pvd_content = pvd_path.read_text()

    cleaned_content = re.sub(
        r'file="[^"]*/([^"]+\.vtu)"',
        r'file="\1"',
        pvd_content,
    )

    pvd_path.write_text(cleaned_content)

    return pvd_path


# =============================================================================
# Data Extraction
# =============================================================================


def parse_data_via_vtk(pvd_filepath: Path):
    """
    Extract pressure, velocity, and cell centers using ParaView.
    """

    data_reader = OpenDataFile(str(pvd_filepath))
    data_reader.UpdatePipeline(1.0)

    client_data = servermanager.Fetch(data_reader)

    vtk_p_array = client_data.GetCellData().GetArray("p")
    vtk_v_array = client_data.GetCellData().GetArray("velocity_liq (m/s)")
    vtk_p_exact_array = client_data.GetCellData().GetArray("pressureExact")
    vtk_v_exact_array = client_data.GetCellData().GetArray("velocityExact")

    num_cells = client_data.GetNumberOfCells()

    pressure = np.array([
        vtk_p_array.GetValue(i)
        for i in range(num_cells)
    ])

    velocity = np.array([
        vtk_v_array.GetTuple(i)
        for i in range(num_cells)
    ])

    pressure_exact = np.array([
        vtk_p_exact_array.GetValue(i)
        for i in range(num_cells)
    ])

    velocity_exact = np.array([
        vtk_v_exact_array.GetTuple(i)
        for i in range(num_cells)
    ])

    cell_centers = []

    for i in range(num_cells):

        cell = client_data.GetCell(i)
        pts = cell.GetPoints()

        num_pts = pts.GetNumberOfPoints()

        coords = np.array([
            pts.GetPoint(j)
            for j in range(num_pts)
        ])

        cell_centers.append(coords.mean(axis=0))

    cell_centers = np.array(cell_centers)

    x = cell_centers[:, 0]
    y = cell_centers[:, 1]

    Delete(data_reader)

    num_nan_p_exact = int(np.isnan(pressure_exact).sum())
    num_nan_v_exact = int(np.isnan(velocity_exact).any(axis=1).sum())
    print(f"  Diagnostic: {num_nan_p_exact} / {num_cells} cells have NaN pressureExact")
    print(f"  Diagnostic: {num_nan_v_exact} / {num_cells} cells have NaN velocityExact")

    return x, y, pressure, velocity, pressure_exact, velocity_exact


# =============================================================================
# Radial Profiles
# =============================================================================


def _bin_radially(r_filtered, values_filtered, bins, n_bins):
    """
    Helper: average a per-cell scalar field into radial bins.
    """

    idx = np.clip(
        np.digitize(r_filtered, bins) - 1,
        0,
        n_bins - 1,
    )

    return np.array([
        values_filtered[idx == i].mean()
        if (idx == i).any()
        else np.nan
        for i in range(n_bins)
    ])


def sorted_exact_profile(cx, cy, pressure_exact, velocity_exact, r1=1.0, r2=2.0):
    """
    Return the analytical fields directly from the per-cell VTK data,
    sorted by radius (no binning). The analytical field is smooth by
    construction, so unlike the numerical solution it doesn't need
    binning to average out angular scatter -- and binning it onto a
    fixed-width grid can leave bins empty wherever the underlying mesh
    is coarser than the bin width, producing spurious gaps.
    """

    r_cells = np.sqrt(cx**2 + cy**2)
    mask = (r_cells >= r1) & (r_cells <= r2)

    r_filtered = r_cells[mask]
    p_exact_filtered = pressure_exact[mask]

    vx_exact = velocity_exact[:, 0]
    vy_exact = velocity_exact[:, 1]
    speed_exact = np.sqrt(vx_exact**2 + vy_exact**2)
    s_exact_filtered = speed_exact[mask]

    order = np.argsort(r_filtered)
    r_sorted = r_filtered[order]
    p_exact_sorted = p_exact_filtered[order]
    s_exact_sorted = s_exact_filtered[order]

    p_exact_sorted = p_exact_sorted - p_exact_sorted[0]

    return r_sorted, s_exact_sorted, p_exact_sorted


def compute_radial_profiles(
    cx,
    cy,
    pressure,
    velocity,
    r1=1.0,
    r2=2.0,
    n_bins=80,
):
    """
    Compute radial averages for the numerical solution, binning to
    average out angular scatter at each radius.
    """

    vx = velocity[:, 0]
    vy = velocity[:, 1]

    speed = np.sqrt(vx**2 + vy**2)

    r_cells = np.sqrt(cx**2 + cy**2)

    mask = (r_cells >= r1) & (r_cells <= r2)

    r_filtered = r_cells[mask]
    p_filtered = pressure[mask]
    s_filtered = speed[mask]

    bins = np.linspace(r1, r2, n_bins + 1)

    bc = 0.5 * (bins[:-1] + bins[1:])

    p_bin = _bin_radially(r_filtered, p_filtered, bins, n_bins)
    s_bin = _bin_radially(r_filtered, s_filtered, bins, n_bins)

    if not np.isnan(p_bin).all():
        p_bin -= p_bin[~np.isnan(p_bin)][0]

    return bc, s_bin, p_bin


# =============================================================================
# Analytical Comparison Plot
# =============================================================================


def generate_analytical_plot(case, r_num, u_num, p_num, r_ref, u_ref, p_ref):
    """
    Generate analytical comparison plot. u_ref/p_ref are the raw
    per-cell analytical fields from the VTK output, sorted by radius
    r_ref (unbinned, so the line is smooth regardless of local mesh
    coarseness). u_num/p_num are the binned numerical solution at
    radii r_num.
    """

    fig, (ax1, ax2) = plt.subplots(
        1,
        2,
        figsize=(12, 5),
        dpi=DPI,
    )

    # Velocity
    ax1.plot(
        r_ref,
        u_ref,
        lw=2,
        color="black",
        label="Analytical",
    )

    ax1.plot(
        r_num,
        u_num,
        "o",
        ms=4,
        mfc="none",
        color="red",
        label="DuMux",
    )

    ax1.set_title("Velocity Profile")
    ax1.set_xlabel("r")
    ax1.set_ylabel("u_theta")
    ax1.legend()
    ax1.grid(True, alpha=0.3)

    # Pressure
    ax2.plot(
        r_ref,
        p_ref,
        lw=2,
        color="black",
        label="Analytical",
    )

    ax2.plot(
        r_num,
        p_num,
        "s",
        ms=4,
        mfc="none",
        color="tab:blue",
        label="DuMux",
    )

    ax2.set_title("Pressure Profile")
    ax2.set_xlabel("r")
    ax2.set_ylabel("p")
    ax2.legend()
    ax2.grid(True, alpha=0.3)

    plt.tight_layout()

    outpath = case["path"] / "analytical_comparison.png"

    fig.savefig(outpath, bbox_inches="tight")

    plt.close(fig)

    print(f"  Saved analytical comparison: {outpath}")


# =============================================================================
# ParaView Rendering
# =============================================================================


import matplotlib.tri as mtri


def render_matplotlib_fields(
    cx,
    cy,
    pressure,
    velocity,
    outdir: Path,
    r1=1.0,
    r2=2.0,
):
    """
    Generate pressure_field.png and velocity_field.png
    using pure matplotlib triangulation.
    """

    vx = velocity[:, 0]
    vy = velocity[:, 1]

    speed = np.sqrt(vx**2 + vy**2)

    triang = mtri.Triangulation(cx, cy)

    # Mask triangles inside inner cylinder
    tri_cx = cx[triang.triangles].mean(axis=1)
    tri_cy = cy[triang.triangles].mean(axis=1)

    tri_r = np.sqrt(tri_cx**2 + tri_cy**2)

    triang.set_mask(tri_r < r1)

    theta = np.linspace(0, 2*np.pi, 400)

    fields = [
        (
            pressure,
            "p",
            "coolwarm",
            "pressure_field.png",
        ),
        (
            speed,
            "|u|",
            "viridis",
            "velocity_field.png",
        ),
    ]

    for field_vals, label, cmap, filename in fields:

        fig, ax = plt.subplots(
            figsize=(7, 6),
            dpi=DPI,
        )

        tc = ax.tricontourf(
            triang,
            field_vals,
            levels=64,
            cmap=cmap,
        )

        fig.colorbar(
            tc,
            ax=ax,
            label=label,
            fraction=0.046,
            pad=0.04,
        )

        # Inner cylinder
        ax.plot(
            r1*np.cos(theta),
            r1*np.sin(theta),
            "k-",
            lw=1.2,
        )

        # Outer cylinder
        ax.plot(
            r2*np.cos(theta),
            r2*np.sin(theta),
            "k-",
            lw=1.2,
        )

        # White fill inside hole
        from matplotlib.patches import Circle

        ax.add_patch(
            Circle(
                (0, 0),
                r1,
                color="white",
                zorder=3,
            )
        )

        ax.plot(
            r1*np.cos(theta),
            r1*np.sin(theta),
            "k-",
            lw=1.2,
            zorder=4,
        )

        ax.set_aspect("equal")

        ax.set_xlabel("x")
        ax.set_ylabel("y")

        ax.set_title(label)

        ax.tick_params(
            which="both",
            direction="in",
        )

        plt.tight_layout()

        outpath = outdir / filename

        fig.savefig(
            outpath,
            bbox_inches="tight",
        )

        plt.close(fig)

        print(f"  Saved: {outpath}")

# =============================================================================
# Main
# =============================================================================


def main():

    parser = argparse.ArgumentParser(
        description="Unified Taylor-Couette post-processing suite"
    )

    parser.add_argument(
        "results_dir",
        nargs="?",
        default=DEFAULT_RESULTS_DIR,
        type=Path,
        help="Path to results directory",
    )

    parser.add_argument(
        "--fields",
        action="store_true",
        help="Generate pressure and velocity field maps",
    )

    args = parser.parse_args()

    results_dir = args.results_dir.resolve()

    if not results_dir.exists():
        raise RuntimeError(
            f"Results directory not found: {results_dir}"
        )

    case = find_case(results_dir)

    print("=" * 80)
    print(f"Results directory: {results_dir}")
    print("=" * 80)

    with tempfile.TemporaryDirectory() as tmp:

        tmpdir = Path(tmp)

        pvd_path = prepare_pvd_for_paraview(
            case["path"],
            tmpdir,
        )

        # -----------------------------------------------------------------
        # Extract data
        # -----------------------------------------------------------------

        cx, cy, pressure, velocity, pressure_exact, velocity_exact = parse_data_via_vtk(
            pvd_path
        )

        # -----------------------------------------------------------------
        # Radial profiles
        # -----------------------------------------------------------------

        r_num, u_num, p_num = compute_radial_profiles(
            cx,
            cy,
            pressure,
            velocity,
        )

        r_ref, u_ref, p_ref = sorted_exact_profile(
            cx,
            cy,
            pressure_exact,
            velocity_exact,
        )

        # -----------------------------------------------------------------
        # Analytical comparison
        # -----------------------------------------------------------------

        generate_analytical_plot(
            case,
            r_num,
            u_num,
            p_num,
            r_ref,
            u_ref,
            p_ref,
        )

        # -----------------------------------------------------------------
        # ParaView field rendering
        # -----------------------------------------------------------------

        if args.fields:

            print("  Rendering pressure field...")

            render_matplotlib_fields(
                cx,
                cy,
                pressure,
                velocity,
                case["path"],
            )


# =============================================================================
# Entry Point
# =============================================================================

if __name__ == "__main__":
    main()
