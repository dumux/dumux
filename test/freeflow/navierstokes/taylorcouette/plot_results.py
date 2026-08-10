#!/usr/bin/env pvpython
"""
plot_results.py
================

Unified post-processing suite for Taylor-Couette benchmark cases.

Features
--------
1. Global convergence plot
2. Optional ParaView 2D field rendering
3. Optional analytical comparison plots
4. Single-configuration filtering
5. Custom results directory support

Examples
--------
# Process ALL cases (convergence plot only)
pvpython3 plot_results.py

# Process a SPECIFIC configuration
pvpython3 plot_results.py --conf 100_0_80_320

# Also generate 2D field maps
pvpython3 plot_results.py --conf 100_0_80_320 --fields

# Also generate analytical comparison plot
pvpython3 plot_results.py --conf 100_0_80_320 --analytical

# Generate everything
pvpython3 plot_results.py --conf 100_0_80_320 --fields --analytical

# Explicit results directory
pvpython3 plot_results.py --conf 100_0_80_320 --fields --analytical /path/to/results
"""

from __future__ import annotations

import argparse
import re
import tempfile
import zipfile
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


def discover_cases(results_dir: Path, conf: str | None = None):
    """
    Discover benchmark case directories.
    """

    cases = []

    for d in sorted(results_dir.iterdir()):

        if not d.is_dir():
            continue

        if conf and d.name != conf:
            continue

        zip_file = d / "test_ff_taylorcouette.zip"

        if not zip_file.exists():
            continue

        tokens = d.name.split("_")

        if len(tokens) != 4:
            continue

        omega1, omega2, nr, ntheta = map(float, tokens)

        cases.append({
            "path": d,
            "zip": zip_file,
            "omega1": omega1,
            "omega2": omega2,
            "nr": int(nr),
            "ntheta": int(ntheta),
        })

    return cases


# =============================================================================
# Extraction Utilities
# =============================================================================


def extract_pvd_and_vtus(zip_path: Path, tmpdir: Path):
    """
    Extract PVD and VTU files from archive and clean paths inside PVD.
    """

    with zipfile.ZipFile(zip_path, "r") as zf:

        for item in zf.namelist():

            filename = Path(item).name

            if item.endswith(".pvd") or item.endswith(".vtu"):

                with zf.open(item) as source, open(tmpdir / filename, "wb") as target:
                    target.write(source.read())

    pvd_files = list(tmpdir.glob("*.pvd"))

    if not pvd_files:
        raise RuntimeError(f"No PVD file found in: {zip_path}")

    pvd_path = pvd_files[0]

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

    num_cells = client_data.GetNumberOfCells()

    pressure = np.array([
        vtk_p_array.GetValue(i)
        for i in range(num_cells)
    ])

    velocity = np.array([
        vtk_v_array.GetTuple(i)
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

    return x, y, pressure, velocity


# =============================================================================
# Analytical Solution
# =============================================================================


def analytical_solution(r, omega1, omega2, r1=1.0, r2=2.0):
    """
    Taylor-Couette analytical solution.
    """

    A = (
        (omega2 * r2**2 - omega1 * r1**2)
        / (r2**2 - r1**2)
    )

    B = (
        ((omega1 - omega2) * r1**2 * r2**2)
        / (r2**2 - r1**2)
    )

    u = A * r + B / r

    p = (
        A**2 * r**2 / 2
        + 2 * A * B * np.log(r)
        - B**2 / (2 * r**2)
    )

    p -= p[0]

    return u, p


# =============================================================================
# Radial Profiles
# =============================================================================


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
    Compute radial averages.
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

    idx = np.clip(
        np.digitize(r_filtered, bins) - 1,
        0,
        n_bins - 1,
    )

    p_bin = np.array([
        p_filtered[idx == i].mean()
        if (idx == i).any()
        else np.nan
        for i in range(n_bins)
    ])

    s_bin = np.array([
        s_filtered[idx == i].mean()
        if (idx == i).any()
        else np.nan
        for i in range(n_bins)
    ])

    if not np.isnan(p_bin).all():
        p_bin -= p_bin[~np.isnan(p_bin)][0]

    return bc, s_bin, p_bin


# =============================================================================
# Analytical Comparison Plot
# =============================================================================


def generate_analytical_plot(case, r_num, u_num, p_num):
    """
    Generate analytical comparison plot.
    """

    u_ref, p_ref = analytical_solution(
        r_num,
        case["omega1"],
        case["omega2"],
    )

    fig, (ax1, ax2) = plt.subplots(
        1,
        2,
        figsize=(12, 5),
        dpi=DPI,
    )

    # Velocity
    ax1.plot(
        r_num,
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
        r_num,
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
        "--conf",
        type=str,
        help="Process only a specific configuration",
    )

    parser.add_argument(
        "--fields",
        action="store_true",
        help="Generate pressure and velocity field maps",
    )

    parser.add_argument(
        "--analytical",
        action="store_true",
        help="Generate analytical comparison plots",
    )

    args = parser.parse_args()

    results_dir = args.results_dir.resolve()

    if not results_dir.exists():
        raise RuntimeError(
            f"Results directory not found: {results_dir}"
        )

    cases = discover_cases(results_dir, args.conf)

    if not cases:
        raise RuntimeError(
            f"No matching cases found inside: {results_dir}"
        )

    metrics_data = []

    print("=" * 80)
    print(f"Results directory: {results_dir}")

    if args.conf:
        print(f"Selected configuration: {args.conf}")

    print("=" * 80)

    # =========================================================================
    # Process Cases
    # =========================================================================

    for case in cases:

        print(f"\nProcessing Case: {case['path'].name}")

        with tempfile.TemporaryDirectory() as tmp:

            tmpdir = Path(tmp)

            pvd_path = extract_pvd_and_vtus(
                case["zip"],
                tmpdir,
            )

            # -----------------------------------------------------------------
            # Extract data
            # -----------------------------------------------------------------

            cx, cy, pressure, velocity = parse_data_via_vtk(
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

            # -----------------------------------------------------------------
            # Analytical comparison
            # -----------------------------------------------------------------

            if args.analytical:

                generate_analytical_plot(
                    case,
                    r_num,
                    u_num,
                    p_num,
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
