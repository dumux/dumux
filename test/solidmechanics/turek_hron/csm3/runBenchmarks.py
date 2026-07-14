#!/usr/bin/env python3
# SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

"""
Instationary CSM Benchmark
Runs DuMuX simulations for CSM3 and compares results against
reference data
"""

import os
import subprocess
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

script_directory = os.path.dirname(os.path.abspath(__file__))
reference_diretory = os.path.normpath(os.path.join(script_directory, "..", "reference_data"))

def compute_csm3_metrics(time, signal, dt_nominal):
    """
    Compute mean, amplitude, and dominant frequency from a transient signal
    during its fully-stabilized limit cycle phase (t >= 8.0s).
    """
    t_arr = np.array(time)
    sig_arr = np.array(signal)

    mask = t_arr >= 8.0
    t_stable = t_arr[mask]
    sig_stable = sig_arr[mask]

    if len(t_stable) < 10:
        return 0.0, 0.0, 0.0

    val_max = np.max(sig_stable)
    val_min = np.min(sig_stable)
    mean_val = (val_max + val_min) / 2.0
    amplitude = (val_max - val_min) / 2.0

    zero_crossings = np.where(np.diff(np.sign(sig_stable - mean_val)))[0]
    if len(zero_crossings) >= 3:
        t_start = t_stable[zero_crossings[0]]
        t_end = t_stable[zero_crossings[-1]]
        num_half_cycles = len(zero_crossings) - 1
        frequency = (num_half_cycles / 2.0) / (t_end - t_start)
    else:
        frequency = 1.0 / (dt_nominal * 100.0)

    return mean_val, amplitude, frequency

###############################################################################
# DuMux Simulation and Data Management
###############################################################################

def run_simulation(name, refinement, dt, vtk):
    """Executes the DuMuX binary for CSM3."""
    binary = "./test_saint_venant_kirchhoff_dynamic"
    if not os.path.exists(binary):
        print(f"Error: Binary {binary} not found.")
        return False
    cmd = [binary, "params.input",
           "-Problem.Name", name,
           "-Grid.Refinement", str(refinement),
           "-TimeLoop.DtInitial", str(dt),
           "-VTKOutput.Every", str(vtk)]
    subprocess.check_call(cmd)
    return True

def get_result(name):
    """Loads DuMuX simulation results from CSV."""
    path = f"{name}_dumux_box.csv"
    if not os.path.exists(path):
        return None
    df = pd.read_csv(path)
    df[["ux", "uy"]] *= 1000.0 # m to mm
    return df

###############################################################################
# Reference Data Access
###############################################################################

def get_reference(level, dt):
    dt_str = str(dt).replace(".", "p")
    filename = f"csm3_l{level}_t{dt_str}.csv"
    path = os.path.join(reference_diretory, filename)

    if not os.path.exists(path):
        return None

    return pd.read_csv(path)

def get_reference_summary():
    """Loads reference data"""
    ref_file = os.path.join(reference_diretory, "csm3.csv")
    if not os.path.exists(ref_file):
        return None
    df = pd.read_csv(ref_file)
    return df

###############################################################################
# Visualizing results
###############################################################################

def print_summary_table(results, output_md=True):
    """Formats and prints a clean comparison table sorted primarily by timestep."""
    table_rows = []
    for entry in results:
        df = entry["data"]
        if df is not None:
            dt_val = entry["dt"]
            nel = int(df["nel"].iloc[0]) if hasattr(df["nel"], "iloc") else int(df["nel"])
            ndof = int(df["ndof"].iloc[0]) if hasattr(df["ndof"], "iloc") else int(df["ndof"])

            ux_mean, ux_amp, ux_freq = compute_csm3_metrics(df["t"], df["ux"], dt_val)
            uy_mean, uy_amp, uy_freq = compute_csm3_metrics(df["t"], df["uy"], dt_val)

            table_rows.append({
                "dt": dt_val,
                "nel": nel,
                "ndofs": ndof,
                "Source": "DuMuX",
                "Discretisation": "Box",
                "ux of A [mm] (Mean +/- Amp [Freq])": f"{ux_mean:+.3f} +/- {ux_amp:.3f} [{ux_freq:.4f}]",
                "uy of A [mm] (Mean +/- Amp [Freq])": f"{uy_mean:+.3f} +/- {uy_amp:.3f} [{uy_freq:.4f}]"
            })

    reference_table_data = get_reference_summary()
    if reference_table_data is not None:
        for _, ref_row in reference_table_data.iterrows():
            table_rows.append({
                "dt": float(ref_row["dt"]),
                "nel": int(ref_row["nel"]),
                "ndofs": int(ref_row["ndof"]),
                "Source": "FeatFlow",
                "Discretisation": "FEM",
                "ux of A [mm] (Mean +/- Amp [Freq])": f"{ref_row['ux_mean']:+.3f} +/- {ref_row['ux_error']:.3f} [{ref_row['ux_frequency']:.4f}]",
                "uy of A [mm] (Mean +/- Amp [Freq])": f"{ref_row['uy_mean']:+.3f} +/- {ref_row['uy_error']:.3f} [{ref_row['uy_frequency']:.4f}]"
            })

    if table_rows:
        df_compare = pd.DataFrame(table_rows)
        df_compare = df_compare.sort_values(by=["dt", "nel", "Source"], ascending=[True, True, True])

        print("\n=== CSM3 Benchmark Table ===")
        print(df_compare.to_string(index=False))

        if output_md:
            md_lines = [
                "| dt | nel | ndofs | Source | Discretisation | ux of A [mm] (Mean +/- Amp [Freq]) | uy of A [mm] (Mean +/- Amp [Freq]) |",
                "| :--- | :--- | :--- | :--- | :--- | :--- |"
            ]
            for _, row in df_compare.iterrows():
                md_lines.append(f"| {row['dt']:.3f} | {row['nel']} | {row['ndofs']} | {row['Source']} | {row['Discretisation']} | {row['ux of A [mm] (Mean +/- Amp [Freq])']} | {row['uy of A [mm] (Mean +/- Amp [Freq])']} |")

            with open("csm_3.md", "w") as f:
                f.write("\n".join(md_lines) + "\n")

def generate_plots(dt, level_d, level_ref, output_filename="csm3_comparison.png"):
    """Compare results over time with a unified global legend and fine reference."""

    ref_fine = get_reference(4, 0.005)
    local_results = get_result(f"csm3_refinement_{level_d}_dt_{dt}")
    reference_results = get_reference(level_ref, dt)

    t_ref, ux_ref, uy_ref = reference_results["t"], reference_results["ux"], reference_results["uy"]
    t_loc, ux_loc, uy_loc = local_results["t"], local_results["ux"], local_results["uy"]
    t_fine, ux_fine, uy_fine = ref_fine["t"], ref_fine["ux"], ref_fine["uy"]

    fig, axs = plt.subplots(2, 2, figsize=(14, 8))
    fig.suptitle("CSM3 Comparison", fontsize=14, fontweight='bold')

    lbl_dumux = f"DuMuX (l={level_d}, $\\Delta t$ = {dt}s)"
    lbl_ref = f"Reference (l={level_ref}, $\\Delta t$ = {dt}s)"
    lbl_fine = "Reference (l=4, $\\Delta t$ = 0.005s)"

    axs[0,0].plot(t_fine, ux_fine, 'g-.', alpha=0.5)
    axs[0,0].plot(t_ref, ux_ref, 'r-', alpha=0.7)
    axs[0,0].plot(t_loc, ux_loc, 'k--')
    axs[0,0].set_ylabel("$u_x$ [mm]")

    l_fine, = axs[0,1].plot(t_fine, uy_fine, 'g-.', alpha=0.5)
    l_ref, = axs[0,1].plot(t_ref, uy_ref, 'b-', alpha=0.7)
    l_loc, = axs[0,1].plot(t_loc, uy_loc, 'k--')
    axs[0,1].set_ylabel("$u_y$ [mm]")

    zoom_configs = [
        (axs[1,0], ux_fine, ux_ref, ux_loc, 'r-', "$u_x$ [mm] (Zoomed)"),
        (axs[1,1], uy_fine, uy_ref, uy_loc, 'b-', "$u_y$ [mm] (Zoomed)")
    ]

    for ax, data_fine, data_ref, data_loc, ref_style, ylabel in zoom_configs:
        ax.plot(t_fine, data_fine, 'g-.', alpha=0.5)
        ax.plot(t_ref, data_ref, ref_style, alpha=0.7)
        ax.plot(t_loc, data_loc, 'k--')
        ax.set_xlim(8.0, 10.0)
        ax.set_ylabel(ylabel)

    for ax in axs.flat:
        ax.set_xlabel("Time [s]")
        ax.grid(True, linestyle=":", alpha=0.6)

    plt.tight_layout()

    fig.legend(
        [l_loc, l_ref, l_fine],
        [lbl_dumux, lbl_ref, lbl_fine],
        loc="center left",
        bbox_to_anchor=(0.78, 0.5),
        fontsize=11,
        frameon=True
    )

    plt.savefig(output_filename, dpi=150)
    plt.close(fig)


###############################################################################
# Actual main
###############################################################################

base_name = "csm3"
levels = [0, 1, 2, 3]
dt = [0.02, 0.01]
vtk_out = [10, 20]

all_results = []
for i in range(len(dt)):
    for j in range(len(levels)):
        name = f"{base_name}_refinement_{levels[j]}_dt_{dt[i]}"
        #run_simulation(name, levels[j], dt[i], vtk_out[i])
        df = get_result(name)
        all_results.append({"level": levels[j], "dt": dt[i], "data": df})

print_summary_table(all_results)

generate_plots(0.01, 2, 2)
