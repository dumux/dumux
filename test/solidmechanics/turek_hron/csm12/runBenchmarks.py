#!/usr/bin/env python3
# SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

"""
Stationary CSM Benchmark
Runs DuMuX simulations for CSM1 and CSM2 and compares results against
reference data
"""

import os
import subprocess
import pandas as pd
import matplotlib.pyplot as plt

b_names = ["csm1", "csm2"]
youngs_moduli = [1.4e6, 5.6e6]

scriptDirectory = os.path.dirname(os.path.abspath(__file__))
referenceDirectory = os.path.normpath(os.path.join(scriptDirectory, "..", "reference_data"))

###############################################################################
# DuMux Simulation and Data Management
###############################################################################

def cleanup_output(case=0, includeResults=False):
    """Finds and deletes benchmark-specific tracking files."""
    extensions = [".vtu", ".pvd"]
    if includeResults:
        extensions.append(".csv")

    for filename in os.listdir("."):
        if any(filename.startswith(b_names[case]) and filename.endswith(extensions)):
            try:
                os.remove(filename)
            except OSError:
                pass

def run_simulation(case=0, refinement=0):
    """Executes the DuMuX binary for stationary CSM cases."""
    binary = "./test_saint_venant_kirchhoff_stationary"
    if not os.path.exists(binary):
        print(f"Error: Binary {binary} not found. Please compile the test.")
        return False
    cmd = [
        binary, "params.input",
        "-Problem.Name", b_names[case],
        "-SpatialParams.YoungsModulus", str(youngs_moduli[case]),
        "-Grid.Refinement", str(refinement)
    ]
    subprocess.check_call(cmd)
    return True

def get_result(case=0):
    """Loads and filters DuMuX simulation results from CSV."""
    path = f"{b_names[case]}_dumux_box.csv"
    if not os.path.exists(path):
        return None
    df = pd.read_csv(path)
    df[["ux", "uy"]] *= 1000.0 # m to mm
    return df

###############################################################################
# Reference Data Access
###############################################################################

def get_reference(case=0):
    """Loads and filters external reference data for a specific case."""
    referenceFile = os.path.join(referenceDirectory, f"{b_names[case]}.csv")
    if not os.path.exists(referenceFile):
        return None
    df = pd.read_csv(referenceFile)
    return df[df["level"].str.endswith("+0")].copy()

###############################################################################
# Visualizing results
###############################################################################

def print_summary_table(case, output_md=True):
    """Formats and prints a summary table of the benchmark results."""
    df = get_result(case)

    if df is not None and not df.empty:
        summary_df = df[["level", "nel", "ndof", "ux", "uy"]].copy()
        summary_df.insert(0, "case", b_names[case])
        summary_df.columns = ["case","level","nel","ndof",
                              "ux [mm]","uy [mm]"]

        print(f"\n Dumux Stationary CSM Benchmark Results ({b_names[case]}):")
        formatters = {"ux [mm]": "{:.5f}".format, "uy [mm]": "{:.5f}".format}
        print(summary_df.to_string(index=False, formatters=formatters))

        if output_md:
            md_filename = f"csm{case + 1}.md"
            with open(md_filename, "w") as f:
                f.write(f"# Dumux Stationary CSM Benchmark Results: {b_names[case].upper()}\n\n")
                f.write(summary_df.to_markdown(index=False, floatfmt=".5f"))
                f.write("\n")

def plot_comparison():
    """Generates a 2x2 grid comparing DuMux results with reference data."""
    fig, axes = plt.subplots(2, 2, figsize=(12, 10), sharex=True)
    components = ["ux", "uy"]
    labels = {"ux": "$u_x$ at A [mm]", "uy": "$u_y$ at A [mm]"}

    for col_idx in [0, 1]:
        df = get_result(col_idx)
        ref_df = get_reference(col_idx)
        if df is None:
            continue

        # Keep all levels to see full convergence path
        # df = df[df["level"] > 0]

        for row_idx, comp in enumerate(components):
            ax = axes[row_idx, col_idx]

            # PLOT: Use standard linear scale for Y to visually inspect convergence
            ax.plot(df["ndof"], df[comp], "o-", label="DuMuX", color="C0", markersize=6)
            if ref_df is not None:
                ax.plot(ref_df["ndof"], ref_df[comp],
                        "s--", label="Reference (TU Dortmund)", color="C1", markersize=5)

            ax.set_xscale("log")
            ax.set_yscale("linear") # Crucial for seeing physical convergence
            ax.grid(True, which="both", ls="-", alpha=0.3)

    # Polishing titles and labels
    for col_idx, case_name in enumerate(b_names):
        e_val = youngs_moduli[col_idx]
        axes[0, col_idx].set_title(
            f"{case_name.upper()} (E = {e_val:.1e} Pa)",
            fontsize=12,
            fontweight="bold",
            pad=15,
        )
        axes[1, col_idx].set_xlabel("Degrees of Freedom (ndof)", fontsize=11)

    # Dynamic Y-labels using clean LaTeX rendering
    axes[0, 0].set_ylabel(labels["ux"], fontsize=12, fontweight="bold")
    axes[1, 0].set_ylabel(labels["uy"], fontsize=12, fontweight="bold")

    handles, labels_legend = axes[0, 0].get_legend_handles_labels()
    if handles:
        fig.legend(
            handles, labels_legend, loc="upper center", bbox_to_anchor=(0.5, 0.96),
            ncol=2, frameon=True, fontsize=11
        )
        fig.tight_layout()
        # Adjusted to make room for top legend safely
        fig.subplots_adjust(top=0.88)

    output_plot = "csm12_comparison.png"
    plt.savefig(output_plot, dpi=300)
    plt.close(fig)

###############################################################################
# Actual main
###############################################################################
for i in [0, 1]:
    #for level in [0, 1, 2, 3, 4]:
        #run_simulation(i, level)
    print_summary_table(i)
plot_comparison()
