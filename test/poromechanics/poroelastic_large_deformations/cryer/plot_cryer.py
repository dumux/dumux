#!/usr/bin/env python3
# SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later
"""
Plot the Cryer sphere consolidation benchmark results.

Produces, for the Kozeny-Carman and constant permeability models of the large-deformation model:
  - normalised pore pressure at the centre vs time
  - volume ratio J = V/V0 at the centre vs time
plus a single-run comparison. The centre pressure of the small-deformation model is added to the
constant-permeability pressure plots. Cryer's solution is taken from the CSV files.

Usage:
  python3 plot_cryer.py [--csv cryer_center_pressure.csv] [--outdir .]

The CSV file is written by main.cc; it contains columns:
  time, t*cv/R0^2, p_f_center, p_analytical, p_f/p0, p_analytical/p0, J_center
"""
import argparse
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import rcParams

rcParams["font.size"] = 10
rcParams["axes.grid"] = True
rcParams["grid.linestyle"] = ":"
rcParams["grid.color"] = "0.7"

LABELS = {
    "p0k_0.0001": r"$p_0/K=0.0001$",
    "p0k_0.25": r"$p_0/K=0.25$",
    "p0k_0.5": r"$p_0/K=0.5$",
}
MARKERS = {"p0k_0.0001": "*", "p0k_0.25": "o", "p0k_0.5": "^"}
COLORS = {"p0k_0.0001": "k", "p0k_0.25": "g", "p0k_0.5": "r"}
SMALL_DEFORMATION_STYLE = {"color": "tab:purple", "linestyle": "-.", "marker": "s",
                           "label": "linear poroelasticity"}


def load_csv(path):
    """Load time-series CSV from main.cc output."""
    data = {}
    with open(path) as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.strip().split(",")
            if len(parts) < 6:
                continue
            vals = [float(p) for p in parts]
            t, tNorm, pfCenter, pAnal, pfRel, pAnalRel = vals[:6]
            J = vals[6] if len(vals) >= 7 else float("nan")
            data.setdefault("t", []).append(t)
            data.setdefault("tNorm", []).append(tNorm)
            data.setdefault("pf", []).append(pfCenter)
            data.setdefault("pAnal", []).append(pAnal)
            data.setdefault("pfRel", []).append(pfRel)
            data.setdefault("pAnalRel", []).append(pAnalRel)
            data.setdefault("J", []).append(J)
    return {k: np.array(v) for k, v in data.items()}


def plot_pressure_figure(datasets, xlabel, ylabel, outfile, small=None):
    """Plot one normalised centre-pressure figure for the three loads."""
    fig, ax = plt.subplots(figsize=(7, 5))
    if small is not None:
        ax.semilogx(small["tNorm"], small["pfRel"], markevery=max(1, len(small["tNorm"]) // 40),
                    linewidth=1, **SMALL_DEFORMATION_STYLE)

    ymax = 1.6
    for key, d in datasets.items():
        tN = d["tNorm"]
        every = max(1, len(tN) // 40)
        ax.semilogx(tN, d["pfRel"],
                    linestyle="-", color=COLORS[key], marker=MARKERS[key],
                    markevery=every, label=LABELS[key], linewidth=1)
        ymax = max(ymax, float(np.nanmax(d["pfRel"])) * 1.05)

    reference = next(iter(datasets.values()))
    ax.semilogx(reference["tNorm"], reference["pAnalRel"], "r-", linewidth=2,
                label="Analytical solution")

    ax.set_xlim([1e-4, 1e0])
    ax.set_ylim([0.0, ymax])
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.legend(loc="center right")
    fig.tight_layout()
    fig.savefig(outfile, dpi=150)
    plt.close(fig)
    print(f"Saved {outfile}")


def plot_volume_figure(datasets, xlabel, ylabel, outfile):
    """Plot one volume-ratio figure: J = V/V0 at the centre for the three loads."""
    fig, ax = plt.subplots(figsize=(7, 5))

    plotted = False
    for key, d in datasets.items():
        if "J" not in d or np.all(np.isnan(d["J"])):
            continue
        tN = d["tNorm"]
        every = max(1, len(tN) // 40)
        ax.semilogx(tN, d["J"],
                    linestyle="-", color=COLORS[key], marker=MARKERS[key],
                    markevery=every, label=LABELS[key], linewidth=1)
        plotted = True

    if not plotted:
        plt.close(fig)
        print(f"  No J data available, skipping {outfile}")
        return

    ax.set_xlim([1e-4, 1e0])
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.legend(loc="lower left")
    fig.tight_layout()
    fig.savefig(outfile, dpi=150)
    plt.close(fig)
    print(f"Saved {outfile}")


def plot_overview(datasets_by_perm, small, outfile, dpi=80):
    """Centre pressure and volume ratio for both permeability models in one figure."""
    colors = {"p0k_0.0001": "tab:blue", "p0k_0.25": "tab:green", "p0k_0.5": "tab:orange"}
    titles = {"Constant": "constant permeability", "KozenyCarman": "Kozeny-Carman permeability"}
    fig, axes = plt.subplots(2, 2, figsize=(8, 5.5), sharex=True)
    for col, perm in enumerate(["Constant", "KozenyCarman"]):
        axP, axJ = axes[0][col], axes[1][col]
        reference = next(iter(datasets_by_perm[perm].values()))
        axP.semilogx(reference["tNorm"], reference["pAnalRel"], "k--", linewidth=1,
                     label="Cryer (1963)", zorder=3)
        if perm == "Constant" and small is not None:
            axP.semilogx(small["tNorm"], small["pfRel"], markersize=3, linewidth=1,
                         **SMALL_DEFORMATION_STYLE)
        for key, d in datasets_by_perm[perm].items():
            axP.semilogx(d["tNorm"], d["pfRel"], color=colors[key], marker=MARKERS[key],
                         markersize=4, linewidth=1, label=LABELS[key])
            axJ.semilogx(d["tNorm"], d["J"], color=colors[key], marker=MARKERS[key],
                         markersize=4, linewidth=1, label=LABELS[key])
        axP.set_title(titles[perm], fontsize=10)
        axP.set_ylim(0.0, 2.0)
        axJ.set_ylim(0.4, 1.05)
        axJ.set_xlim(1e-3, 1.0)
        axJ.set_xlabel(r"$t \cdot c_v / R_0^2$")
        if col == 0:
            axP.set_ylabel(r"$p / p_0$ at centre")
            axJ.set_ylabel(r"$J$ at centre")
    axes[0][0].legend(fontsize=8, loc="lower left")
    fig.tight_layout()
    fig.savefig(outfile, dpi=dpi)
    plt.close(fig)
    print(f"Saved {outfile}")


def plot_simulation_vs_analytical(sim_csv, outfile, label="simulation"):
    """Plot a single simulation CSV (from main.cc output) vs the analytical solution."""
    data = load_csv(sim_csv)
    if not data:
        print(f"  Could not load {sim_csv}")
        return

    tN = data["tNorm"]
    pRel = data["pfRel"]

    fig, ax = plt.subplots(figsize=(7, 5))
    ax.semilogx(tN, pRel, "k-o", markevery=1, markersize=5, label=label)
    ax.semilogx(tN, data["pAnalRel"], "r-", linewidth=2, label="Analytical solution (Cryer 1963)")
    ax.axhline(1.0, color="0.5", linestyle="--", linewidth=0.8)
    ax.set_xlim([max(1e-4, tN.min() * 0.5), max(1e0, tN.max() * 2)])
    ax.set_ylim([0.0, max(1.65, pRel.max() * 1.1)])
    ax.set_xlabel(r"$t \cdot c_v / R_0^2$")
    ax.set_ylabel(r"$p(0,t) / p_0$")
    ax.legend()
    fig.tight_layout()
    fig.savefig(outfile, dpi=150)
    plt.close(fig)
    print(f"Saved {outfile}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--csv", default="cryer_center_pressure.csv",
                    help="Direct CSV written by the simulation (single run)")
    ap.add_argument("--csv_prefix", default="cryer_center_pressure",
                    help="Prefix for parametric CSV files; expects <prefix>_<tag>.csv")
    ap.add_argument("--outdir", default=".")
    args = ap.parse_args()

    csv_path = args.csv
    if not os.path.isabs(csv_path):
        csv_path = os.path.join(args.outdir, csv_path)
    if os.path.exists(csv_path):
        plot_simulation_vs_analytical(
            csv_path,
            outfile=os.path.join(args.outdir, "fig_simulation_vs_analytical.png"),
            label="DuMux"
        )
    elif os.path.exists(args.csv):
        plot_simulation_vs_analytical(
            args.csv,
            outfile=os.path.join(args.outdir, "fig_simulation_vs_analytical.png"),
            label="DuMux"
        )
    else:
        print(f"Single-run CSV not found: {args.csv}")

    p0k_values = [0.0001, 0.25, 0.5]
    perm_models = ["KozenyCarman", "Constant"]
    fig_files = {
        "KozenyCarman": "pressure_kozenycarman.png",
        "Constant": "pressure_constant.png",
    }
    fig_files_vol = {
        "KozenyCarman": "volume_kozenycarman.png",
        "Constant": "volume_constant.png",
    }

    small = None
    for fname in (f"{args.csv_prefix}_small_deformation.csv",
                  os.path.join(args.outdir, f"{args.csv_prefix}_small_deformation.csv")):
        if os.path.exists(fname):
            small = load_csv(fname)
            break

    all_datasets = {}
    for perm in perm_models:
        datasets = {}
        for p0k in p0k_values:
            key = f"p0k_{p0k}"
            fname = f"{args.csv_prefix}_{perm}_p0k{p0k}.csv"
            if not os.path.exists(fname):
                fname = os.path.join(args.outdir, fname)
            if not os.path.exists(fname):
                continue
            datasets[key] = load_csv(fname)

        if not datasets:
            continue
        all_datasets[perm] = datasets

        plot_pressure_figure(
            datasets,
            xlabel=r"$t \cdot c_v / R_0^2$",
            ylabel=r"$p / p_0$",
            outfile=os.path.join(args.outdir, fig_files[perm]),
            small=small if perm == "Constant" else None,
        )

        plot_volume_figure(
            datasets,
            xlabel=r"$t \cdot c_v / R_0^2$",
            ylabel=r"$J = V/V_0$ at center",
            outfile=os.path.join(args.outdir, fig_files_vol[perm]),
        )

    if all(all_datasets.get(perm) for perm in perm_models):
        plot_overview(all_datasets, small, os.path.join(args.outdir, "cryer_overview.png"))

    print("Done.")


if __name__ == "__main__":
    main()
