#!/usr/bin/env python3
"""
Plot the Cryer sphere consolidation benchmark results.

Produces, for the Kozeny-Carman and constant permeability models:
  - normalised pore pressure at the centre vs time
  - volume ratio J = V/V0 at the centre vs time
plus the analytical reference and a single-run comparison.

Usage:
  python3 plot_cryer.py [--csv cryer_center_pressure.csv] [--outdir .]

The CSV file is written by main.cc; it contains columns:
  time, t*cv/R0^2, p_f_center, p_analytical, p_f/p0, p_analytical/p0, J_center
"""
import argparse
import os
import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import rcParams

# Figure style
rcParams["font.size"] = 10
rcParams["axes.grid"] = True
rcParams["grid.linestyle"] = ":"
rcParams["grid.color"] = "0.7"

LABELS = {
    "p0k_0.0001": r"$p_0/K=0.0001$",
    "p0k_0.25":   r"$p_0/K=0.25$",
    "p0k_0.5":    r"$p_0/K=0.5$",
}
MARKERS = {"p0k_0.0001": "*", "p0k_0.25": "o", "p0k_0.5": "^"}
COLORS  = {"p0k_0.0001": "k", "p0k_0.25":  "g", "p0k_0.5":  "r"}


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
            J = vals[6] if len(vals) >= 7 else float("nan")  # volume ratio (optional)
            data.setdefault("t", []).append(t)
            data.setdefault("tNorm", []).append(tNorm)
            data.setdefault("pf", []).append(pfCenter)
            data.setdefault("pAnal", []).append(pAnal)
            data.setdefault("pfRel", []).append(pfRel)
            data.setdefault("pAnalRel", []).append(pAnalRel)
            data.setdefault("J", []).append(J)
    return {k: np.array(v) for k, v in data.items()}


def plot_pressure_figure(datasets, analytical_tNorm, analytical_pRel,
                         xlabel, ylabel, outfile):
    """Plot one normalised centre-pressure figure for the three loads."""
    fig, ax = plt.subplots(figsize=(7, 5))

    ymax = 1.6
    for key, d in datasets.items():
        tN = d["tNorm"]
        # thin out markers (show ~10 per decade)
        every = max(1, len(tN) // 40)
        ax.semilogx(tN, d["pfRel"],
                    linestyle="-", color=COLORS[key], marker=MARKERS[key],
                    markevery=every, label=LABELS[key], linewidth=1)
        ymax = max(ymax, float(np.nanmax(d["pfRel"])) * 1.05)

    # Analytical solution (red line, no markers)
    if analytical_tNorm is not None:
        ax.semilogx(analytical_tNorm, analytical_pRel, "r-", linewidth=2,
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


def compute_analytical_curve(K, G, alphaB, Sp, kappa0, muf, R0, nRoots=120):
    """Compute analytical Cryer solution p(0,t)/p0 over normalised time."""
    m    = (K + 4.0/3.0*G) / (4.0*G) * (1.0 + K*Sp / (alphaB**2))
    cv   = kappa0 * (K + 4.0/3.0*G) / muf
    beta = alphaB**2 + (K + 4.0/3.0*G)*Sp

    # Find roots of (1 - m*z^2)*sin(z) - z*cos(z) = 0
    def f(z):
        return (1.0 - m*z**2)*np.sin(z) - z*np.cos(z)

    roots = []
    root_tol = 1e-6
    n_sub = 128
    for k in range(400):
        interval_lo, interval_hi = k*np.pi + 1e-9, (k+1)*np.pi - 1e-9
        step = (interval_hi - interval_lo) / n_sub

        lo = interval_lo
        flo = f(lo)
        for s in range(n_sub):
            hi = interval_hi if s == n_sub - 1 else lo + step
            fhi = f(hi)

            if flo*fhi > 0:
                lo = hi
                flo = fhi
                continue

            a, b = lo, hi
            fa = flo
            for _ in range(200):
                mid = 0.5*(a + b)
                fm = f(mid)
                if abs(fm) < 1e-14 or (b - a) < 1e-14:
                    a = mid
                    b = mid
                    break
                if fa*fm <= 0:
                    b = mid
                else:
                    a = mid
                    fa = fm

            root = 0.5*(a + b)
            if root > root_tol and (not roots or abs(root - roots[-1]) > root_tol):
                roots.append(root)
                if len(roots) >= nRoots:
                    break

            lo = hi
            flo = fhi

        if len(roots) >= nRoots:
            break

    # Normalised time axis: t*cv/R0^2 from 1e-4 to 1e1
    tNorm = np.logspace(-4, 1, 300)
    pRel  = np.zeros_like(tNorm)
    for i, tn in enumerate(tNorm):
        s = 0.0
        for z in roots:
            num = np.sin(z) - z
            den = m*z*np.cos(z) + (2.0*m - 1.0)*np.sin(z)
            if abs(den) < 1e-14:
                continue
            exp_arg = -z**2 * tn / beta   # tn = cv*t/R0^2
            if exp_arg < -700:
                continue
            s += num/den * np.exp(exp_arg)
        pRel[i] = 2.0*m*s

    return tNorm, pRel


def plot_analytical_only(aT, aP, outfile, K, G, alphaB, Sp, kappa0, muf, R0):
    """Standalone plot of the analytical Cryer solution."""
    fig, ax = plt.subplots(figsize=(7, 5))
    ax.semilogx(aT, aP, "r-", linewidth=2, label="Analytical solution (Cryer 1963)")
    ax.axhline(1.0, color="0.5", linestyle="--", linewidth=0.8)
    ax.set_xlim([1e-4, 1e0])
    ax.set_ylim([0.0, 1.65])
    ax.set_xlabel(r"$t \cdot c_v / R_0^2$")
    ax.set_ylabel(r"$p(0,t) / p_0$")
    ax.legend()
    fig.tight_layout()
    fig.savefig(outfile, dpi=150)
    plt.close(fig)
    print(f"Saved {outfile}")


def plot_simulation_vs_analytical(sim_csv, aT, aP, outfile, label="simulation"):
    """Plot a single simulation CSV (from main.cc output) vs the analytical solution."""
    data = load_csv(sim_csv)
    if not data:
        print(f"  Could not load {sim_csv}")
        return

    tN   = np.array(data["tNorm"])
    pRel = np.array(data["pfRel"])

    fig, ax = plt.subplots(figsize=(7, 5))
    ax.semilogx(tN, pRel, "k-o", markevery=1, markersize=5, label=label)
    ax.semilogx(aT, aP,   "r-",  linewidth=2, label="Analytical solution (Cryer 1963)")
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
    ap.add_argument("--csv",        default="cryer_center_pressure.csv",
                    help="Direct CSV written by the simulation (single run)")
    ap.add_argument("--csv_prefix", default="cryer_center_pressure",
                    help="Prefix for parametric CSV files; expects <prefix>_<tag>.csv")
    ap.add_argument("--outdir", default=".")
    # Material parameters (must match the simulation)
    ap.add_argument("--K",      type=float, default=1e6)
    ap.add_argument("--G",      type=float, default=1.5e6)
    ap.add_argument("--alphaB", type=float, default=1.0)
    ap.add_argument("--Sp",     type=float, default=0.0)
    ap.add_argument("--kappa0", type=float, default=1e-12)
    ap.add_argument("--muf",    type=float, default=1e-3)
    ap.add_argument("--R0",     type=float, default=1.0)
    args = ap.parse_args()

    cv = args.kappa0 * (args.K + 4.0/3.0*args.G) / args.muf
    print(f"Parameters: K={args.K:.3g} Pa, G={args.G:.3g} Pa, cv={cv:.3g} m²/s")
    m = (args.K + 4.0/3.0*args.G) / (4.0*args.G) * (1.0 + args.K*args.Sp/args.alphaB**2)
    print(f"  m={m:.4f},  t_end = R0²/cv = {args.R0**2/cv:.2f} s")

    # Analytical curve
    aT, aP = compute_analytical_curve(args.K, args.G, args.alphaB, args.Sp,
                                       args.kappa0, args.muf, args.R0)

    # Standalone analytical plot
    plot_analytical_only(aT, aP,
                          os.path.join(args.outdir, "fig_analytical.png"),
                          args.K, args.G, args.alphaB, args.Sp,
                          args.kappa0, args.muf, args.R0)

    # Single-run comparison plot
    csv_path = args.csv
    if not os.path.isabs(csv_path):
        csv_path = os.path.join(args.outdir, csv_path)
    if os.path.exists(csv_path):
        plot_simulation_vs_analytical(
            csv_path, aT, aP,
            outfile=os.path.join(args.outdir, "fig_simulation_vs_analytical.png"),
            label="DuMux"
        )
    elif os.path.exists(args.csv):
        plot_simulation_vs_analytical(
            args.csv, aT, aP,
            outfile=os.path.join(args.outdir, "fig_simulation_vs_analytical.png"),
            label="DuMux"
        )
    else:
        print(f"Single-run CSV not found: {args.csv}")

    # Parametric comparison (multiple loads / permeability models)
    p0k_values = [0.0001, 0.25, 0.5]
    perm_models = ["KozenyCarman", "Constant"]
    fig_files = {
        "KozenyCarman": "pressure_kozenycarman.png",
        "Constant":     "pressure_constant.png",
    }
    fig_files_vol = {
        "KozenyCarman": "volume_kozenycarman.png",
        "Constant":     "volume_constant.png",
    }

    for perm in perm_models:
        datasets = {}
        for p0k in p0k_values:
            key  = f"p0k_{p0k}"
            fname = f"{args.csv_prefix}_{perm}_p0k{p0k}.csv"
            if not os.path.exists(fname):
                fname = os.path.join(args.outdir, fname)
            if not os.path.exists(fname):
                continue
            datasets[key] = load_csv(fname)

        if not datasets:
            continue

        # normalised centre pressure vs analytical
        plot_pressure_figure(
            datasets, aT, aP,
            xlabel=r"$t \cdot c_v / R_0^2$",
            ylabel=r"$p / p_0$",
            outfile=os.path.join(args.outdir, fig_files[perm]),
        )

        # volume ratio J = V/V_0 at the centre
        plot_volume_figure(
            datasets,
            xlabel=r"$t \cdot c_v / R_0^2$",
            ylabel=r"$J = V/V_0$ at center",
            outfile=os.path.join(args.outdir, fig_files_vol[perm]),
        )

    print("Done.")


if __name__ == "__main__":
    main()
