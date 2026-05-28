#!/usr/bin/env python3
"""
Spatial convergence study for the Boussinesq one-sided Rayleigh–Bénard test.

Doubles the grid in each direction starting from BASE_NX x BASE_NY, collects
the Sherwood-number CSV written by each run, and produces a two-panel log–log
figure:
  left  – Sh(t) curves for every resolution
  right – L² error of Sh(t) vs. finest-grid reference as a function of h = H/N_y

Usage
-----
    python3 convergence.py <path/to/test_boussinesq_onesided_rb_box> [n_levels]

    n_levels defaults to 3  (100x50, 200x100, 400x200)
"""

import subprocess
import sys
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

# ── configuration ─────────────────────────────────────────────────────────────
BINARY    = Path(sys.argv[1]) if len(sys.argv) > 1 else Path("./test_boussinesq_onesided_rb_box")
N_LEVELS  = int(sys.argv[2]) if len(sys.argv) > 2 else 3

PARAMS    = (Path(__file__).parent / "params.input").resolve()
OUTDIR    = Path(__file__).parent / "convergence_runs"
PLOT_FILE = OUTDIR / "convergence.pdf"

BASE_NX, BASE_NY = 100, 50
resolutions = [(BASE_NX * 2**k, BASE_NY * 2**k) for k in range(N_LEVELS)]


# ── helpers ───────────────────────────────────────────────────────────────────
def run(nx, ny):
    name   = f"boussinesq_{nx}x{ny}"
    rundir = OUTDIR / name
    rundir.mkdir(parents=True, exist_ok=True)
    csv    = rundir / f"{name}_sherwood.csv"

    if csv.exists():
        print(f"  skip  {nx:4d}×{ny:3d}  (cached)")
        return csv

    print(f"  run   {nx:4d}×{ny:3d} …", flush=True)
    result = subprocess.run(
        [
            str(BINARY.resolve()),
            str(PARAMS),
            "-Grid.Cells",               f"{nx} {ny}",
            "-Problem.Name",             name,
            "-Output.VtkOutputInterval", "1e30",   # suppress mid-run VTK
        ],
        cwd=str(rundir),
        capture_output=True,
        text=True,
    )
    if result.returncode != 0:
        print(f"\nSimulation failed for {nx}×{ny}:")
        print(result.stderr[-3000:])
        sys.exit(1)
    return csv


def load(csv):
    data = np.genfromtxt(csv, delimiter=",", skip_header=1)
    return data[:, 0], data[:, 2]   # time, F = Sh/Ra


# ── run all levels ─────────────────────────────────────────────────────────────
OUTDIR.mkdir(parents=True, exist_ok=True)
print(f"Binary : {BINARY.resolve()}")
print(f"Params : {PARAMS}")
print(f"Levels : {N_LEVELS}  →  {resolutions}")
print()

datasets = []
for nx, ny in resolutions:
    csv     = run(nx, ny)
    t, Sh   = load(csv)
    datasets.append({"nx": nx, "ny": ny, "h": 1.0 / ny, "t": t, "Sh": Sh})  # Sh here is F = Sh/Ra

# ── convergence metric: L² error vs. finest reference ─────────────────────────
t_ref, Sh_ref = datasets[-1]["t"], datasets[-1]["Sh"]

for d in datasets[:-1]:
    Sh_i    = np.interp(t_ref, d["t"], d["Sh"])
    d["err"] = np.sqrt(np.trapz((Sh_i - Sh_ref) ** 2, t_ref))

# estimated convergence order from last two coarse levels
if len(datasets) >= 3:
    d0, d1 = datasets[-3], datasets[-2]
    p_est = np.log(d0["err"] / d1["err"]) / np.log(d0["h"] / d1["h"])
    print(f"\nEstimated convergence order (last two coarse levels): p ≈ {p_est:.2f}")


# ── matplotlib ────────────────────────────────────────────────────────────────
plt.rcParams.update(
    {
        "text.usetex":      True,
        "font.family":      "serif",
        "font.size":        11,
        "axes.labelsize":   13,
        "axes.titlesize":   12,
        "legend.fontsize":  10,
        "xtick.labelsize":  10,
        "ytick.labelsize":  10,
        "figure.dpi":       150,
    }
)

colors = plt.cm.plasma(np.linspace(0.15, 0.85, len(datasets)))

fig, (ax_sh, ax_cv) = plt.subplots(1, 2, figsize=(13, 5))
fig.suptitle(
    r"Boussinesq one-sided RB -- spatial convergence"
    rf" ($\mathrm{{Ra}} = 500$,\ $T = 10$)",
    fontsize=13,
    y=1.02,
)

# ── left panel: Sh(t) ─────────────────────────────────────────────────────────
for d, c in zip(datasets, colors):
    ax_sh.plot(
        d["t"], d["Sh"],
        color=c, lw=1.6,
        label=rf"${d['nx']}\times{d['ny']}$",
    )

ax_sh.set_xscale("log")
ax_sh.set_yscale("log")
ax_sh.set_xlabel(r"$t$")
ax_sh.set_ylabel(r"$F(t) = \mathrm{Sh}(t)\,/\,\mathrm{Ra}$")
ax_sh.set_title(r"Dimensionless dissolution flux")
ax_sh.legend(title=r"Grid $N_x \times N_y$", framealpha=0.9)
ax_sh.grid(True, which="both", ls=":", alpha=0.45)
ax_sh.xaxis.set_minor_locator(ticker.LogLocator(subs="auto"))
ax_sh.yaxis.set_minor_locator(ticker.LogLocator(subs="auto"))

# ── right panel: convergence ───────────────────────────────────────────────────
h_vals  = [d["h"]   for d in datasets[:-1]]
err_vals = [d["err"] for d in datasets[:-1]]

ax_cv.loglog(
    h_vals, err_vals,
    "o-", color="steelblue", ms=7, lw=2,
    label=r"$\|F_h - F_{\mathrm{ref}}\|_{L^2(0,T)}$",
    zorder=3,
)

# reference slopes anchored at the coarsest data point
h_arr = np.array(h_vals)
slope_styles = [
    (1, "--", "dimgray",  r"$\mathcal{O}(h^1)$"),
    (2, ":",  "firebrick", r"$\mathcal{O}(h^2)$"),
]
for order, ls, col, lbl in slope_styles:
    ref = err_vals[0] * (h_arr / h_arr[0]) ** order
    ax_cv.loglog(h_arr, ref, ls=ls, color=col, lw=1.4, label=lbl)

ax_cv.set_xlabel(r"$h = H \,/\, N_y$")
ax_cv.set_ylabel(r"$L^2$ error in $F = \mathrm{Sh}\,/\,\mathrm{Ra}$")
ax_cv.set_title(r"Spatial convergence")
ax_cv.legend(framealpha=0.9)
ax_cv.grid(True, which="both", ls=":", alpha=0.45)
ax_cv.xaxis.set_minor_locator(ticker.LogLocator(subs="auto"))
ax_cv.yaxis.set_minor_locator(ticker.LogLocator(subs="auto"))
ax_cv.invert_xaxis()   # coarse → fine left → right

fig.tight_layout()
fig.savefig(PLOT_FILE, bbox_inches="tight")
print(f"\nPlot saved → {PLOT_FILE}")
plt.show()