#!/usr/bin/env python3
# SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

"""
Onset of one-sided convection: growth of a perturbation of the diffusive boundary layer, DuMux compared with the
linearized problem of linear_stability.py.

Usage: python3 run_benchmark.py ./test_1p2c_onesidedconvection_box [output directory of the figures]

Prints the errors of the amplification and writes the figures onesidedconvection_profile.png (the perturbation profile at
several times) and onesidedconvection_growthrate.png (the growth rate of the maximum, with the markers of Riaz et al. 2006).
"""

import csv
import os
import shutil
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from linear_stability import LinearizedProblem, dominant_mode

RAYLEIGH, INITIAL_TIME, AMPLITUDE = 500.0, 0.01, 1e-6  # AMPLITUDE: Problem.Amplitude
GRIDS = {8: ("8 32", 1.03**2), 16: ("16 64", 1.03), 32: ("32 128", 1.03**0.5)}  # cells: (cells in z, grading towards the top)
PROFILE_TIMES = [0.1, 0.25, 0.5, 1.0, 2.0]
PROFILE_COLORS = ["#86b6ef", "#5598e7", "#2a78d6", "#1c5cab", "#0d366b"]  # one hue from light to dark
DEPTH_SHOWN = 0.25
GROWTH_RATE_WIDTHS = {4: "0 0.25", 8: "0 0.125"}  # k / pi: box width, k = pi H / L; until t = 4
GROWTH_RATE_COLORS = {4: "#2a78d6", 8: "#eb6834"}
RIAZ_DATA = "riaz2006_fig4_ivp.csv"  # markers of the initial-value solution in Figure 4 of Riaz et al. (2006)
INK, INK_SECONDARY, SURFACE = "#0b0b0b", "#52514e", "#ffffff"


def run(executable, name, n, dt, *args):
    """Run the benchmark on the grid of n x 5 n cells and return the time series and the final amplitude profile"""
    cells_z, grading = GRIDS[n]
    subprocess.run(
        [executable, "params.input", "-Problem.Name", name, "-Grid.Cells0", str(n), "-Grid.Cells1", cells_z,
         "-Grid.Grading1", f"1 {-grading:.15g}", "-TimeLoop.DimensionlessTimeStepSize", str(dt), "-Vtk.Precision", "Float64", *args],
        capture_output=True, text=True, check=True, env=dict(os.environ, OMP_NUM_THREADS="1"),
    )
    return (np.genfromtxt(f"{name}.csv", delimiter=",", names=True),
            np.genfromtxt(f"{name}_profile.csv", delimiter=",", names=True))


def new_axes(height):
    if shutil.which("latex") and shutil.which("dvipng"):
        plt.rcParams.update({"text.usetex": True, "font.family": "serif", "font.serif": ["Computer Modern Roman"], "font.size": 12})
    else:
        plt.rcParams.update({"font.family": "serif", "mathtext.fontset": "cm", "font.size": 12})
    plt.rcParams.update({"text.color": INK, "axes.labelcolor": INK, "xtick.color": INK_SECONDARY, "ytick.color": INK_SECONDARY,
                         "axes.edgecolor": INK_SECONDARY})
    fig, ax = plt.subplots(figsize=(8.0, height), facecolor=SURFACE)
    ax.set_facecolor(SURFACE)
    ax.grid(color="#e4e3df", lw=0.6)
    ax.set_axisbelow(True)
    ax.spines[["top", "right"]].set_visible(False)
    return fig, ax


def legend(ax, entries):
    """entries: (label, marker, open marker?) in the neutral color; no marker means a line"""
    handles = [plt.Line2D([], [], color=INK_SECONDARY, lw=0 if marker else 2.0, marker=marker, markersize=6,
                          markerfacecolor=SURFACE if hollow else INK_SECONDARY, markeredgewidth=1.3) for _, marker, hollow in entries]
    ax.legend(handles, [label for label, _, _ in entries], frameon=False, fontsize=10, loc="lower right")


def plot_profiles(profiles, reference, filename):
    """Shape of the perturbation over the depth at several times: linearized problem and DuMux"""
    fig, ax = new_axes(5.0)
    depth = np.linspace(0.0, DEPTH_SHOWN, 1201)
    fine = np.linspace(0.0, 1.0, 20001)
    normalize = lambda t: np.abs(reference(t)(fine)).max()  # every profile is divided by its maximum

    initial = reference(INITIAL_TIME)(depth) / normalize(INITIAL_TIME)  # thin and dotted: the initial perturbation
    ax.plot(initial, depth, color=INK_SECONDARY, lw=1.0, ls=":")
    ax.annotate(r"$t_0$", (1.0, depth[np.argmax(initial)]), xytext=(7, 0), textcoords="offset points", va="center",
                fontsize=10, color=INK_SECONDARY)
    for t, color in zip(PROFILE_TIMES, PROFILE_COLORS):
        curve = reference(t)(depth) / normalize(t)
        ax.plot(curve, depth, color=color, lw=2.0, solid_capstyle="round")
        p = profiles[t]
        shown = np.flatnonzero(p["depth"] <= DEPTH_SHOWN)
        marked = shown[np.unique(np.linspace(0, len(shown) - 1, 26).astype(int))]
        ax.plot(p["amplitude"][marked] / np.abs(p["amplitude"]).max(), p["depth"][marked], linestyle="none", marker="o",
                markersize=4.5, markerfacecolor=SURFACE, markeredgecolor=color, markeredgewidth=1.3)
        ax.annotate(rf"$t = {t:g}$", (1.0, depth[np.argmax(curve)]), xytext=(7, 0), textcoords="offset points", va="center",
                    fontsize=10, color=INK_SECONDARY)

    ax.set_ylim(DEPTH_SHOWN, -0.012)  # room above the top boundary at eta = 0, so that nothing is clipped there
    ax.set_xlim(-0.05, 1.25)
    ax.set_xlabel(r"normalized amplitude of the perturbation $\hat c / \max|\hat c|$")
    ax.set_ylabel(r"depth $\eta = d / H$")
    legend(ax, [("linear stability solver", None, False), (r"DuMux, $16 \times 80$ cells", "o", True)])
    fig.tight_layout()
    fig.savefig(filename, dpi=120, facecolor=SURFACE)


def plot_growth_rates(results, filename):
    """Growth rate d ln max|c| / dt of the maximum of the perturbation: linearized problem, Riaz et al. (2006) and DuMux"""
    with open(RIAZ_DATA, newline="", encoding="utf-8") as f:
        riaz = list(csv.DictReader(line for line in f if not line.startswith("#")))
    fig, ax = new_axes(4.6)
    times = np.linspace(INITIAL_TIME, 4.0, 3000)
    print("\nGrowth rate d ln max|c| / dt at t = 1, 2, 4 (Riaz et al. interpolated between their markers):")
    for k, color in GROWTH_RATE_COLORS.items():
        problem = LinearizedProblem(RAYLEIGH, k * np.pi, INITIAL_TIME, times[-1])
        solver = np.gradient(np.log(problem.max_norm(problem.solve(dominant_mode(RAYLEIGH, INITIAL_TIME), times))), times)
        series = results[f"k {k} pi"][0]
        dumux = np.gradient(np.log(series["maxAmplitude"]), series["time"])
        t_riaz = np.array([float(row["t"]) for row in riaz if int(row["k_over_pi"]) == k])
        rate_riaz = np.array([float(row["growth_rate"]) for row in riaz if int(row["k_over_pi"]) == k])

        ax.plot(times, solver, color=color, lw=2.0, solid_capstyle="round")
        ax.plot(t_riaz, rate_riaz, linestyle="none", marker="v", markersize=6, markerfacecolor=SURFACE, markeredgecolor=color,
                markeredgewidth=1.3)
        ax.plot(series["time"][::40], dumux[::40], linestyle="none", marker="o", markersize=4.0, color=color)
        ax.annotate(rf"$k = {k}\pi$", (times[-1], solver[-1]), xytext=(6, 6), textcoords="offset points", fontsize=11,
                    color=INK_SECONDARY)

        late = series["time"] > 0.8
        error = dumux - np.interp(series["time"], times, solver)
        error_riaz = (rate_riaz - np.interp(t_riaz, times, solver))[t_riaz > 0.8]
        print(f"k = {k} pi: linearized problem " + ", ".join(f"{np.interp(t, times, solver):.3f}" for t in (1, 2, 4))
              + "; DuMux " + ", ".join(f"{np.interp(t, series['time'], dumux):.3f}" for t in (1, 2, 4))
              + "; Riaz " + ", ".join(f"{np.interp(t, t_riaz, rate_riaz):.3f}" for t in (1, 2, 3.5)))
        print(f"  for t > 0.8, difference to the linearized problem: DuMux largest {np.abs(error[late]).max():.3f}, rms "
              f"{np.sqrt(np.mean(error[late]**2)):.3f}; Riaz largest {np.abs(error_riaz).max():.3f}, rms "
              f"{np.sqrt(np.mean(error_riaz**2)):.3f}")

    ax.set_ylim(-0.5, 2.5)
    ax.set_xlim(0.0, times[-1] + 0.05)
    ax.set_xlabel(r"$t$ $[\phi H / U]$")
    ax.set_ylabel(r"growth rate $\sigma = \mathrm{d} \ln \max|\hat c| / \mathrm{d}t$")
    legend(ax, [("linear stability solver", None, False), ("Riaz et al. (2006), Fig. 4, IVP", "v", True),
                (r"DuMux, $16 \times 80$ cells", "o", False)])
    fig.tight_layout()
    fig.savefig(filename, dpi=120, facecolor=SURFACE)


def main():
    if len(sys.argv) < 2:
        sys.exit(__doc__)
    executable = os.path.abspath(sys.argv[1])
    image_dir = sys.argv[2] if len(sys.argv) > 2 else "."

    # the grid refinement with the time step proportional to the square of the cell size, so that the error of the first-order
    # time integration decreases like the second-order error in space, and the profiles at the times of the figure and the
    # growth rates for two wavenumbers, on 16 x 80 cells with the time step of the test
    jobs = {f"grid {n}": (f"osc_grid_{n}", n, 4e-3*(min(GRIDS)/n)**2, ()) for n in GRIDS}
    jobs.update({f"t {t:g}": (f"osc_t_{t:g}", 16, 2e-3, ("-TimeLoop.DimensionlessEndTime", str(t))) for t in PROFILE_TIMES})
    jobs.update({f"k {k} pi": (f"osc_k{k}pi", 16, 2e-3, ("-Grid.Positions0", width, "-TimeLoop.DimensionlessEndTime", "4"))
                 for k, width in GROWTH_RATE_WIDTHS.items()})
    with ThreadPoolExecutor(max_workers=len(jobs)) as pool:
        futures = {key: pool.submit(run, executable, name, n, dt, *args) for key, (name, n, dt, args) in jobs.items()}
        results = {key: future.result() for key, future in futures.items()}

    # the amplification A(t) / A(t0) of the L2 norm of the perturbation
    times = results[f"grid {max(GRIDS)}"][0]["time"]
    problem = LinearizedProblem(RAYLEIGH, 8 * np.pi, INITIAL_TIME, times[-1])
    c = problem.solve(dominant_mode(RAYLEIGH, INITIAL_TIME), times)
    log_reference = np.log(problem.l2_norm(c) / problem.l2_norm(c[0]))

    def error(key, t):
        series = results[key][0]
        return np.interp(t, series["time"], np.log(series["amplitude"] / series["amplitude"][0])
                         - np.interp(series["time"], times, log_reference))

    print("Reference ln A(t)/A(t0) at t = 0.5, 1, 2: " + ", ".join(f"{np.interp(t, times, log_reference):.5f}" for t in (0.5, 1, 2)))
    print(f"Minimum of A(t)/A(t0): {np.exp(log_reference.min()):.3f} at t = {times[np.argmin(log_reference)]:.3f}")
    print("\nError of ln A under refinement with the time step proportional to the square of the cell size")
    print("| Cells | time step | t = 0.5 | t = 1 | t = 2 | rate |")
    print("|------|------|------|------|------|------|")
    previous = None
    for n in GRIDS:
        current = error(f"grid {n}", 2)
        rate = f"{np.log2(abs(previous / current)):.2f}" if previous else "-"
        previous = current
        print(f"| ${n} \\times {5 * n}$ | {4e-3*(min(GRIDS)/n)**2:g} | " + " | ".join(f"{error(f'grid {n}', t):.4f}" for t in (0.5, 1, 2)) + f" | {rate} |")
    print(f"Error at the setting of the test (16 x 80, time step 2e-3): {error('t 2', 2):.4f}")

    # the linearized perturbation at the times of the figure, for the same initial perturbation (amplitude 1)
    all_times = np.array([INITIAL_TIME] + PROFILE_TIMES)
    coefficients = dict(zip(all_times, problem.solve(dominant_mode(RAYLEIGH, INITIAL_TIME), all_times)))
    reference = lambda t: (lambda depth: problem.profile(coefficients[t], depth))

    print("\nAmplitude profile, 16 x 80 cells, time step 2e-3: maximum of the amplitude c / eps (linearized problem, DuMux) and\n"
          "largest difference of the normalized profiles:")
    profiles = {t: results[f"t {t:g}"][1] for t in PROFILE_TIMES}
    for t, p in profiles.items():
        peak = np.abs(reference(t)(np.linspace(0.0, 1.0, 20001))).max()
        difference = np.abs(p["amplitude"] / np.abs(p["amplitude"]).max() - reference(t)(p["depth"]) / peak).max()
        print(f"t = {t:g}: {peak:.4f}, {np.abs(p['amplitude']).max() / AMPLITUDE:.4f}, {difference:.1e}")

    plot_profiles(profiles, reference, f"{image_dir}/onesidedconvection_profile.png")
    plot_growth_rates(results, f"{image_dir}/onesidedconvection_growthrate.png")


if __name__ == "__main__":
    main()
