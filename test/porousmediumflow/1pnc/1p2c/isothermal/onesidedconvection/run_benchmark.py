#!/usr/bin/env python3
# SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

"""
One-sided convection benchmark: amplification of a finger mode under grid and time-step refinement,
compared with the linearized problem solved by linear_stability.py, and the diffusive boundary layer
without perturbation.

Usage: python3 run_benchmark.py ./test_1p2c_onesidedconvection_box [output directory of the figures]
"""

import os
import subprocess
import sys
import xml.etree.ElementTree as ET
from concurrent.futures import ThreadPoolExecutor

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import BarycentricInterpolator

from linear_stability import LinearizedProblem, dominant_mode

RAYLEIGH, WAVENUMBER, INITIAL_TIME = 500.0, 8 * np.pi, 0.01
PROBES = (0.5, 1.0, 2.0)

# cells per half wavelength, cells in the lower and the upper zone, grading towards the top
GRIDS = {"8 x 40": (8, "8 32", 1.03**2), "16 x 80": (16, "16 64", 1.03), "32 x 160": (32, "32 128", 1.03**0.5)}
TIME_STEPS = [2e-3, 1e-3, 5e-4]


def run(executable, name, cells, dt, *args):
    """Run the benchmark and return the time series and the final amplitude profile"""
    nx, nz, grading = cells
    env = dict(os.environ, VECLIB_MAXIMUM_THREADS="1", OMP_NUM_THREADS="1")
    subprocess.run(
        [executable, "params.input", "-Problem.Name", name, "-Grid.Cells0", str(nx), "-Grid.Cells1", nz,
         "-Grid.Grading1", f"1 {-grading:.15g}", "-TimeLoop.DimensionlessTimeStepSize", str(dt),
         "-Vtk.Precision", "Float64", *args],
        capture_output=True, text=True, check=True, env=env,
    )
    series = np.genfromtxt(f"{name}.csv", delimiter=",", names=True)
    profile = np.genfromtxt(f"{name}_profile.csv", delimiter=",", names=True)
    return series, profile


def read_vtu(filename):
    """Return the points and the point data of an ASCII VTU file"""
    piece = ET.parse(filename).getroot().find("UnstructuredGrid/Piece")
    to_array = lambda array: np.array(array.text.split(), dtype=float).reshape(
        -1, int(array.get("NumberOfComponents", "1")))
    points = to_array(piece.find("Points/DataArray"))
    point_data = {array.get("Name"): to_array(array) for array in piece.find("PointData")}
    return points, point_data


def tensor_grid(points, values):
    """Arrange vertex values of a tensor-product grid as an array (z, x)"""
    x, i = np.unique(np.round(points[:, 0], 12), return_inverse=True)
    z, j = np.unique(np.round(points[:, 1], 12), return_inverse=True)
    field = np.empty((len(z), len(x)))
    field[j, i] = values
    return x, z, field


def plot_mode(vtu, profile, reference, filename):
    """Perturbation of the horizontal mean with streamlines, and the amplitude profile against the reference"""
    points, data = read_vtu(vtu)
    x, z, massfraction = tensor_grid(points, data["X^Solute_liq"][:, 0])
    perturbation = massfraction - massfraction.mean(axis=1, keepdims=True)
    perturbation /= np.abs(perturbation).max()
    _, _, u = tensor_grid(points, data["velocity_liq (m/s)"][:, 0])
    _, _, w = tensor_grid(points, data["velocity_liq (m/s)"][:, 1])
    shown = z >= 0.6

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(7, 4.2), layout="constrained",
                                   gridspec_kw={"width_ratios": [1, 1.4]})
    image = ax1.contourf(x, 1.0 - z[shown], perturbation[shown], levels=np.linspace(-1, 1, 21), cmap="RdBu_r")
    xu = np.linspace(x[0], x[-1], 41)
    depth = np.linspace(0.0, 0.4, 129)
    resample = lambda f: np.array([np.interp(xu, x, row) for row in
                                   np.array([np.interp(1.0 - depth, z, col) for col in f.T]).T])
    ax1.streamplot(xu, depth, resample(u), -resample(w), color="k", linewidth=0.5, density=0.7, arrowsize=0.6)
    ax1.set_xlim(x[0], x[-1])
    ax1.set_ylim(0.4, 0.0)
    ax1.set_aspect("equal")
    ax1.set_xlabel(r"$x / H$")
    ax1.set_ylabel(r"depth $d / H$")
    fig.colorbar(image, ax=ax1, label=r"$(X - \bar X) / \max|X - \bar X|$", ticks=[-1, 0, 1], shrink=0.8)

    fine = np.linspace(0.0, 0.4, 801)
    ax2.plot(reference(fine), fine, "k-", lw=2, label="linearized problem")
    ax2.plot(profile["amplitude"] / np.abs(profile["amplitude"]).max(), profile["depth"], "C1--",
             label="DuMux, 32 x 160")
    ax2.set_ylim(0.4, 0.0)
    ax2.set_xlabel(r"amplitude $\hat c / \max|\hat c|$ at $t = 2$")
    ax2.set_ylabel(r"depth $d / H$")
    ax2.legend(fontsize=8, loc="lower right")
    fig.savefig(filename, dpi=150)


def main():
    if len(sys.argv) < 2:
        sys.exit(__doc__)
    executable = os.path.abspath(sys.argv[1])
    image_dir = sys.argv[2] if len(sys.argv) > 2 else "."

    jobs = {f"grid {g}": (f"osc_grid_{cells[0]}", cells, 1e-3, ()) for g, cells in GRIDS.items()}
    jobs.update({f"dt {dt:g}": (f"osc_dt_{dt:g}", GRIDS["16 x 80"], dt, ()) for dt in TIME_STEPS if dt != 1e-3})
    jobs["diffusion"] = ("osc_diffusion", GRIDS["16 x 80"], 1e-3, ("-Problem.Amplitude", "0"))
    with ThreadPoolExecutor(max_workers=len(jobs)) as pool:
        futures = {key: pool.submit(run, executable, name, cells, dt, *args)
                   for key, (name, cells, dt, args) in jobs.items()}
        results = {key: future.result() for key, future in futures.items()}
    results["dt 0.001"] = results["grid 16 x 80"]

    times = results["dt 0.0005"][0]["time"]
    problem = LinearizedProblem(RAYLEIGH, WAVENUMBER)
    c = problem.solve(dominant_mode(RAYLEIGH, INITIAL_TIME), INITIAL_TIME, times)
    log_reference = np.log(problem.l2_norm(c) / problem.l2_norm(c[0]))

    def error(key):
        series = results[key][0]
        return np.log(series["amplitude"] / series["amplitude"][0]) - np.interp(series["time"], times, log_reference)

    def errors_at(key):
        return [np.interp(t, results[key][0]["time"], error(key)) for t in PROBES]

    print("Reference ln A(t)/A(t0) at t = 0.5, 1, 2: " + ", ".join(f"{np.interp(t, times, log_reference):.5f}" for t in PROBES))
    print(f"Minimum of A(t)/A(t0): {np.exp(log_reference.min()):.3f} at t = {times[np.argmin(log_reference)]:.3f}")

    time_error = [2.0 * (a - b) for a, b in zip(errors_at("dt 0.001"), errors_at("dt 0.0005"))]
    print("\nGrid refinement, time step 1e-3: error of ln A")
    print("| Cells | t = 0.5 | t = 1 | t = 2 | t = 2, without time error | rate |")
    print("|------|------|------|------|------|------|")
    previous = None
    for name in GRIDS:
        e = errors_at(f"grid {name}")
        spatial = e[2] - time_error[2]
        rate = f"{np.log2(abs(previous / spatial)):.2f}" if previous else "-"
        previous = spatial
        label = name.replace(" x ", r" \times ")
        print(f"| ${label}$ | {e[0]:.4f} | {e[1]:.4f} | {e[2]:.4f} | {spatial:.4f} | {rate} |")

    print("\nTime-step refinement, 16 x 80 cells: error of ln A")
    print("| Time step | t = 0.5 | t = 1 | t = 2 |")
    print("|------|------|------|------|")
    for dt in TIME_STEPS:
        e = errors_at(f"dt {dt:g}")
        print(f"| ${dt:g}$ | {e[0]:.4f} | {e[1]:.4f} | {e[2]:.4f} |")
    print("Implicit Euler error at time step 1e-3, assuming first order: " + ", ".join(f"{v:.4f}" for v in time_error))

    print("\nLargest error of ln A before t = 0.1: " + ", ".join(
        f"{name}: {error(f'grid {name}')[results[f'grid {name}'][0]['time'] < 0.1].max():.4f}" for name in GRIDS))

    end = problem.solve(dominant_mode(RAYLEIGH, INITIAL_TIME), INITIAL_TIME, np.array([INITIAL_TIME, times[-1]]))[-1]
    interpolant = BarycentricInterpolator(problem.z, end)
    scale = np.abs(interpolant(np.linspace(0.0, 1.0, 20001))).max()
    reference_profile = lambda depth: interpolant(depth) / scale
    print("Largest difference of the normalized amplitude profile at t = 2: " + ", ".join(
        f"{name}: {np.abs(p['amplitude'] / np.abs(p['amplitude']).max() - reference_profile(p['depth'])).max():.1e}"
        for name, p in ((n, results[f"grid {n}"][1]) for n in GRIDS)))

    series = results["diffusion"][0]
    relative = series["flux"] / series["exactFlux"] - 1.0
    late = series["time"] > 0.1
    print(f"\nWithout perturbation, 16 x 80 cells, time step 1e-3: largest horizontal variation in a row "
          f"{series['maxRowDeviation'].max():.1e}; relative error of the flux after the first step {relative[1]:.3f}, "
          f"largest for t > 0.1 {np.abs(relative[late]).max():.1e}, at t = 2 {abs(relative[-1]):.1e}")

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 3.8))
    ax1.semilogy(times, np.exp(log_reference), "k-", lw=2, label="linearized problem")
    ax2.axhline(0.0, color="0.6", lw=0.8)
    for name, color in zip(GRIDS, ("C0", "C1", "C2")):
        s = results[f"grid {name}"][0]
        ax1.semilogy(s["time"], s["amplitude"] / s["amplitude"][0], "--", color=color, label=f"DuMux, {name}")
        ax2.plot(s["time"], error(f"grid {name}"), color=color, label=name)
    ax1.set_xlabel(r"$t$ $[\phi H / U]$")
    ax1.set_ylabel(r"amplification $A(t) / A(t_0)$")
    ax1.legend(fontsize=8)
    ax2.set_xlabel(r"$t$ $[\phi H / U]$")
    ax2.set_ylabel(r"error of $\ln A$")
    ax2.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(f"{image_dir}/onesidedconvection_growth.png", dpi=150)

    plot_mode("osc_grid_32-00001.vtu", results["grid 32 x 160"][1], reference_profile,
              f"{image_dir}/onesidedconvection_mode.png")


if __name__ == "__main__":
    main()
