#!/usr/bin/env python3
# SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later
"""Run Mandel's problem and plot the solution against the analytical solution.

Usage:
  python3 plot_results.py <executable> [output directory]

Runs the executable with params.input of the current directory and writes mandel_pressure.png
(pressure profiles and the pressure close to the center over time) and mandel_displacement.png
(horizontal displacement of the plate) to the output directory (default: current directory).
Requires numpy, matplotlib and pyvista.
"""

import configparser
import subprocess
import sys
from pathlib import Path

PROFILE_TIMES = [10.0, 100.0, 1000.0, 5000.0, 20000.0]
NAME = "mandel_plot"


def readParameters(fileName):
    parser = configparser.ConfigParser(inline_comment_prefixes=("#",))
    parser.optionxform = str
    parser.read(fileName)
    value = lambda group, key: float(parser[group][key].split()[0])
    return {
        "F": value("Problem", "Force"),
        "E": value("Problem", "EModulus"),
        "nu": value("Problem", "Nu"),
        "alpha": value("Problem", "BiotCoefficient"),
        "phi0": value("Problem", "InitialPorosity"),
        "Kf": value("Problem", "FluidBulkModulus"),
        "k": value("Problem", "Permeability"),
        "mu": value("Component", "LiquidDynamicViscosity"),
        "a": float(parser["Grid"]["UpperRight"].split()[0]),
    }


class AnalyticalSolution:
    """Series solution of Mandel's problem (Abousleiman et al., 1996)."""

    def __init__(self, params, numTerms=240):
        import numpy as np

        self.F, self.a = params["F"], params["a"]
        E, nu, alpha, phi0 = params["E"], params["nu"], params["alpha"], params["phi0"]
        self.G = E/(2.0*(1.0 + nu))
        K = E/(3.0*(1.0 - 2.0*nu))
        M = 1.0/(phi0/params["Kf"] + (alpha - phi0)*(1.0 - alpha)/K)
        Ku = K + alpha**2*M
        B = alpha*M/Ku
        self.nu, self.nuu = nu, (3.0*Ku - 2.0*self.G)/(2.0*(3.0*Ku + self.G))
        self.c = params["k"]/params["mu"]*M*(K + 4.0/3.0*self.G)/(Ku + 4.0/3.0*self.G)
        self.pPlus = self.F*B*(1.0 + self.nuu)/(3.0*self.a)

        ratio = (1.0 - nu)/(self.nuu - nu)
        f = lambda x: np.tan(x)/x - ratio
        delta = 1e-8
        lower = np.array([delta] + [(n - 0.5)*np.pi + delta for n in range(1, numTerms)])
        upper = np.array([(n + 0.5)*np.pi - delta for n in range(numTerms)])
        for _ in range(100):
            mid = 0.5*(lower + upper)
            positive = f(mid) > 0.0
            upper = np.where(positive, mid, upper)
            lower = np.where(positive, lower, mid)
        self.beta = 0.5*(lower + upper)

    def _decay(self, t):
        import numpy as np
        return np.exp(-self.beta**2*self.c*t/self.a**2)

    def pressure(self, x, t):
        import numpy as np
        b = self.beta
        coefficients = np.sin(b)/(b - np.sin(b)*np.cos(b))*self._decay(t)
        return 2.0*self.pPlus*np.sum(coefficients[:, None]*(np.cos(np.outer(b, x)/self.a) - np.cos(b)[:, None]), axis=0)

    def displacementX(self, x, t):
        import numpy as np
        b = self.beta
        denominator = b - np.sin(b)*np.cos(b)
        uniformStrainSum = np.sum(np.sin(b)*np.cos(b)/denominator*self._decay(t))
        xSum = np.sum((np.cos(b)/denominator*self._decay(t))[:, None]*np.sin(np.outer(b, x)/self.a), axis=0)
        return (self.nu/2.0 - self.nuu*uniformStrainSum)*self.F/(self.G*self.a)*x + self.F/self.G*xSum


def run(executable, checkPoints):
    subprocess.run(
        [
            str(executable), "params.input",
            "-TimeLoop.TEnd", str(checkPoints[-1]),
            "-TimeLoop.MaxTimeStepSize", "100",
            "-TimeLoop.CheckPoints", " ".join(f"{t:g}" for t in checkPoints),
            "-Vtk.OutputName", NAME,
        ],
        check=True, stdout=subprocess.DEVNULL,
    )


def readOutput(pvd):
    """Return a dictionary mapping output time to mesh."""
    import pyvista as pv
    reader = pv.get_reader(pvd)
    meshes = {}
    for t in reader.time_values:
        reader.set_active_time_value(t)
        meshes[t] = reader.read()[0]
    return meshes


def main():
    if len(sys.argv) < 2:
        sys.stderr.write(__doc__)
        sys.exit(1)

    import numpy as np
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    executable = Path(sys.argv[1]).resolve()
    outdir = Path(sys.argv[2]).resolve() if len(sys.argv) > 2 else Path.cwd()

    historyTimes = list(np.geomspace(1.0, PROFILE_TIMES[-1], 40))
    checkPoints = sorted(set(historyTimes) | set(PROFILE_TIMES))
    run(executable, checkPoints)

    exact = AnalyticalSolution(readParameters("params.input"))
    flow = readOutput(f"{NAME}_flow.pvd")
    mech = readOutput(f"{NAME}_mech.pvd")
    outputTime = lambda t: min(flow, key=lambda s: abs(s - t))

    colors = plt.rcParams["axes.prop_cycle"].by_key()["color"]
    xFine = np.linspace(0.0, exact.a, 400)
    legendTitle = "lines: analytical, markers: DuMux"

    # pressure in the cells closest to y = b/2 (the pressure is independent of y)
    mesh = flow[outputTime(0.0)]
    centers = mesh.cell_centers().points
    rowY = centers[np.argmin(np.abs(centers[:, 1] - 0.5*mesh.bounds[3])), 1]
    row = np.where(np.isclose(centers[:, 1], rowY))[0]
    row = row[np.argsort(centers[row, 0])]
    xCells = centers[row, 0]

    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))
    for i, t in enumerate(PROFILE_TIMES):
        p = np.asarray(flow[outputTime(t)].cell_data["p"]).ravel()[row]
        axes[0].plot(xFine/exact.a, exact.pressure(xFine, t)/exact.pPlus, color=colors[i], label=f"t = {t:g} s")
        axes[0].plot(xCells/exact.a, p/exact.pPlus, "o", color=colors[i], markerfacecolor="none", markersize=4)
    axes[0].set_xlabel("$x/a$")
    axes[0].set_xlim(0, 1)
    axes[0].set_ylim(bottom=0)

    times = np.array(sorted(t for t in flow if t > 0.0))
    pNumerical = [np.asarray(flow[t].cell_data["p"]).ravel()[row[0]] for t in times]
    tFine = np.geomspace(times[0], times[-1], 400)
    pExact = [exact.pressure(np.array([xCells[0]]), t)[0] for t in tFine]
    axes[1].semilogx(tFine, np.array(pExact)/exact.pPlus, color=colors[0], label=f"$x = {xCells[0]:g}$ m")
    axes[1].semilogx(times, np.array(pNumerical)/exact.pPlus, "o", color=colors[0], markerfacecolor="none", markersize=4)
    axes[1].set_xlabel("$t$ [s]")
    for ax in axes:
        ax.set_ylabel("$p/p^+$")
        ax.grid(alpha=0.3)
        ax.legend(title=legendTitle, fontsize=9, title_fontsize=9)
    fig.tight_layout()
    fig.savefig(outdir / "mandel_pressure.png", dpi=150)

    # horizontal displacement of the vertices on the plate
    mesh = mech[outputTime(0.0)]
    points = mesh.points
    top = np.where(np.isclose(points[:, 1], mesh.bounds[3]))[0]
    top = top[np.argsort(points[top, 0])]

    fig, ax = plt.subplots(figsize=(7, 4.5))
    for i, t in enumerate(PROFILE_TIMES):
        ux = np.asarray(mech[outputTime(t)].point_data["u"])[top, 0]
        ax.plot(xFine/exact.a, exact.displacementX(xFine, t), color=colors[i], label=f"t = {t:g} s")
        ax.plot(points[top, 0]/exact.a, ux, "o", color=colors[i], markerfacecolor="none", markersize=4)
    ax.set_xlabel("$x/a$")
    ax.set_ylabel("$u_x$ [m]")
    ax.set_xlim(0, 1)
    ax.grid(alpha=0.3)
    ax.legend(title=legendTitle, fontsize=9, title_fontsize=9)
    fig.tight_layout()
    fig.savefig(outdir / "mandel_displacement.png", dpi=150)

    print(f"Wrote {outdir / 'mandel_pressure.png'} and {outdir / 'mandel_displacement.png'}")


if __name__ == "__main__":
    main()
