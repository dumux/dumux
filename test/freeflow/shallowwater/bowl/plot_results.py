#!/usr/bin/env python3
# SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later
"""Oscillation in a paraboloid bowl: profiles and errors with respect to Thacker's exact solution.

Runs the executable on successively refined grids with a time step size proportional to the cell
size. The exact solution is computed here independently of the DuMux implementation. With --second-order,
it also runs the executable with the second-order flux, with implicit Euler and with BDF2, and compares
the schemes.

Usage:
  python3 plot_results.py <executable> [--second-order <executable>] [output directory of the figures]
"""

import math
import os
import re
import subprocess
import sys
import xml.etree.ElementTree as ET

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
from plotstyle import EXACT, plt, series

CELLS = [40, 80, 160, 320]
SECOND_ORDER_CELLS = [40, 80, 160]
SECOND_ORDER_SCHEMES = {"implicit Euler": ("bowl_secondorder", "ImplicitEuler"), "BDF2": ("bowl_secondorder_bdf2", "BDF2")}
# time step size times number of cells per direction
TIME_STEP_SIZE_TIMES_CELLS = 4.0
PROFILE_TIMES = {"surface": 13.0, "discharge": 11.0}
SPLIT_TIMES = [7.0, 13.0]
BAND = 0.2

ERRORS = re.compile(r"At t = (\S+) s:\nErrors \(h, u, v\) on \d+ elements:\n\s*L1: (\S+) \S+ \S+\n\s*L2: \S+ \S+ \S+\n"
                    r"Norms of the exact solution \(h, u, v\):\n\s*L1: (\S+) \S+ \S+")


def parameters(filename="params.input"):
    values = dict(re.findall(r"^\s*(\w+)\s*=\s*([^#\n]+)", open(filename).read(), re.MULTILINE))
    return (float(values["BowlDepthAtCenter"]), float(values["BowlParaboloidRadius"]),
            float(values["BowlInitialWaterElevationAtCenter"]), float(values["Gravity"]))


def thacker(x, y, t, D0, L, eta0, g):
    """Water depth and velocity of the radially symmetric oscillation (Thacker 1981, Sec. 4)"""
    omega = math.sqrt(8.0*g*D0)/L
    A = ((D0 + eta0)**2 - D0**2)/((D0 + eta0)**2 + D0**2)
    r2 = (np.asarray(x)**2 + np.asarray(y)**2)/L**2
    c = 1.0 - A*math.cos(omega*t)
    depth = D0*(math.sqrt(1.0 - A**2)/c - r2*(1.0 - A**2)/c**2)
    factor = 0.5*omega*A*math.sin(omega*t)/c
    wet = depth > 0.0
    return np.where(wet, depth, 0.0), np.where(wet, factor*x, 0.0), np.where(wet, factor*y, 0.0)


def shoreline_radius(t, D0, L, eta0, g):
    omega = math.sqrt(8.0*g*D0)/L
    A = ((D0 + eta0)**2 - D0**2)/((D0 + eta0)**2 + D0**2)
    return L*math.sqrt((1.0 - A*math.cos(omega*t))/math.sqrt(1.0 - A**2))


def read_cells(filename):
    """Cell centers and cell data of an ASCII VTU file"""
    piece = ET.parse(filename).getroot().find("UnstructuredGrid/Piece")
    values = lambda array: np.array(array.text.split(), dtype=float)
    points = values(piece.find("Points/DataArray")).reshape(-1, 3)
    cells = {array.get("Name"): array for array in piece.find("Cells")}
    connectivity = values(cells["connectivity"]).astype(int)
    offsets = values(cells["offsets"]).astype(int)
    centers = np.array([points[connectivity[b:e]].mean(axis=0)
                        for b, e in zip(np.concatenate([[0], offsets[:-1]]), offsets)])
    return centers, {array.get("Name"): values(array) for array in piece.find("CellData")}


def vtu_at(name, time):
    """The output file of the given time step"""
    entries = re.findall(r'timestep="([^"]+)"[^>]*file="([^"]+)"', open(f"{name}.pvd").read())
    t, filename = min(((float(t), f) for t, f in entries), key=lambda e: abs(e[0] - time))
    if abs(t - time) > 1e-6:
        raise RuntimeError(f"No output at t = {time} s in {name}.pvd")
    return filename


def run(executable, cells, timeStepReduction=1, prefix="bowl", args=()):
    """Run the bowl and return the times and the relative L1 errors of the water depth"""
    name = f"{prefix}_{cells}" if timeStepReduction == 1 else f"{prefix}_{cells}_dt{timeStepReduction}"
    dt = TIME_STEP_SIZE_TIMES_CELLS/cells/timeStepReduction
    output = subprocess.run(
        [executable, "params.input", *args, "-Problem.Name", name, "-Grid.Cells0", str(cells), "-Grid.Cells1", str(cells),
         "-TimeLoop.MaxTimeStepSize", str(dt), "-TimeLoop.DtInitial", str(dt), "-Vtk.Precision", "Float64"],
        check=True, capture_output=True, text=True,
    ).stdout
    with open(f"{name}.log", "w") as log:
        log.write(output)
    t, l1, exactL1 = np.array(ERRORS.findall(output), dtype=float).T
    if not np.all(np.isfinite(l1)):
        raise RuntimeError(f"Non-finite errors for {cells} cells")
    return name, t, l1/exactL1


def row_values(name, cells, time):
    """Cell centers, water depth, velocity in x-direction and bed in the row of cells next to y = 0"""
    centers, data = read_cells(vtu_at(name, time))
    row = np.abs(centers[:, 1] - 2.0/cells) < 1e-9
    order = np.argsort(centers[row, 0])
    return centers[row][order], data["waterDepth"][row][order], data["velocityX"][row][order], data["bedSurface"][row][order]


def peak_discharge_ratio(name, cells, time, D0, L, eta0, g):
    """Largest discharge in the row next to y = 0 relative to the exact one"""
    centers, h, u, _ = row_values(name, cells, time)
    hExact, uExact, _ = thacker(centers[:, 0], centers[:, 1], time, D0, L, eta0, g)
    return np.abs(h*u).max()/np.abs(hExact*uExact).max()


def executable_path(name):
    return name if "/" in name else "./" + name


def main():
    args = sys.argv[1:]
    secondOrder = None
    if "--second-order" in args:
        i = args.index("--second-order")
        secondOrder = executable_path(args[i + 1])
        del args[i:i + 2]
    if not args:
        sys.exit(__doc__)
    executable = executable_path(args[0])
    imageDir = args[1] if len(args) > 1 else "."

    D0, L, eta0, g = parameters()
    period = 2.0*math.pi/(math.sqrt(8.0*g*D0)/L)
    print(f"Oscillation period T = {period:.6f} s\n")

    results = [run(executable, cells) for cells in CELLS]

    print("| cells | time step size (s) | rel. L1 error h, t = T | rate | rel. L1 error h, t = 2T | rate |")
    print("|--:|--:|--:|--:|--:|--:|")
    previous = None
    for cells, (name, t, error) in zip(CELLS, results):
        current = (error[np.argmin(np.abs(t - period))], error[-1])
        rate = [f"{math.log(p/c, 2):.2f}" for p, c in zip(previous, current)] if previous else ["-", "-"]
        print(f"| {cells} x {cells} | {TIME_STEP_SIZE_TIMES_CELLS/cells:.4f} "
              f"| {current[0]:.3e} | {rate[0]} | {current[1]:.3e} | {rate[1]} |")
        previous = current

    # damping of the oscillation: largest discharge in the row next to y = 0 relative to the exact one
    time = PROFILE_TIMES["discharge"]
    print(f"\n| cells | time step size (s) | peak discharge / exact peak discharge at t = {time:g} s | logarithm |")
    print("|--:|--:|--:|--:|")
    reducedTimeStep = run(executable, CELLS[0], timeStepReduction=8)
    for cells, timeStepReduction, (name, *_) in [*((c, 1, r) for c, r in zip(CELLS, results)),
                                                (CELLS[0], 8, reducedTimeStep)]:
        ratio = peak_discharge_ratio(name, cells, time, D0, L, eta0, g)
        print(f"| {cells} x {cells} | {TIME_STEP_SIZE_TIMES_CELLS/cells/timeStepReduction:.4f} "
              f"| {ratio:.3f} | {math.log(ratio):.2f} |")

    # error near the shoreline, in the wet interior and on the dry bed
    print(f"\n| t (s) | cells | within {BAND} m of the shoreline | rate | wet interior | rate | dry bed | rate |")
    print("|--:|--:|--:|--:|--:|--:|--:|--:|")
    for time in SPLIT_TIMES:
        previous = None
        radius = shoreline_radius(time, D0, L, eta0, g)
        for cells, (name, *_) in zip(CELLS, results):
            centers, data = read_cells(vtu_at(name, time))
            hExact, _, _ = thacker(centers[:, 0], centers[:, 1], time, D0, L, eta0, g)
            r = np.hypot(centers[:, 0], centers[:, 1])
            contribution = np.abs(data["waterDepth"] - hExact)/np.abs(hExact).sum()
            current = [contribution[np.abs(r - radius) < BAND].sum(), contribution[r <= radius - BAND].sum(),
                       contribution[r >= radius + BAND].sum()]
            rate = [f"{math.log(p/c, 2):.2f}" for p, c in zip(previous, current)] if previous else ["-"]*3
            print(f"| {time:g} | {cells} x {cells} | " + " | ".join(f"{c:.3e} | {q}" for c, q in zip(current, rate)) + " |")
            previous = current

    fig, (axS, axQ, axE) = plt.subplots(1, 3, figsize=(12, 4.6), layout="constrained")
    xFine = np.linspace(-2.0, 2.0, 2001)
    bed = D0*xFine**2/L**2
    for ax, quantity in [(axS, "surface"), (axQ, "discharge")]:
        time = PROFILE_TIMES[quantity]
        h, u, _ = thacker(xFine, 0.0*xFine, time, D0, L, eta0, g)
        if quantity == "surface":
            ax.fill_between(xFine, bed, -0.002, color="0.85", lw=0, label="bed")
            ax.plot(xFine, h + bed, label="exact", **EXACT)
        else:
            ax.plot(xFine, h*u, label="exact", **EXACT)
        for i, (cells, (name, *_)) in enumerate(zip(CELLS, results)):
            centers, hs, us, zs = row_values(name, cells, time)
            value = hs + zs if quantity == "surface" else hs*us
            ax.plot(centers[:, 0], value, label=f"{cells} x {cells}", **series(i, cells))
        ax.set_xlabel("$x$ (m)")
        ax.set_title(f"$t = {time:g}$ s")
    axS.set_ylim(-0.002, 0.05)
    axS.set_xlim(-2.0, 2.0)
    axS.set_ylabel(r"free surface $h + z_b$ (m)")
    axQ.set_ylabel(r"discharge $h u$ (m$^2$/s)")
    handles, labels = axS.get_legend_handles_labels()
    fig.legend(handles, labels, loc="outside lower center", ncol=len(labels))

    for i, (cells, (name, t, error)) in enumerate(zip(CELLS, results)):
        axE.semilogy(t/period, error, label=f"{cells} x {cells}", **series(i, len(t)))
    axE.set_xlabel("$t / T$")
    axE.set_ylabel(r"relative $L^1$ error of $h$")
    fig.savefig(f"{imageDir}/shallowwater_bowl.png")

    if secondOrder:
        compare_second_order(secondOrder, results, reducedTimeStep, imageDir, D0, L, eta0, g, period)


def compare_second_order(executable, firstOrder, firstOrderReducedTimeStep, imageDir, D0, L, eta0, g, period):
    """Errors and damping of the second-order flux with implicit Euler and BDF2, compared to the first-order flux"""
    results, reduced = {}, {}
    for scheme, (prefix, parameter) in SECOND_ORDER_SCHEMES.items():
        args = ["-TimeLoop.Scheme", parameter]
        results[scheme] = [run(executable, cells, prefix=prefix, args=args) for cells in SECOND_ORDER_CELLS]
        reduced[scheme] = run(executable, SECOND_ORDER_CELLS[0], timeStepReduction=8, prefix=prefix, args=args)

        print(f"\nSecond-order flux, {scheme}:\n")
        print("| cells | time step size (s) | rel. L1 error h, t = T | rate | rel. L1 error h, t = 2T | rate |")
        print("|--:|--:|--:|--:|--:|--:|")
        previous = None
        for cells, (name, t, error) in zip(SECOND_ORDER_CELLS, results[scheme]):
            current = (error[np.argmin(np.abs(t - period))], error[-1])
            rate = [f"{math.log(p/c, 2):.2f}" for p, c in zip(previous, current)] if previous else ["-", "-"]
            firstOrderSteps = re.findall(r"Time steps taken with BDF1: (\d+)", open(f"{name}.log").read())
            print(f"| {cells} x {cells} | {TIME_STEP_SIZE_TIMES_CELLS/cells:.4f} "
                  f"| {current[0]:.3e} | {rate[0]} | {current[1]:.3e} | {rate[1]} |"
                  + (f" BDF1 steps: {firstOrderSteps[0]}" if firstOrderSteps else ""))
            previous = current

    time = PROFILE_TIMES["discharge"]
    print(f"\n| cells | time step size (s) | first order | second order, implicit Euler | second order, BDF2 |")
    print("|--:|--:|--:|--:|--:|")
    rows = [(cells, 1, firstOrder[i][0], {s: results[s][i][0] for s in results})
            for i, cells in enumerate(SECOND_ORDER_CELLS)]
    rows.append((SECOND_ORDER_CELLS[0], 8, firstOrderReducedTimeStep[0], {s: reduced[s][0] for s in reduced}))
    for cells, timeStepReduction, firstOrderName, names in rows:
        ratios = [peak_discharge_ratio(n, cells, time, D0, L, eta0, g) for n in [firstOrderName, *names.values()]]
        print(f"| {cells} x {cells} | {TIME_STEP_SIZE_TIMES_CELLS/cells/timeStepReduction:.4f} | "
              + " | ".join(f"{r:.3f}" for r in ratios) + " |")

    cells = SECOND_ORDER_CELLS[-1]
    runs = [("first order", firstOrder[CELLS.index(cells)])]
    runs += [(f"second order, {scheme}", results[scheme][-1]) for scheme in results]
    fig, (axS, axQ, axE) = plt.subplots(1, 3, figsize=(12, 4.6), layout="constrained")
    xFine = np.linspace(-2.0, 2.0, 2001)
    bed = D0*xFine**2/L**2
    for ax, quantity in [(axS, "surface"), (axQ, "discharge")]:
        time = PROFILE_TIMES[quantity]
        h, u, _ = thacker(xFine, 0.0*xFine, time, D0, L, eta0, g)
        if quantity == "surface":
            ax.fill_between(xFine, bed, -0.002, color="0.85", lw=0, label="bed")
            ax.plot(xFine, h + bed, label="exact", **EXACT)
        else:
            ax.plot(xFine, h*u, label="exact", **EXACT)
        for i, (label, (name, *_)) in enumerate(runs):
            centers, hs, us, zs = row_values(name, cells, time)
            value = hs + zs if quantity == "surface" else hs*us
            ax.plot(centers[:, 0], value, label=label, **series(i, cells))
        ax.set_xlabel("$x$ (m)")
        ax.set_title(f"{cells} x {cells}, $t = {time:g}$ s")
    axS.set_ylim(-0.002, 0.05)
    axS.set_xlim(-2.0, 2.0)
    axS.set_ylabel(r"free surface $h + z_b$ (m)")
    axQ.set_ylabel(r"discharge $h u$ (m$^2$/s)")
    handles, labels = axS.get_legend_handles_labels()
    fig.legend(handles, labels, loc="outside lower center", ncol=len(labels))

    for i, (label, (name, t, error)) in enumerate(runs):
        axE.semilogy(t/period, error, label=label, **series(i, len(t)))
    axE.set_xlabel("$t / T$")
    axE.set_ylabel(r"relative $L^1$ error of $h$")
    fig.savefig(f"{imageDir}/shallowwater_bowl_secondorder.png")


if __name__ == "__main__":
    main()
