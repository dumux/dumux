#!/usr/bin/env python3
# SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later
"""Wet dam break: profiles at the final time and relative L1 errors under grid refinement.

The exact solution (Stoker) is computed here independently of the DuMux implementation. For
comparison, the script also solves the problem with a textbook first-order Godunov scheme
(exact Riemann solver, explicit Euler) with the same cell and time step sizes. With --second-order,
it also runs the executable with the second-order flux and compares both schemes.

Usage:
  python3 plot_results.py <executable> [--second-order <executable>] [output directory of the figures]
"""

import math
import os
import sys
import xml.etree.ElementTree as ET

import numpy as np
from scipy.optimize import brentq

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
from plotstyle import EXACT, cell_ticks, plt, series, slope_triangle
from convergencetest import TIME_STEP_SIZE_100, rates, run

CELLS = [100, 200, 400, 800, 1600]
PROFILE_CELLS = [100, 400, 1600]
DEPTH_LEFT, DEPTH_RIGHT, GATE, LENGTH, GRAVITY, TIME = 4.0, 1.0, 10.0, 20.0, 9.81, 1.0


def stoker(x, t):
    """Depth and velocity of the dam break on a wet bed, left rarefaction and right shock"""
    cl = math.sqrt(GRAVITY*DEPTH_LEFT)

    def mismatch(hm):
        shock = (hm - DEPTH_RIGHT)*math.sqrt(0.5*GRAVITY*(hm + DEPTH_RIGHT)/(hm*DEPTH_RIGHT))
        return 2.0*(cl - math.sqrt(GRAVITY*hm)) - shock

    hm = brentq(mismatch, DEPTH_RIGHT, DEPTH_LEFT, xtol=1e-14)
    um = 2.0*(cl - math.sqrt(GRAVITY*hm))
    shockSpeed = hm*um/(hm - DEPTH_RIGHT)
    xi = (np.asarray(x) - GATE)/t
    fan = (xi > -cl) & (xi <= um - math.sqrt(GRAVITY*hm))
    h = np.where(xi <= -cl, DEPTH_LEFT, np.where(fan, (2.0*cl - xi)**2/(9.0*GRAVITY),
                 np.where(xi <= shockSpeed, hm, DEPTH_RIGHT)))
    u = np.where(fan, 2.0/3.0*(cl + xi), np.where((xi > -cl) & (xi <= shockSpeed), um, 0.0))
    return h, u, hm, um, shockSpeed


def godunov_flux(hl, ul, hr, ur):
    """Flux of the exact solution of the Riemann problem at x/t = 0 for wet states"""
    g = GRAVITY
    cl, cr = np.sqrt(g*hl), np.sqrt(g*hr)
    f = lambda h, hk, ck: np.where(h <= hk, 2.0*(np.sqrt(g*h) - ck), (h - hk)*np.sqrt(0.5*g*(h + hk)/(h*hk)))
    df = lambda h, hk: np.where(h <= hk, np.sqrt(g/h), np.sqrt(0.5*g*(h + hk)/(h*hk))
                                - 0.25*g*(h - hk)/(np.sqrt(0.5*g*(h + hk)/(h*hk))*h*h))
    hs = (0.5*(cl + cr) - 0.25*(ur - ul))**2/g
    for _ in range(30):
        hs = np.maximum(hs - (f(hs, hl, cl) + f(hs, hr, cr) + ur - ul)/(df(hs, hl) + df(hs, hr)), 1e-12)
    us = 0.5*(ul + ur) + 0.5*(f(hs, hr, cr) - f(hs, hl, cl))
    cs = np.sqrt(g*hs)
    h, u = hs.copy(), us.copy()
    left = us >= 0.0
    keepL = left & np.where(hs > hl, ul - cl*np.sqrt(0.5*(hs + hl)*hs)/hl >= 0.0, ul - cl >= 0.0)
    fanL = left & (hs <= hl) & (ul - cl < 0.0) & (us - cs > 0.0)
    keepR = ~left & np.where(hs > hr, ur + cr*np.sqrt(0.5*(hs + hr)*hs)/hr <= 0.0, ur + cr <= 0.0)
    fanR = ~left & (hs <= hr) & (ur + cr > 0.0) & (us + cs < 0.0)
    h = np.select([keepL, fanL, keepR, fanR], [hl, ((ul + 2*cl)/3)**2/g, hr, ((-ur + 2*cr)/3)**2/g], h)
    u = np.select([keepL, fanL, keepR, fanR], [ul, (ul + 2*cl)/3, ur, (ur - 2*cr)/3], u)
    return h*u, h*u*u + 0.5*g*h*h


def godunov(cells, dt):
    """First-order Godunov scheme with explicit Euler and reflecting walls"""
    dx = LENGTH/cells
    x = (np.arange(cells) + 0.5)*dx
    h, q = np.where(x < GATE, DEPTH_LEFT, DEPTH_RIGHT), np.zeros(cells)
    for _ in range(round(TIME/dt)):
        u = q/h
        hb, ub = np.concatenate([[h[0]], h, [h[-1]]]), np.concatenate([[-u[0]], u, [-u[-1]]])
        fh, fq = godunov_flux(hb[:-1], ub[:-1], hb[1:], ub[1:])
        h, q = h - dt/dx*(fh[1:] - fh[:-1]), q - dt/dx*(fq[1:] - fq[:-1])
    return x, h


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

    results = [run(executable, cells, ["-Vtk.Precision", "Float64"]) for cells in CELLS]
    _, _, hm, um, shockSpeed = stoker(0.0, TIME)
    print(f"Stoker: h_m = {hm:.6f} m, u_m = {um:.6f} m/s, shock speed {shockSpeed:.6f} m/s\n")

    reference = []
    for cells in CELLS:
        x, h = godunov(cells, TIME_STEP_SIZE_100*100/cells)
        hExact, *_ = stoker(x, TIME)
        reference.append(np.abs(h - hExact).sum()/np.abs(hExact).sum())

    rh, ru, rg = rates([r["h"] for r in results]), rates([r["u"] for r in results]), rates(reference)
    print("| cells | time step size (s) | rel. L1 error h | rate | rel. L1 error u | rate "
          "| rel. L1 error h, explicit Godunov | rate |")
    print("|--:|--:|--:|--:|--:|--:|--:|--:|")
    for i, (cells, r) in enumerate(zip(CELLS, results)):
        rate = lambda values: f"{values[i - 1]:.2f}" if i > 0 else "-"
        print(f"| {cells} | {TIME_STEP_SIZE_100*100/cells:.3e} | {r['h']:.3e} | {rate(rh)} "
              f"| {r['u']:.3e} | {rate(ru)} | {reference[i]:.3e} | {rate(rg)} |")

    # contributions to the relative L1 error of the water depth from the vicinity of the waves
    cl = math.sqrt(GRAVITY*DEPTH_LEFT)
    waves = {"bore": GATE + shockSpeed*TIME, "head of the rarefaction": GATE - cl*TIME,
             "tail of the rarefaction": GATE + (um - math.sqrt(GRAVITY*hm))*TIME}
    print("\n| cells | max deviation of exactWaterDepth (m) | "
          + " | ".join(f"within 1 m of the {w}" for w in waves) + " | elsewhere |")
    print("|--:|--:|" + "--:|"*(len(waves) + 1))
    for cells in CELLS:
        centers, data = read_cells(f"dambreak_{cells}-00001.vtu")
        h, *_ = stoker(centers[:, 0], TIME)
        contribution = np.abs(data["waterDepth"] - h)/np.abs(h).sum()
        near = {w: np.abs(centers[:, 0] - position) < 1.0 for w, position in waves.items()}
        elsewhere = ~np.any(list(near.values()), axis=0)
        print(f"| {cells} | {np.abs(data['exactWaterDepth'] - h).max():.1e} | "
              + " | ".join(f"{contribution[m].sum():.3e}" for m in near.values())
              + f" | {contribution[elsewhere].sum():.3e} |")

    fig, (axH, axU, axE) = plt.subplots(1, 3, figsize=(12, 4.0), layout="constrained")
    xFine = np.linspace(0.0, LENGTH, 4001)
    hExact, uExact, *_ = stoker(xFine, TIME)
    axH.plot(xFine, hExact, label="exact", **EXACT)
    axU.plot(xFine, uExact, label="exact", **EXACT)
    for i, cells in enumerate(PROFILE_CELLS):
        centers, data = read_cells(f"dambreak_{cells}-00001.vtu")
        order = np.argsort(centers[:, 0])
        axH.plot(centers[order, 0], data["waterDepth"][order], label=f"{cells} cells", **series(i, cells))
        axU.plot(centers[order, 0], data["velocityX"][order], label=f"{cells} cells", **series(i, cells))
    axH.set_xlabel("$x$ (m)")
    axH.set_ylabel("water depth $h$ (m)")
    axU.set_xlabel("$x$ (m)")
    axU.set_ylabel("velocity $u$ (m/s)")
    axH.legend()

    axE.loglog(CELLS, [r["h"] for r in results], label=r"$h$, DuMux", **series(0))
    axE.loglog(CELLS, [r["u"] for r in results], label=r"$u$, DuMux", **series(1))
    axE.loglog(CELLS, reference, label=r"$h$, explicit Godunov", **series(2))
    slope_triangle(axE, CELLS[-2], CELLS[-1], 0.6*reference[-2], 1)
    cell_ticks(axE, CELLS)
    axE.set_xlabel("cells along the channel")
    axE.set_ylabel(r"relative $L^1$ error at $t = 1$ s")
    axE.legend()
    fig.savefig(f"{imageDir}/shallowwater_dambreak.png")

    if secondOrder:
        compare_second_order(secondOrder, results, imageDir)


def compare_second_order(executable, firstOrder, imageDir):
    """Errors, over- and undershoots of the second-order flux, and profiles of both schemes"""
    results = [run(executable, cells, ["-Vtk.Precision", "Float64"], name="dambreak_secondorder") for cells in CELLS]
    rh, ru = rates([r["h"] for r in results]), rates([r["u"] for r in results])
    print("\nSecond-order flux:\n")
    print("| cells | rel. L1 error h | rate | rel. L1 error u | rate | error h / first order "
          "| max h - h_L (m) | h_R - min h (m) |")
    print("|--:|--:|--:|--:|--:|--:|--:|--:|")
    for i, (cells, r, r1) in enumerate(zip(CELLS, results, firstOrder)):
        rate = lambda values: f"{values[i - 1]:.2f}" if i > 0 else "-"
        _, data = read_cells(f"dambreak_secondorder_{cells}-00001.vtu")
        print(f"| {cells} | {r['h']:.3e} | {rate(rh)} | {r['u']:.3e} | {rate(ru)} | {r['h']/r1['h']:.2f} "
              f"| {data['waterDepth'].max() - DEPTH_LEFT:.1e} | {DEPTH_RIGHT - data['waterDepth'].min():.1e} |")

    cells = PROFILE_CELLS[1]
    fig, (axH, axU, axE) = plt.subplots(1, 3, figsize=(12, 4.0), layout="constrained")
    xFine = np.linspace(0.0, LENGTH, 4001)
    hExact, uExact, *_ = stoker(xFine, TIME)
    axH.plot(xFine, hExact, label="exact", **EXACT)
    axU.plot(xFine, uExact, label="exact", **EXACT)
    for i, (name, label) in enumerate([("dambreak", "first order"), ("dambreak_secondorder", "second order")]):
        centers, data = read_cells(f"{name}_{cells}-00001.vtu")
        order = np.argsort(centers[:, 0])
        axH.plot(centers[order, 0], data["waterDepth"][order], label=label, **series(i, cells))
        axU.plot(centers[order, 0], data["velocityX"][order], label=label, **series(i, cells))
    for ax in (axH, axU):
        ax.set_xlabel("$x$ (m)")
        ax.set_title(f"{cells} cells")
    axH.set_ylabel("water depth $h$ (m)")
    axU.set_ylabel("velocity $u$ (m/s)")
    axH.legend()

    axE.loglog(CELLS, [r["h"] for r in firstOrder], label=r"$h$, first order", **series(0))
    axE.loglog(CELLS, [r["u"] for r in firstOrder], label=r"$u$, first order", **series(2))
    axE.loglog(CELLS, [r["h"] for r in results], label=r"$h$, second order", **series(1))
    axE.loglog(CELLS, [r["u"] for r in results], label=r"$u$, second order", **series(3))
    slope_triangle(axE, CELLS[-2], CELLS[-1], 0.6*results[-2]["h"], 1)
    cell_ticks(axE, CELLS)
    axE.set_xlabel("cells along the channel")
    axE.set_ylabel(r"relative $L^1$ error at $t = 1$ s")
    axE.legend(loc="upper center", bbox_to_anchor=(0.5, -0.2), ncol=2, fontsize=11)
    fig.savefig(f"{imageDir}/shallowwater_dambreak_secondorder.png")


if __name__ == "__main__":
    main()
