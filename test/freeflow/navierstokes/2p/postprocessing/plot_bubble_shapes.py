#!/usr/bin/env python3
# SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later
"""
Overlay the phase-field zero level set (phi = 0, the bubble interface) at a sequence of times
in ONE figure - a shape-evolution diagnostic for the rising-bubble CH-NS test.

The simulation uses the LEFT half domain [0,0.5]x[0,2] with a symmetry plane at x=0.5 (the
bubble centre). Each contour is therefore MIRRORED about x=0.5 to show the full bubble.

Usage:
    plot_bubble_shapes.py <run.pvd> <out.png> "<title>" [t0 t1 t2 ...]
If no times are given, defaults to 0,0.5,...,3.0. Colour encodes time (viridis).

Requires: matplotlib, numpy, meshio.
"""
import sys, os, xml.etree.ElementTree as ET
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.tri as mtri
from matplotlib.cm import ScalarMappable
from matplotlib.colors import Normalize
import meshio

SYM_X = 0.5  # symmetry plane / mirror axis

def pvd_index(pvd):
    """Return sorted [(time, vtu_path), ...] from a ParaView .pvd collection."""
    base = os.path.dirname(os.path.abspath(pvd))
    out = []
    for ds in ET.parse(pvd).getroot().iter("DataSet"):
        out.append((float(ds.get("timestep")), os.path.join(base, ds.get("file"))))
    return sorted(out)

def nearest(index, t):
    return min(index, key=lambda it: abs(it[0] - t))

def contour_segments(vtu, level=0.0):
    """phi=level segments as a list of (N,2) arrays, in the simulated half-domain."""
    m = meshio.read(vtu)
    x, y = m.points[:, 0], m.points[:, 1]
    phi = np.asarray(m.point_data["phi"]).ravel()
    tris = m.cells_dict.get("triangle")
    if tris is None:
        return []
    tri = mtri.Triangulation(x, y, tris)
    # extract on a THROWAWAY figure so the main axes is never touched/cleared
    figx, axx = plt.subplots()
    cs = axx.tricontour(tri, phi, levels=[level])
    segs = [np.asarray(seg) for seg in cs.allsegs[0] if len(seg) > 1]
    plt.close(figx)
    return segs

def main():
    pvd, out, title = sys.argv[1], sys.argv[2], sys.argv[3]
    times = [float(a) for a in sys.argv[4:]] or [0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0]
    index = pvd_index(pvd)

    fig, ax = plt.subplots(figsize=(5.2, 8.4))
    norm = Normalize(vmin=min(times), vmax=max(times))
    cmap = plt.get_cmap("viridis")

    for t in times:
        tt, vtu = nearest(index, t)
        col = cmap(norm(t))
        for seg in contour_segments(vtu):
            xs, ys = seg[:, 0], seg[:, 1]
            ax.plot(xs, ys, color=col, lw=1.6)                 # simulated half
            ax.plot(2 * SYM_X - xs, ys, color=col, lw=1.6)     # mirror about x=0.5
    ax.axvline(SYM_X, color="0.8", lw=0.7, ls=":")             # symmetry plane

    ax.set_aspect("equal")
    ax.set_xlim(0.15, 0.85)
    ax.set_xlabel("x"); ax.set_ylabel("y")
    ax.set_title(title, fontsize=11)
    ax.grid(alpha=0.25)
    sm = ScalarMappable(norm=norm, cmap=cmap); sm.set_array([])
    cb = fig.colorbar(sm, ax=ax, fraction=0.045, pad=0.02); cb.set_label("time")
    fig.tight_layout()
    fig.savefig(out, dpi=140, bbox_inches="tight")
    print("wrote", out, "| times:", ", ".join(f"{t:g}" for t in times))

if __name__ == "__main__":
    main()
