#!/usr/bin/env python3
# SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later
"""Draw a modern schematic of the Cryer benchmark setup and boundary conditions.

By symmetry only one octant of the sphere is modelled; the figure shows a meridional
quarter-disk cross-section with colour-coded boundary conditions:

  - applied compressive surface load p_0           (coral arrows on the arc)
  - drained sphere surface, p_f = 0                (blue arc)
  - symmetry planes: zero normal displacement      (teal roller supports)
    and no flow

The inset shows the actual 3D tetrahedral mesh (sphere_eighth.msh) used in the
simulation; if the mesh file is not found it falls back to a smooth shaded octant.

Usage:
  python3 make_setup_svg.py [--out cryer_setup.svg] [--mesh sphere_eighth.msh]
"""
import argparse
import os

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.patches import Wedge, FancyArrowPatch, Circle
from matplotlib.lines import Line2D
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
ap = argparse.ArgumentParser()
ap.add_argument("--out", default="cryer_setup.svg")
ap.add_argument("--mesh", default=os.path.join(SCRIPT_DIR, "sphere_eighth.msh"))
args = ap.parse_args()

# palette
LOAD = "#E8633C"   # mechanical load
DRAIN = "#2D7DD2"  # drained surface
SYMM = "#149E96"   # symmetry planes
INK = "#22303B"
# plain sans text everywhere (no mathtext); Unicode subscripts where available
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 11})

R = 1.0
fig, ax = plt.subplots(figsize=(8.4, 7.0))
ax.set_aspect("equal")
ax.set_axis_off()

# ---- body: radial gradient (suggests the sphere), clipped to the quarter disk ----
n = 500
gx, gy = np.meshgrid(np.linspace(0, R, n), np.linspace(0, R, n))
rad = np.sqrt(gx**2 + gy**2) / R
cmap = LinearSegmentedColormap.from_list("body", ["#eef6fd", "#bcdcf6", "#8fc1ec"])
im = ax.imshow(np.clip(rad, 0, 1), extent=[0, R, 0, R], origin="lower",
               cmap=cmap, vmin=0, vmax=1.05, zorder=0, interpolation="bilinear")
clip = Wedge((0, 0), R, 0, 90, width=R, transform=ax.transData,
             facecolor="none", edgecolor="none")
ax.add_patch(clip)
im.set_clip_path(clip)

a = np.linspace(0, np.pi / 2, 200)
# soft shadow under the body
ax.fill(np.r_[0, R * np.cos(a), 0] + 0.012, np.r_[0, R * np.sin(a), 0] - 0.012,
        color="#000000", alpha=0.06, zorder=-1)

# ---- drained spherical surface (blue arc) ----
ax.plot(R * np.cos(a), R * np.sin(a), color=DRAIN, lw=3.2, zorder=3,
        solid_capstyle="round")

# ---- applied load: coral arrows pressing onto the arc ----
D = 0.20
angs = np.linspace(0.06, np.pi / 2 - 0.06, 11)
ax.plot((R + D) * np.cos(angs), (R + D) * np.sin(angs), color=LOAD, lw=1.0,
        alpha=0.45, zorder=2)
for th in angs:
    nx, ny = np.cos(th), np.sin(th)
    ax.add_patch(FancyArrowPatch(((R + D) * nx, (R + D) * ny), (R * nx, R * ny),
                 arrowstyle="-|>", mutation_scale=13, lw=1.8, color=LOAD,
                 shrinkA=0, shrinkB=0, zorder=4))

# ---- symmetry planes: teal edges + roller supports ----
ax.plot([0, R], [0, 0], color=SYMM, lw=3.0, zorder=3, solid_capstyle="round")
ax.plot([0, 0], [0, R], color=SYMM, lw=3.0, zorder=3, solid_capstyle="round")


def rollers(fixed, offset):
    """Row of roller circles + support rail just outside a symmetry edge."""
    r = 0.028
    pos = np.linspace(0.22, 0.86, 4) * R
    for s in pos:
        if offset[0]:   # left edge (x = 0), rollers to the left
            c = (fixed + offset[0] * (r + 0.012), s)
        else:           # bottom edge (y = 0), rollers below
            c = (s, fixed + offset[1] * (r + 0.012))
        ax.add_patch(Circle(c, r, facecolor="white", edgecolor=SYMM, lw=1.4, zorder=4))
    if offset[0]:
        xr = fixed + offset[0] * (2 * r + 0.024)
        ax.plot([xr, xr], [pos.min() - 0.04, pos.max() + 0.04], color=SYMM, lw=1.6, zorder=4)
    else:
        yr = fixed + offset[1] * (2 * r + 0.024)
        ax.plot([pos.min() - 0.04, pos.max() + 0.04], [yr, yr], color=SYMM, lw=1.6, zorder=4)


rollers(0.0, (-1, 0))   # left edge
rollers(0.0, (0, -1))   # bottom edge

# ---- centre marker ----
ax.add_patch(Circle((0, 0), 0.022, facecolor=INK, edgecolor="white", lw=1.0, zorder=6))

# ---- labels ----
ax.text(0.84 * R, 1.12 * R, "p₀", color=LOAD, fontsize=22, fontweight="bold",
        ha="center", va="center", zorder=6)
ax.annotate("drained surface\np = 0", xy=(R * np.cos(0.62), R * np.sin(0.62)),
            xytext=(1.02 * R, 0.30 * R), color=DRAIN, fontsize=11, ha="left", va="center",
            arrowprops=dict(arrowstyle="-", color=DRAIN, lw=1.0, alpha=0.7))
ax.text(0.52 * R, -0.135, "symmetry  ·  uₙ = 0  ·  no flow", color=SYMM,
        fontsize=10.5, ha="center", va="top")

# ---- legend ----
handles = [
    Line2D([0], [0], color=LOAD, lw=2.4, marker=">", markersize=7, label="applied load  p₀"),
    Line2D([0], [0], color=DRAIN, lw=3.0, label="drained surface  (p = 0)"),
    Line2D([0], [0], color=SYMM, lw=3.0, marker="o", markerfacecolor="white",
           markersize=7, label="symmetry:  uₙ = 0,  no flow"),
]
ax.legend(handles=handles, loc="lower right", bbox_to_anchor=(1.16, -0.02),
          frameon=False, fontsize=10.5, handlelength=1.8)

ax.set_xlim(-0.20, R + 0.55)
ax.set_ylim(-0.26, R + 0.30)

fig.tight_layout()


# ---- inset: the actual 3D mesh used in the simulation -----------------------
def load_msh_surface(path):
    """Parse a gmsh msh2 file -> (nodes dict, list of (i,j,k,phys) triangles)."""
    nodes, tris = {}, []
    with open(path) as f:
        lines = f.read().splitlines()
    i = 0
    while i < len(lines):
        ln = lines[i]
        if ln.startswith("$Nodes"):
            cnt = int(lines[i + 1]); i += 2
            for k in range(cnt):
                p = lines[i + k].split()
                nodes[int(p[0])] = (float(p[1]), float(p[2]), float(p[3]))
            i += cnt
        elif ln.startswith("$Elements"):
            cnt = int(lines[i + 1]); i += 2
            for k in range(cnt):
                p = lines[i + k].split()
                if int(p[1]) == 2:                      # 3-node triangle
                    ntags = int(p[2])
                    phys = int(p[3]) if ntags >= 1 else 0
                    a3 = list(map(int, p[3 + ntags:3 + ntags + 3]))
                    tris.append((a3[0], a3[1], a3[2], phys))
            i += cnt
        else:
            i += 1
    return nodes, tris


iax = fig.add_axes([0.63, 0.62, 0.30, 0.32], projection="3d")
iax.set_axis_off()
iax.view_init(elev=20, azim=40)
try:
    nodes, tris = load_msh_surface(args.mesh)
    sphere_polys, sphere_norm, flat_polys = [], [], []
    for (a3, b3, c3, phys) in tris:
        tri = [nodes[a3], nodes[b3], nodes[c3]]
        if phys == 4:                                   # sphere surface
            sphere_polys.append(tri)
            sphere_norm.append(np.mean(tri, axis=0))
        else:                                           # symmetry planes
            flat_polys.append(tri)
    # translucent so the back of the mesh shines through
    iax.add_collection3d(Poly3DCollection(
        flat_polys, facecolor="#dfeee9", edgecolor=SYMM, linewidths=0.18, alpha=0.32))
    iax.add_collection3d(Poly3DCollection(
        sphere_polys, facecolor="#bcdcf6", edgecolor=DRAIN, linewidths=0.2, alpha=0.4))
    # a few load arrows on the meshed surface
    for nrm in sphere_norm[:: max(1, len(sphere_norm) // 7)]:
        nrm = np.array(nrm); u = nrm / np.linalg.norm(nrm)
        tail, vec = u * 1.28, -0.26 * u
        iax.quiver(*tail, *vec, color=LOAD, lw=1.3, arrow_length_ratio=0.5)
    iax.text2D(0.5, -0.02, "3D mesh (one octant)", transform=iax.transAxes,
               ha="center", color=INK, fontsize=9)
except (OSError, KeyError, ValueError):
    # fallback: smooth shaded octant sphere
    t = np.linspace(0, np.pi / 2, 30); p = np.linspace(0, np.pi / 2, 30)
    T, P = np.meshgrid(t, p)
    iax.plot_surface(np.sin(T) * np.cos(P), np.sin(T) * np.sin(P), np.cos(T),
                     color="#bcdcf6", alpha=0.9, linewidth=0, antialiased=True)
    iax.text2D(0.5, -0.02, "one octant", transform=iax.transAxes,
               ha="center", color=INK, fontsize=9)
iax.set_xlim(0, 1); iax.set_ylim(0, 1); iax.set_zlim(0, 1)
iax.set_box_aspect((1, 1, 1))

out_svg = args.out
out_png = args.out.rsplit(".", 1)[0] + ".png"
fig.savefig(out_svg, bbox_inches="tight")
fig.savefig(out_png, dpi=150, bbox_inches="tight")
plt.close(fig)
print(f"Saved {out_svg} and {out_png}")
