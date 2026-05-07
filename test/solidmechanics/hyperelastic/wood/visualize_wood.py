#!/usr/bin/env python3
# SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

import argparse
import configparser
import math
import multiprocessing as mp
import os
import pathlib
import sys
import xml.etree.ElementTree as ET

import numpy as np
import pyvista as pv

EARLYWOOD = np.array([0.82, 0.63, 0.38], dtype=np.float32)
LATEWOOD = np.array([0.46, 0.28, 0.12], dtype=np.float32)
SECTION_LINE = np.array([0.12, 0.08, 0.04], dtype=np.float32)
WET_HUE = np.array([0.15, 0.55, 0.95], dtype=np.float32)


def read_pvd(path: pathlib.Path):
    tree = ET.parse(path)
    entries = []
    for ds in tree.getroot().iter("DataSet"):
        t = float(ds.attrib.get("timestep", 0.0))
        entries.append((t, path.parent / ds.attrib["file"]))
    entries.sort(key=lambda e: e[0])
    return entries


def read_pith(path: pathlib.Path) -> np.ndarray:
    cfg = configparser.ConfigParser()
    cfg.read(str(path))
    return np.array([float(cfg["SpatialParams"]["PithX"]),
                     float(cfg["SpatialParams"]["PithY"])], dtype=float)


def wood_colors(r: np.ndarray, t: np.ndarray, moisture,
                ring_width: float, section_spacing: float) -> np.ndarray:
    # Annual rings: smooth earlywood -> latewood transition within each ring.
    ew_frac = 0.6
    phase = np.mod(r / ring_width, 1.0)
    s = np.where(
        phase < ew_frac,
        0.5 * (1.0 - np.cos(math.pi * phase / ew_frac)),
        0.5 * (1.0 + np.cos(math.pi * (phase - ew_frac) / (1.0 - ew_frac))),
    ).astype(np.float32)[:, None]
    colors = (1.0 - s) * EARLYWOOD + s * LATEWOOD

    # Radial section lines (constant angle from pith, perpendicular to rings).
    sec_phase = np.mod(t / section_spacing, 1.0)
    centered = np.minimum(sec_phase, 1.0 - sec_phase)
    line = np.power(np.exp(-np.square(centered / 0.12)), 1.6).astype(np.float32)
    line = np.clip(0.4 * line, 0.0, 1.0)[:, None]
    colors = (1.0 - line) * colors + line * SECTION_LINE

    if moisture is not None:
        wet = np.clip((np.asarray(moisture).reshape(-1) - 0.05) / 0.25, 0.0, 1.0)[:, None]
        colors = (1.0 - wet) * colors + wet * WET_HUE
    return np.clip(colors, 0.0, 1.0)


def build_surface(mesh: pv.UnstructuredGrid, pith: np.ndarray,
                  warp_scale: float, ring_width: float,
                  section_spacing: float, subdivisions: int) -> pv.PolyData:
    ref_xy = mesh.points[:, :2].copy()
    d2 = np.asarray(mesh.point_data["d"])
    d3 = np.column_stack([d2[:, :2], np.zeros(len(d2))])
    mesh.point_data["displacement"] = d3
    warped = mesh.warp_by_vector("displacement", factor=warp_scale)
    warped.point_data["ReferencePoints"] = np.column_stack([ref_xy, np.zeros(len(ref_xy))])

    surface = warped.extract_surface().triangulate()
    if subdivisions > 0:
        surface = surface.subdivide(subdivisions, subfilter="linear")

    delta = surface.point_data["ReferencePoints"][:, :2] - pith
    r = np.linalg.norm(delta, axis=1)
    theta = np.arctan2(delta[:, 1], delta[:, 0])
    theta_ref = np.arctan2(delta[:, 1].mean(), delta[:, 0].mean())
    # Iso-contours of t are radial rays (constant theta from pith).
    dtheta = np.mod(theta - theta_ref + math.pi, 2.0 * math.pi) - math.pi
    t_coord = float(r.mean()) * dtheta

    moisture = surface.point_data.get("m", None)
    surface.point_data["WoodColor"] = wood_colors(
        r, t_coord, moisture, ring_width, section_spacing
    )
    return surface


def interpolate_mesh(a: pv.UnstructuredGrid, b: pv.UnstructuredGrid,
                     alpha: float) -> pv.UnstructuredGrid:
    if a.n_points != b.n_points or a.n_cells != b.n_cells:
        return a
    out = a.copy(deep=True)
    out.points = (1.0 - alpha) * a.points + alpha * b.points
    for name in set(a.point_data) & set(b.point_data):
        x, y = np.asarray(a.point_data[name]), np.asarray(b.point_data[name])
        if x.shape == y.shape and np.issubdtype(x.dtype, np.number):
            out.point_data[name] = (1.0 - alpha) * x + alpha * y
    return out


def time_brackets(times: np.ndarray, num_frames: int):
    """Yield (target_time, i0, i1, alpha) for uniform target times."""
    targets = np.linspace(times[0], times[-1], num_frames)
    for t in targets:
        i = int(np.searchsorted(times, t, side="left"))
        if i <= 0:
            yield float(t), 0, 0, 0.0
        elif i >= len(times):
            last = len(times) - 1
            yield float(t), last, last, 0.0
        else:
            dt = float(times[i] - times[i - 1])
            alpha = float((t - times[i - 1]) / dt) if dt > 0.0 else 0.0
            yield float(t), i - 1, i, alpha


def render_worker(spec: dict):
    a = pv.read(spec["vtu_a"])
    if spec["alpha"] > 1e-12 and spec["vtu_a"] != spec["vtu_b"]:
        mesh = interpolate_mesh(a, pv.read(spec["vtu_b"]), spec["alpha"])
    else:
        mesh = a
    pith = np.asarray(spec["pith"], dtype=float)
    surface = build_surface(
        mesh, pith, spec["warp_scale"], spec["ring_width"],
        spec["section_spacing"], spec["subdivisions"],
    )
    plotter = pv.Plotter(off_screen=True, window_size=(1024, 768))
    plotter.set_background("white")
    plotter.camera_position = spec["camera"]
    plotter.add_mesh(surface, scalars="WoodColor", rgb=True,
                     smooth_shading=True, show_scalar_bar=False)
    plotter.add_text(f"t = {spec['time'] / 86400.0:.2f} days",
                     position="upper_left", font_size=12, color="black")
    img = plotter.screenshot(None, return_img=True)
    plotter.close()
    return spec["index"], img


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--pvd", default="test_wood.pvd")
    p.add_argument("--input", default="params.input")
    p.add_argument("--gif", default="wood.gif")
    p.add_argument("--fps", type=int, default=10)
    p.add_argument("--frames", type=int, default=100,
                   help="Number of frames uniformly spaced in physical time "
                        "(default: number of snapshots in the PVD)")
    p.add_argument("--ring-width", type=float, default=0.04)
    p.add_argument("--section-spacing", type=float, default=0.02)
    p.add_argument("--warp-scale", type=float, default=1.0)
    p.add_argument("--subdivisions", type=int, default=3)
    p.add_argument("-j", "--processes", type=int,
                   default=max(1, (os.cpu_count() or 1) - 1))
    args = p.parse_args()

    pvd_path = pathlib.Path(args.pvd)
    entries = read_pvd(pvd_path)
    if not entries:
        sys.exit(f"No data in {pvd_path}")

    input_path = pathlib.Path(args.input)
    if not input_path.exists():
        input_path = pvd_path.parent / args.input
    pith = read_pith(input_path)

    bounds = pv.read(str(entries[0][1])).bounds
    cx = 0.5 * (bounds[0] + bounds[1])
    cy = 0.5 * (bounds[2] + bounds[3])
    cam_dist = max(bounds[1] - bounds[0], bounds[3] - bounds[2]) * 1.6
    camera = [(cx, cy, cam_dist), (cx, cy, 0.0), (0.0, 1.0, 0.0)]

    times = np.array([t for t, _ in entries], dtype=float)
    num_frames = args.frames if args.frames is not None else len(entries)
    specs = [{
        "index": i, "time": t,
        "vtu_a": str(entries[i0][1]),
        "vtu_b": str(entries[i1][1]),
        "alpha": alpha,
        "pith": pith.tolist(),
        "warp_scale": args.warp_scale,
        "ring_width": args.ring_width,
        "section_spacing": args.section_spacing,
        "subdivisions": args.subdivisions,
        "camera": camera,
    } for i, (t, i0, i1, alpha) in enumerate(time_brackets(times, num_frames))]

    n = max(1, args.processes)
    print(f"Rendering {len(specs)} frames with {n} workers -> {args.gif}")
    images = [None] * len(specs)
    with mp.get_context("spawn").Pool(processes=n) as pool:
        for idx, img in pool.imap_unordered(render_worker, specs):
            images[idx] = img
            print(f"  frame {idx+1}/{len(specs)}")

    import imageio
    imageio.mimsave(args.gif, images, duration=1.0 / args.fps, loop=0)
    print(f"GIF saved to {args.gif}")


if __name__ == "__main__":
    main()
