#!/usr/bin/env python3
"""Plot phi=0 interface snapshots from Stefan-tube VTU output."""

import argparse
import colorsys
import os
import tempfile
import xml.etree.ElementTree as ET
from pathlib import Path


def data_array(parent, name):
    for array in parent.iter("DataArray"):
        if array.attrib.get("Name") == name:
            return array
    raise KeyError(name)


def floats(array):
    return [float(v) for v in array.text.split()]


def ints(array):
    return [int(v) for v in array.text.split()]


def read_pvd(path):
    root = ET.parse(path).getroot()
    datasets = []
    for data_set in root.iter("DataSet"):
        time = float(data_set.attrib["timestep"])
        file_name = data_set.attrib["file"]
        datasets.append((time, (path.parent / file_name).resolve()))
    return sorted(datasets)


def select_snapshots(datasets, count):
    if len(datasets) <= count:
        return datasets
    indices = {
        round(i*(len(datasets) - 1)/(count - 1))
        for i in range(count)
    }
    return [datasets[i] for i in sorted(indices)]


def zero_segments(vtu_path):
    root = ET.parse(vtu_path).getroot()
    points = floats(data_array(root, "Coordinates"))
    phi = floats(data_array(root, "phi"))
    connectivity = ints(data_array(root, "connectivity"))
    offsets = ints(data_array(root, "offsets"))
    cell_types = ints(data_array(root, "types"))

    coords = [
        (points[i], points[i + 1])
        for i in range(0, len(points), 3)
    ]

    segments = []
    start = 0
    for offset, cell_type in zip(offsets, cell_types):
        cell = connectivity[start:offset]
        start = offset
        if cell_type != 5 or len(cell) != 3:
            continue

        crossings = []
        for a, b in ((0, 1), (1, 2), (2, 0)):
            ia, ib = cell[a], cell[b]
            pa, pb = phi[ia], phi[ib]
            if pa == 0.0 and pb == 0.0:
                continue
            if pa == 0.0:
                crossings.append(coords[ia])
            elif pb == 0.0:
                crossings.append(coords[ib])
            elif pa*pb < 0.0:
                t = pa/(pa - pb)
                xa, ya = coords[ia]
                xb, yb = coords[ib]
                crossings.append((xa + t*(xb - xa), ya + t*(yb - ya)))

        unique = []
        for point in crossings:
            if not any((abs(point[0] - p[0]) < 1e-12 and abs(point[1] - p[1]) < 1e-12) for p in unique):
                unique.append(point)
        if len(unique) == 2:
            segments.append(unique)

    return segments


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("pvd")
    parser.add_argument("--out", default="stefantube_interfaces.svg")
    parser.add_argument("--count", type=int, default=7)
    args = parser.parse_args()

    datasets = select_snapshots(read_pvd(Path(args.pvd)), args.count)

    mpl_config = Path(tempfile.gettempdir()) / "dumux-matplotlib-cache"
    mpl_config.mkdir(parents=True, exist_ok=True)
    os.environ.setdefault("MPLCONFIGDIR", str(mpl_config))

    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.collections import LineCollection
    from matplotlib.lines import Line2D

    plt.rcParams["svg.fonttype"] = "none"
    fig, ax = plt.subplots(figsize=(7.0, 3.2), constrained_layout=True)

    n = max(1, len(datasets) - 1)
    handles = []
    all_segments = []
    for i, (time, vtu_path) in enumerate(datasets):
        r = i/n
        saturation = 0.35 + 0.65*r
        value = 0.95 - 0.15*r
        alpha = 0.20 + 0.75*r
        linewidth = 0.8 + 1.8*r
        color = colorsys.hsv_to_rgb(0.60, saturation, value)
        segments = zero_segments(vtu_path)
        ax.add_collection(LineCollection(segments, colors=[color], linewidths=linewidth, alpha=alpha))
        all_segments.extend(segments)
        handles.append(Line2D([0], [0], color=color, alpha=alpha, lw=linewidth, label=f"{time:g}"))

    if all_segments:
        xs = [point[0] for segment in all_segments for point in segment]
        ys = [point[1] for segment in all_segments for point in segment]
        x_pad = max(0.001, 0.18*(max(xs) - min(xs)))
        y_pad = max(0.002, 0.06*(max(ys) - min(ys)))
        ax.set_xlim(min(xs) - x_pad, max(xs) + x_pad)
        ax.set_ylim(min(ys) - y_pad, max(ys) + y_pad)
    else:
        ax.autoscale()

    ax.set_aspect("auto")
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_title(r"$\phi=0$ meniscus snapshots")
    ax.grid(True, color="0.90", linewidth=0.6)
    ax.legend(handles=handles, title="t", loc="center left", bbox_to_anchor=(1.02, 0.5), frameon=False)
    fig.savefig(args.out, dpi=200)


if __name__ == "__main__":
    main()
