#!/usr/bin/env python3
# SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later
"""Build, run, and visualize the Buckley-Leverett two-phase test.

The script expects PyVista and Matplotlib for post-processing:

    python3 -m pip install pyvista matplotlib

It builds the DuMux target, runs the executable from the CMake binary
directory, samples the final VTU result along y = 37.5 m, and writes
line-plot and last-frame images.
"""

from __future__ import annotations

import argparse
import glob
import subprocess
import sys
from pathlib import Path


TARGET = "test_2p_buckleyleverett_tpfa"
OUTPUT_NAME = "buckleyleverett"
CASE_REL_DIR = Path("test/porousmediumflow/2p/buckleyleverett")
DEFAULT_POINT_1 = (0.0, 37.5, 0.0)
DEFAULT_POINT_2 = (100.0, 37.5, 0.0)


def source_root() -> Path:
    return Path(__file__).resolve().parents[4]


def default_build_dir(root: Path) -> Path:
    return root / "build-cmake"


def run(command: list[str], cwd: Path | None = None) -> None:
    print(f"+ {' '.join(command)}")
    subprocess.run(command, cwd=cwd, check=True)

def build_target(build_dir: Path) -> None:
    run(["make", TARGET], cwd=build_dir)


def output_dir(build_dir: Path) -> Path:
    return build_dir / CASE_REL_DIR


def remove_old_outputs(case_build_dir: Path, problem_name: str) -> None:
    for pattern in (f"{problem_name}-*.vtu", f"{problem_name}.pvd", f"{problem_name}-*.pvtu"):
        for file_name in glob.glob(str(case_build_dir / pattern)):
            Path(file_name).unlink()


def run_case(case_build_dir: Path, problem_name: str) -> None:
    executable = case_build_dir / TARGET
    if not executable.exists():
        raise FileNotFoundError(f"Could not find executable: {executable}")

    run([str(executable), "params.input", "-Problem.Name", problem_name], cwd=case_build_dir)


def last_vtu(case_build_dir: Path, problem_name: str) -> Path:
    files = sorted(case_build_dir.glob(f"{problem_name}-*.vtu"))
    if not files:
        raise FileNotFoundError(f"No VTU output found for '{problem_name}' in {case_build_dir}")
    return files[-1]


def require_plot_modules():
    try:
        import matplotlib.pyplot as plt
        import pyvista as pv
    except ImportError as error:
        raise SystemExit(
            "Post-processing requires PyVista and Matplotlib. "
            "Install them with: python3 -m pip install pyvista matplotlib"
        ) from error

    return pv, plt


def make_line_plot(vtu_file: Path, image_file: Path, point_1, point_2, samples: int) -> None:
    pv, plt = require_plot_modules()

    mesh = pv.read(vtu_file)
    mesh = mesh.cell_data_to_point_data(pass_cell_data=True)

    line = mesh.sample_over_line(point_1, point_2, resolution=samples - 1)
    x = line.points[:, 0]

    missing = [name for name in ("S_aq", "Sw_exact") if name not in line.array_names]
    if missing:
        raise KeyError(f"Missing sampled field(s): {', '.join(missing)}")

    fig, ax = plt.subplots(figsize=(9, 4.8), constrained_layout=True)
    ax.plot(x, line["S_aq"], label="Sw", linewidth=2)
    ax.plot(x, line["Sw_exact"], label="Sw_exact", linewidth=2)
    ax.set_xlabel("x [m]")
    ax.set_ylabel("saturation [-]")
    ax.set_xlim(min(point_1[0], point_2[0]), max(point_1[0], point_2[0]))
    ax.grid(True, alpha=0.3)
    ax.legend()
    fig.savefig(image_file, dpi=200)
    plt.close(fig)


def make_frame_image(vtu_file: Path, image_file: Path, show: bool) -> None:
    pv, _ = require_plot_modules()

    mesh = pv.read(vtu_file)
    plotter = pv.Plotter(off_screen=not show, window_size=(650, 400))
    plotter.add_mesh(
        mesh,
        scalars="S_aq",
        show_edges=False,
        cmap="coolwarm",
        scalar_bar_args={
            "title": "S_w",
            "vertical": True,
            "position_x": 0.85,
            "position_y": 0.20,
            "width": 0.08,
            "height": 0.60},
    )
    plotter.view_xy()
    plotter.camera.zoom(1.2)


    if show:
        plotter.show(screenshot=str(image_file))
    else:
        plotter.show(screenshot=str(image_file), auto_close=True)


def parse_point(values: list[float]) -> tuple[float, float, float]:
    if len(values) != 3:
        raise argparse.ArgumentTypeError("Point requires exactly three coordinates")
    return tuple(values)


def main() -> int:
    root = source_root()

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--build-dir", type=Path, default=default_build_dir(root),
                        help="CMake build directory")
    parser.add_argument("--problem-name", default=TARGET,
                        help="Problem.Name used for the simulation output")
    parser.add_argument("--lineplot", type=Path, default=None,
                        help="Output path for the line plot PNG")
    parser.add_argument("--frame", type=Path, default=None,
                        help="Output path for the final-frame PNG")
    parser.add_argument("--samples", type=int, default=501,
                        help="Number of samples on the line")
    parser.add_argument("--point-1", nargs=3, type=float, default=DEFAULT_POINT_1,
                        metavar=("X", "Y", "Z"), help="First line-sampling point")
    parser.add_argument("--point-2", nargs=3, type=float, default=DEFAULT_POINT_2,
                        metavar=("X", "Y", "Z"), help="Second line-sampling point")
    parser.add_argument("--show", action="store_true",
                        help="Open an interactive PyVista window for the final frame")
    parser.add_argument("--skip-run", action="store_true",
                        help="Skip build and execution and post-process the latest existing VTU")
    args = parser.parse_args()

    build_dir = args.build_dir.resolve()
    case_build_dir = output_dir(build_dir)
    problem_name = args.problem_name

    if args.samples < 2:
        raise SystemExit("--samples must be at least 2")

    lineplot_file = args.lineplot or case_build_dir / f"{OUTPUT_NAME}_lineplot.png"
    frame_file = args.frame or case_build_dir / f"{OUTPUT_NAME}_sw.png"

    if not args.skip_run:
        build_target(build_dir)
        remove_old_outputs(case_build_dir, problem_name)
        run_case(case_build_dir, problem_name)

    final_vtu = last_vtu(case_build_dir, problem_name)
    print(f"Using final VTU file: {final_vtu}")

    make_line_plot(final_vtu, lineplot_file, tuple(args.point_1), tuple(args.point_2), args.samples)
    print(f"Wrote line plot: {lineplot_file}")

    make_frame_image(final_vtu, frame_file, args.show)
    print(f"Wrote final-frame image: {frame_file}")

    return 0


if __name__ == "__main__":
    sys.exit(main())
