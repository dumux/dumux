#!/usr/bin/env python3
"""Run the Stefan-tube evaporation benchmark and plot the QoIs."""

import argparse
import csv
import subprocess
import sys
from pathlib import Path


def main():
    here = Path(__file__).resolve().parent
    repo = here.parents[4]

    parser = argparse.ArgumentParser()
    parser.add_argument("--build-dir", type=Path, default=repo / "build-cmake")
    parser.add_argument("--t-end", type=float)
    parser.add_argument("overrides", nargs=argparse.REMAINDER)
    args = parser.parse_args()

    exe = args.build_dir / "test/freeflow/navierstokes/2p/stefantube/test_ff_stokes_2p_stefantube"
    params = here / "params.input"
    if not exe.exists():
        subprocess.run(
            ["cmake", "--build", str(args.build_dir), "--target", "test_ff_stokes_2p_stefantube", "-j2"],
            check=True,
        )

    cmd = [str(exe), str(params)]
    if args.t_end is not None:
        cmd += ["-TimeLoop.TEnd", str(args.t_end)]
    if args.overrides:
        extra = args.overrides[1:] if args.overrides[0] == "--" else args.overrides
        cmd += extra

    subprocess.run(cmd, cwd=args.build_dir, check=True)

    csv_path = args.build_dir / "test_ff_stokes_2p_stefantube_qoi.csv"
    plot_path = args.build_dir / "test_ff_stokes_2p_stefantube_qoi.svg"
    subprocess.run([sys.executable, str(here / "plot_qoi.py"), str(csv_path), "--out", str(plot_path)], check=True)

    with csv_path.open() as f:
        rows = list(csv.DictReader(f))
    first = next((r for r in rows if float(r["time"]) > 0.0), rows[0])
    print(
        "first ell QoI: "
        f"t={float(first['time']):.6g}, "
        f"ell={float(first['ell']):.8g}, "
        f"ell_tube={float(first['ell_tube']):.8g}, "
        f"err={float(first['ell_error']):.3g}"
    )
    print(f"plot: {plot_path}")


if __name__ == "__main__":
    main()
