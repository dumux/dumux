#!/usr/bin/env python3
"""Run the coupled epsilon/h Stefan convergence check."""

import argparse
import csv
import subprocess
import sys
from pathlib import Path


CASES = [
    (2.0e-2, 1.25e-3, 5, 0.01),
    (1.0e-2, 6.25e-4, 6, 0.01),
    (5.0e-3, 3.125e-4, 6, 0.01),
    (2.5e-3, 1.5625e-4, 7, 0.005),
    (1.25e-3, 7.8125e-5, 8, 0.0025),
]


def case_name(eps):
    return f"test_ff_stokes_2p_stefan_conv_eps_{eps:.4g}".replace(".", "p").replace("-", "m")


def main():
    here = Path(__file__).resolve().parent
    repo = here.parents[4]

    parser = argparse.ArgumentParser()
    parser.add_argument("--build-dir", type=Path, default=repo / "build-cmake")
    parser.add_argument("--t-end", type=float, default=0.05)
    parser.add_argument("--k", type=float, default=500.0)
    parser.add_argument("--mobility-factor", type=float, default=1e-3)
    args = parser.parse_args()

    exe = args.build_dir / "test/freeflow/navierstokes/2p/stefan/test_ff_stokes_2p_stefan"
    params = here / "params.input"
    if not exe.exists():
        subprocess.run(
            ["cmake", "--build", str(args.build_dir), "--target", "test_ff_stokes_2p_stefan", "-j2"],
            check=True,
        )

    summary = []
    for eps, h_min, max_level, refine_tol in CASES:
        mobility = args.mobility_factor*eps
        name = case_name(eps)
        cmd = [
            str(exe),
            str(params),
            "-Problem.Name",
            name,
            "-Problem.EvaporationRateCoeff",
            str(args.k),
            "-Problem.InterfaceThickness",
            str(eps),
            "-Problem.Mobility",
            str(mobility),
            "-Adaptive.MinElementSize",
            str(h_min),
            "-Adaptive.MaxLevel",
            str(max_level),
            "-Adaptive.InitMaxLevel",
            str(max_level),
            "-Adaptive.RefineTolerance",
            str(refine_tol),
            "-Adaptive.CoarsenTolerance",
            "0.0",
            "-TimeLoop.TEnd",
            str(args.t_end),
            "-TimeLoop.OutputInterval",
            str(args.t_end),
        ]
        print(f"running eps={eps:g}, h_min={h_min:g}")
        completed = subprocess.run(
            cmd,
            cwd=args.build_dir,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
        )
        if completed.returncode:
            print(completed.stdout[-4000:])
            return completed.returncode

        with (args.build_dir / f"{name}_qoi.csv").open() as f:
            last = list(csv.DictReader(f))[-1]
        summary.append({
            "eps": eps,
            "h_min": h_min,
            "mobility": mobility,
            "xi_error": float(last["xi_error"]),
            "xi_time_l2": float(last["xi_time_l2"]),
            "vapor_l2_rel": float(last["vapor_l2_rel"]),
            "evap_rate": float(last["evap_rate"]),
            "evap_rate_stefan": float(last["evap_rate_stefan"]),
        })

    out = args.build_dir / "test_ff_stokes_2p_stefan_convergence.csv"
    with out.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=summary[0].keys())
        writer.writeheader()
        writer.writerows(summary)

    print("eps,h_min,mobility,xi_time_l2,xi_error,vapor_l2_rel")
    for row in summary:
        print(
            f"{row['eps']:.6g},{row['h_min']:.6g},{row['mobility']:.6g},"
            f"{row['xi_time_l2']:.6g},{row['xi_error']:.6g},"
            f"{row['vapor_l2_rel']:.6g}"
        )
    print(f"summary: {out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
