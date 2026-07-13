#!/usr/bin/env python3
"""Run the coupled epsilon/h Stefan-tube convergence check."""

import argparse
import csv
import subprocess
import sys
from pathlib import Path


CASES = [
    # eps, h_min = eps/4, refinement tolerance
    (2.5e-3, 6.25e-4, 0.01),
    (1.25e-3, 3.125e-4, 0.005),
    (6.25e-4, 1.5625e-4, 0.0025),
    (3.125e-4, 7.8125e-5, 0.00125),
]

MAX_LEVEL_CAP = 15


def case_name(eps):
    return f"test_ff_stokes_2p_stefantube_conv_eps_{eps:.4g}".replace(".", "p").replace("-", "m")


def main():
    here = Path(__file__).resolve().parent
    repo = here.parents[4]

    parser = argparse.ArgumentParser()
    parser.add_argument("--build-dir", type=Path, default=repo / "build-cmake")
    parser.add_argument("--t-end", type=float, default=0.10)
    parser.add_argument("--k", type=float, default=5000.0)
    parser.add_argument("--mobility-factor", type=float, default=1e-3)
    parser.add_argument("--rho-liquid", type=float, default=100.0)
    args = parser.parse_args()

    exe = args.build_dir / "test/freeflow/navierstokes/2p/stefantube/test_ff_stokes_2p_stefantube"
    params = here / "params.input"
    if not exe.exists():
        subprocess.run(
            ["cmake", "--build", str(args.build_dir), "--target", "test_ff_stokes_2p_stefantube", "-j2"],
            check=True,
        )

    summary = []
    for eps, h_min, refine_tol in CASES:
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
            str(MAX_LEVEL_CAP),
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
            rows = list(csv.DictReader(f))
        first = rows[0]
        last = rows[-1]
        elapsed = float(last["time"]) - float(first["time"])
        avg_rate = args.rho_liquid*(float(last["ell"]) - float(first["ell"]))/elapsed
        avg_rate_ref = args.rho_liquid*(float(last["ell_tube"]) - float(first["ell_tube"]))/elapsed
        summary.append({
            "eps": eps,
            "h_min": h_min,
            "mobility": mobility,
            "ell_error": float(last["ell_error"]),
            "ell_time_l2": float(last["ell_time_l2"]),
            "vapor_l2_rel": float(last["vapor_l2_rel"]),
            "evap_rate": float(last["evap_rate"]),
            "evap_rate_tube": float(last["evap_rate_tube"]),
            "evap_rate_rel_error": abs(float(last["evap_rate"]) - float(last["evap_rate_tube"]))
                                   /abs(float(last["evap_rate_tube"])),
            "evap_rate_ratio": float(last["evap_rate"])/float(last["evap_rate_tube"]),
            "avg_rate": avg_rate,
            "avg_rate_tube": avg_rate_ref,
            "avg_rate_rel_error": abs(avg_rate - avg_rate_ref)/abs(avg_rate_ref),
            "avg_rate_ratio": avg_rate/avg_rate_ref,
        })

    out = args.build_dir / "test_ff_stokes_2p_stefantube_convergence.csv"
    with out.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=summary[0].keys())
        writer.writeheader()
        writer.writerows(summary)

    print("eps,h_min,mobility,ell_time_l2,ell_error,evap_rate_rel_error,evap_rate_ratio,avg_rate_rel_error,avg_rate_ratio,vapor_l2_rel")
    for row in summary:
        print(
            f"{row['eps']:.6g},{row['h_min']:.6g},{row['mobility']:.6g},"
            f"{row['ell_time_l2']:.6g},{row['ell_error']:.6g},"
            f"{row['evap_rate_rel_error']:.6g},{row['evap_rate_ratio']:.6g},"
            f"{row['avg_rate_rel_error']:.6g},{row['avg_rate_ratio']:.6g},"
            f"{row['vapor_l2_rel']:.6g}"
        )
    print(f"summary: {out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
