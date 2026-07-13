#!/usr/bin/env python3
"""Plot Stefan-tube gas-length and evaporation-rate QoIs."""

import argparse
import csv
import os
import tempfile
from pathlib import Path


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("csv")
    parser.add_argument("--out", default="stefantube_qoi.svg")
    args = parser.parse_args()

    with open(args.csv) as f:
        rows = list(csv.DictReader(f))

    t = [float(r["time"]) for r in rows]
    ell = [float(r["ell"]) for r in rows]
    ell_ref = [float(r["ell_tube"]) for r in rows]
    evap = [float(r["evap_rate"]) for r in rows]
    evap_ref = [float(r["evap_rate_tube"]) for r in rows]

    mpl_config = Path(tempfile.gettempdir()) / "dumux-matplotlib-cache"
    mpl_config.mkdir(parents=True, exist_ok=True)
    os.environ.setdefault("MPLCONFIGDIR", str(mpl_config))

    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    plt.rcParams["svg.fonttype"] = "none"
    fig, (ax_ell, ax_rate) = plt.subplots(
        2, 1, figsize=(6.4, 5.0), sharex=True, constrained_layout=True
    )

    ax_ell.plot(t, ell, label=r"phase-field $\ell$", color="#1f77b4", linewidth=2.0)
    ax_ell.plot(t, ell_ref, label=r"tube $\ell$", color="#d62728", linewidth=2.0, linestyle="--")
    ax_ell.set_ylabel("gas length")
    ax_ell.grid(True, color="0.88", linewidth=0.8)
    ax_ell.legend(frameon=False)

    ax_rate.plot(t, evap, label="phase-field rate", color="#1f77b4", linewidth=2.0)
    ax_rate.plot(t, evap_ref, label="tube rate", color="#d62728", linewidth=2.0, linestyle="--")
    ax_rate.set_xlabel("time")
    ax_rate.set_ylabel("evaporation rate")
    ax_rate.grid(True, color="0.88", linewidth=0.8)
    ax_rate.legend(frameon=False)
    fig.savefig(args.out, dpi=200)


if __name__ == "__main__":
    main()
