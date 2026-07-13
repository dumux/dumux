#!/usr/bin/env python3
"""Plot Stefan interface-position QoI."""

import argparse
import csv
import os
import tempfile
from pathlib import Path


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("csv")
    parser.add_argument("--out", default="stefan_xi.svg")
    args = parser.parse_args()

    with open(args.csv) as f:
        rows = list(csv.DictReader(f))

    t = [float(r["time"]) for r in rows]
    xi = [float(r["xi"]) for r in rows]
    xi_ref = [float(r["xi_stefan"]) for r in rows]

    mpl_config = Path(tempfile.gettempdir()) / "dumux-matplotlib-cache"
    mpl_config.mkdir(parents=True, exist_ok=True)
    os.environ.setdefault("MPLCONFIGDIR", str(mpl_config))

    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    plt.rcParams["svg.fonttype"] = "none"
    fig, ax = plt.subplots(figsize=(6.4, 3.6), constrained_layout=True)
    ax.plot(t, xi, label=r"phase-field $\xi$", color="#1f77b4", linewidth=2.0)
    ax.plot(t, xi_ref, label=r"Stefan $\xi$", color="#d62728", linewidth=2.0, linestyle="--")
    ax.set_xlabel("time")
    ax.set_ylabel("interface position")
    ax.grid(True, color="0.88", linewidth=0.8)
    ax.legend(frameon=False)
    fig.savefig(args.out, dpi=200)


if __name__ == "__main__":
    main()
