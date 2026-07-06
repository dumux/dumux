#!/usr/bin/env python3
# SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later
"""
Plot the DuMux PQ1Bubble Cahn-Hilliard / Navier-Stokes rising-bubble QoIs against the
Hysing et al. (2009) two-dimensional bubble-dynamics benchmark reference data.

It builds a single 3-panel figure  (rise velocity | circularity | center of mass)  for one
Hysing case, overlaying one or more DuMux runs on the two reference groups (TP2D, MooNMD).

--------------------------------------------------------------------------------------------
USAGE
    # from a run LOG (the raw stdout of test_ff_stokes_2p_pq1bubble_simplex):
    python3 plot_hysing_benchmark.py --case 1 --out case1.png  run1.log:my-label

    # from an already-parsed CSV (see --dump-csv), or mixing several runs:
    python3 plot_hysing_benchmark.py --case 2 --out case2.png  a.log:coarse  b.csv:fine

    # just parse a log into a QoI csv (no plot):
    python3 plot_hysing_benchmark.py --dump-csv run1.log run1.csv

Each positional argument is  <path>:<label>  where <path> ends in .log (parsed) or .csv.
Reference data is read from  ../references/  (resolved relative to this script), so the
script is self-contained and can be run from anywhere.

Requires: matplotlib, numpy  (pip install matplotlib numpy).
--------------------------------------------------------------------------------------------
The QoIs printed by the test each time step (see main.cc / problem.hh):
  rise_velocity(_sharp)  : V_c = <u_y> averaged over the bubble; sharp = benchmark def. (phi>0)
  circularity            : c = 2 sqrt(pi A)/P from the phi=0 contour
  centroid_y(_sharp)     : y_c = center of mass;  sharp = benchmark def. (phi>0)
Benchmark reference targets (Case 1): peak V_c=0.2417, min circ=0.9013, y_c(t=3)=1.0813.
"""
import sys, os, re, argparse

ANSI = re.compile(r"\x1b\[[0-9;]*m")
_pat = {
    "time":  re.compile(r"time:\s*([-\d.eE+]+),\s*time step size"),
    "cy":    re.compile(r"centroid_y\s*=\s*([-\d.eE+]+).*centroid_y_sharp\s*=\s*([-\d.eE+]+).*volume\s*=\s*([-\d.eE+]+)"),
    "circ":  re.compile(r"(?<!_)circularity\s*=\s*([-\d.eE+]+).*perimeter\s*=\s*([-\d.eE+]+).*area\s*=\s*([-\d.eE+]+)"),
    "eps":   re.compile(r"eps_eff\s*=\s*([-\d.eE+]+)"),
    "maxv":  re.compile(r"max\|v\|\s*=\s*([-\d.eE+]+)"),
    "rv":    re.compile(r"rise_velocity\s*=\s*([-\d.eE+]+).*rise_velocity_sharp\s*=\s*([-\d.eE+]+)"),
}
COLS = ["time", "area", "circ", "cy", "cy_sharp", "vol", "rv", "rv_sharp", "maxv", "eps"]


def parse_log(fn):
    rows, cur = [], {}
    for raw in open(fn, errors="replace"):
        line = ANSI.sub("", raw)
        m = _pat["time"].search(line)
        if m:
            cur = {"time": float(m.group(1))}; continue
        if not cur:
            continue
        if (m := _pat["cy"].search(line)):
            cur["cy"], cur["cy_sharp"], cur["vol"] = map(float, m.groups()); continue
        if (m := _pat["circ"].search(line)):
            cur["circ"], cur["perim"], cur["area"] = map(float, m.groups()); continue
        if (m := _pat["eps"].search(line)):
            cur["eps"] = float(m.group(1)); continue
        if (m := _pat["maxv"].search(line)):
            cur["maxv"] = float(m.group(1)); continue
        if (m := _pat["rv"].search(line)):
            cur["rv"], cur["rv_sharp"] = float(m.group(1)), float(m.group(2))
            rows.append(cur); cur = {}
    return rows


def load_csv(fn):
    import csv
    return list(csv.DictReader(open(fn)))


def load_run(path):
    """Return dict-of-lists keyed by COLS for a .log or .csv path."""
    rows = parse_log(path) if path.endswith(".log") else load_csv(path)
    out = {}
    for k in set(COLS) | {"circ", "rv", "rv_sharp", "cy", "cy_sharp"}:
        vals = []
        for r in rows:
            v = r.get(k)
            if v not in (None, ""):
                vals.append(float(v))
        if vals:
            out[k] = vals
    return out


def write_csv(rows, fn):
    with open(fn, "w") as f:
        f.write(",".join(COLS) + "\n")
        for r in rows:
            f.write(",".join(f"{r.get(c):.8g}" if isinstance(r.get(c), float) else "" for c in COLS) + "\n")


def load_ref(fn):
    t, circ, comy, rv = [], [], [], []
    for l in open(fn):
        if l.startswith("#") or not l.strip():
            continue
        a = l.split()
        t.append(float(a[0])); circ.append(float(a[1])); comy.append(float(a[2])); rv.append(float(a[3]))
    return t, circ, comy, rv


def clip(t, *ys, tmax=3.0):
    idx = [i for i, tt in enumerate(t) if tt <= tmax + 1e-9]
    return ([t[i] for i in idx], *[[y[i] for i in idx] for y in ys])


def main():
    ap = argparse.ArgumentParser(description="Plot Hysing rising-bubble QoIs vs reference.")
    ap.add_argument("--case", choices=["1", "2"], help="Hysing benchmark case")
    ap.add_argument("--out", default="hysing.png", help="output PNG path")
    ap.add_argument("--tmax", type=float, default=3.0, help="clip plots to this end time")
    ap.add_argument("--refdir", default=None, help="reference-data directory (default ../references)")
    ap.add_argument("--dump-csv", nargs=2, metavar=("LOG", "CSV"), help="parse LOG to CSV and exit")
    ap.add_argument("runs", nargs="*", help="<path.log|path.csv>:<label> entries")
    args = ap.parse_args()

    if args.dump_csv:
        rows = parse_log(args.dump_csv[0]); write_csv(rows, args.dump_csv[1])
        print(f"parsed {len(rows)} rows -> {args.dump_csv[1]}"); return
    if not args.case or not args.runs:
        ap.error("need --case and at least one run (unless --dump-csv)")

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    here = os.path.dirname(os.path.abspath(__file__))
    refdir = args.refdir or os.path.join(here, "..", "references")

    if args.case == "1":
        refs = [("hysing2009_case1_group1_TP2D.txt",  "Ref: TP2D (grp1)",   "#111111", "-"),
                ("hysing2009_case1_group3_MooNMD.txt", "Ref: MooNMD (grp3)", "#555555", "--")]
        target = dict(peakV=0.2417, minC=0.9013, comy=1.0813)
        title = "Hysing Case 1  (rho 1000/100, mu 10/1, sigma=24.5)"
    else:
        refs = [("hysing2009_case2_group1_TP2D.txt",  "Ref: TP2D (grp1, no break-up)",  "#111111", "-"),
                ("hysing2009_case2_group3_MooNMD.txt", "Ref: MooNMD (grp3, break-up)",   "#555555", "--")]
        target = dict(peakV=None, minC=None, comy=None)
        title = "Hysing Case 2  (rho 1000/1, mu 10/0.1, sigma=1.96)"

    fig, ax = plt.subplots(1, 3, figsize=(16.5, 5.0))
    run_colors = ["#d62728", "#1f77b4", "#2ca02c", "#9467bd", "#ff7f0e"]

    for fn, lab, col, ls in refs:
        p = os.path.join(refdir, fn)
        if not os.path.exists(p):
            print(f"WARN missing reference {p}", file=sys.stderr); continue
        t, circ, comy, rv = clip(*load_ref(p), tmax=args.tmax)
        ax[0].plot(t, rv, col, ls=ls, lw=2.0, label=lab, zorder=1)
        ax[1].plot(t, circ, col, ls=ls, lw=2.0, label=lab, zorder=1)
        ax[2].plot(t, comy, col, ls=ls, lw=2.0, label=lab, zorder=1)

    for i, entry in enumerate(args.runs):
        path, _, lab = entry.partition(":")
        lab = lab or os.path.basename(path)
        d = load_run(path); c = run_colors[i % len(run_colors)]
        ax[0].plot(d["time"], d["rv_sharp"], c, lw=1.7, label=f"{lab} (sharp)", zorder=3)
        ax[0].plot(d["time"], d["rv"], c, lw=1.0, ls=":", alpha=0.7, label=f"{lab} (diffuse)", zorder=3)
        ax[1].plot(d["time"], d["circ"], c, lw=1.7, label=lab, zorder=3)
        ax[2].plot(d["time"], d.get("cy_sharp", d["cy"]), c, lw=1.7, label=f"{lab} (sharp)", zorder=3)
        ax[2].plot(d["time"], d["cy"], c, lw=1.0, ls=":", alpha=0.7, label=f"{lab} (diffuse)", zorder=3)
        ip = max(range(len(d["rv_sharp"])), key=lambda k: d["rv_sharp"][k])
        ax[0].plot(d["time"][ip], d["rv_sharp"][ip], "o", color=c, ms=5, zorder=4)
        print(f"{lab}: peak V_c(sharp)={d['rv_sharp'][ip]:.4f}@t={d['time'][ip]:.3f}  "
              f"min circ={min(d['circ']):.4f}  y_c(end)={d.get('cy_sharp',d['cy'])[-1]:.4f}")

    if target["peakV"]:
        ax[0].axhline(target["peakV"], color="0.7", lw=0.8, ls=":")
        ax[0].annotate(f"ref peak {target['peakV']:.4f}", (0.05, target["peakV"]),
                       fontsize=8, va="bottom", color="0.4")
        ax[1].axhline(target["minC"], color="0.7", lw=0.8, ls=":")

    for a, ttl, yl in zip(ax, ["Rise velocity  V_c(t)", "Circularity  c(t)", "Center of mass  y_c(t)"],
                          ["V_c", "circularity", "y_c"]):
        a.set_xlabel("time"); a.set_ylabel(yl); a.set_title(ttl)
        a.set_xlim(0, args.tmax); a.grid(alpha=0.3); a.legend(fontsize=7.5, loc="best")

    fig.suptitle(title + "   —   DuMux PQ1Bubble CH-NS (Aland-Voigt Model 1) vs reference",
                 fontsize=12, y=1.00)
    fig.tight_layout()
    fig.savefig(args.out, dpi=130, bbox_inches="tight")
    print("wrote", args.out)


if __name__ == "__main__":
    main()
