#!/usr/bin/env python3
# SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

"""
DFG flow around a cylinder benchmark (Schäfer & Turek, 1996): evaluation and plots

Reads the benchmark indicators of the runs of benchmark.py (runs/2d<case>/<name>/<name>.json):

* Case 2D-1 (--case 1, default): plots the drag coefficient, the lift coefficient and the pressure
  difference over the number of unknowns together with the FeatFlow reference values and the
  bounds of the paper (results/convergence_2d1.png). The runs on the default meshes of the
  regression tests (<scheme>_default) are shown as unfilled markers.
* Case 2D-2 (--case 2): evaluates the last full period of the lift coefficient after --tstart.
* Case 2D-3 (--case 3): evaluates the maxima of the drag and lift coefficients and the pressure
  difference at t = 8 s.

For the cases 2D-2 and 2D-3, results/case_2d_<case>.md compares the quantities with the references,
results/case_2d_<case>.png shows the time series and results/indicators_2d_<case>.png the
benchmark indicators over the number of unknowns, together with the FeatFlow results of all
levels and time step sizes and the bounds of the paper.

The results are compared with the FeatFlow results (reference/featflow.json,
written by reference/fetch_featflow.py) and the bounds of the paper (reference/1996.json).
If the FeatFlow time series have been downloaded with reference/fetch_featflow.py, the time series
of the finest FeatFlow run is plotted as well.

Usage (in the build directory test/freeflow/navierstokes/flowaroundcylinder):
    python3 visualize.py --case 1
    python3 visualize.py --case 2 --tstart 5
"""

import argparse
import glob
import json
import os
import re
import sys

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.legend_handler import HandlerTuple
from matplotlib.patches import Patch
from matplotlib.ticker import NullFormatter

from benchmark import RESULTS_DIR, SCHEMES, TIME_SCHEMES, runsDir

plt.rcParams.update({"mathtext.fontset": "cm", "font.family": "serif", "font.size": 11})

REFERENCE_DIR = os.path.join(os.path.dirname(os.path.realpath(__file__)), "reference")

# minimum time between two minima of the lift coefficient in case 2D-2 (period about 0.33 s)
MIN_PERIOD = 0.2
# cylinder diameter [m] and mean inflow velocity [m/s] of the cases 2D-2 and 2D-3
DIAMETER, MEAN_VELOCITY = 0.1, 1.0

LABELS = {
    "diamond_tpfa": "Diamond / CCTpfa",
    "pq1bubble_box": "PQ1Bubble / Box",
    "pq1bubble_diamond": "PQ1Bubble / Diamond (quadrilaterals)",
    "hybrid_pq1bubble_box": "PQ1Bubble (hybrid) / Box",
}
# the marker shows the grid cell shape (also of the FeatFlow results), the color the scheme
CELL_MARKERS = {"triangles": "^", "quadrilaterals": "s"}
CELL_SHAPES = {
    "diamond_tpfa": "triangles",
    "pq1bubble_box": "triangles",
    "pq1bubble_diamond": "quadrilaterals",
    "hybrid_pq1bubble_box": "triangles",
}
MARKERS = {scheme: CELL_MARKERS[shape] for scheme, shape in CELL_SHAPES.items()}
COLORS = {scheme: f"C{i}" for i, scheme in enumerate(SCHEMES)}
# further runs of the same scheme in the time-dependent cases
RUN_LINESTYLES = ["-", "--", ":", "-."]
RUN_FILLSTYLES = ["full", "none", "left", "right"]
SYMBOLS = {"cDrag": r"$C_D$", "cLift": r"$C_L$", "pDiff": r"$\Delta p$ [Pa]"}

REFERENCE_STYLE = {"color": "k", "linestyle": "-", "linewidth": 1.0}
BOUNDS_STYLE = {"color": "0.6", "alpha": 0.4, "linewidth": 0.0}
FINEST_STEP = r"$\Delta t = 1/1600$ s"
FEATFLOW_LABEL = r"FeatFlow $Q_2/P_1^\mathrm{disc}$ (level 6, " + FINEST_STEP + ")"
FEATFLOW_SCHEMES = {
    "q2p1": r"$Q_2/P_1^\mathrm{disc}$",
    "q1p0": r"$\tilde{Q}_1/Q_0$",
}
# line styles of the FeatFlow time series
FEATFLOW_SERIES_STYLES = {
    "q2p1": {"color": "k", "linestyle": "--"},
    "q1p0": {"color": "k", "linestyle": ":"},
}

# FeatFlow data sets of the schemes per case (see reference/fetch_featflow.py): for case 2D-2,
# the Q1~/Q0 results are only available in the old data set, whose force evaluation is of lower
# order (the drag is about 2 % too low on level 2 and 0.2 % on level 6)
FEATFLOW_DATASETS = {
    2: {"q2p1": "new data", "q1p0": "old data, lower-order force evaluation"},
    3: {"q2p1": "", "q1p0": ""},
}
FEATFLOW_STYLES = {
    "q2p1": {"marker": "s", "linestyle": "-", "linewidth": 1.0, "markersize": 4},
    "q1p0": {"marker": "s", "linestyle": "--", "linewidth": 1.0, "markersize": 4, "mfc": "none"},
}

# benchmark indicators (key, axis label) of the cases 2D-2 and 2D-3 with bounds in the paper
INDICATORS = {
    2: [
        ("CDMax", r"$C_{D,\max}$"),
        ("CLMax", r"$C_{L,\max}$"),
        ("dP", r"$\Delta p(t_0 + 1/(2f))$ [Pa]"),
    ],
    3: [
        ("CDMax", r"$C_{D,\max}$"),
        ("CLMax", r"$C_{L,\max}$"),
        ("dP", r"$\Delta p(t = 8\,\mathrm{s})$ [Pa]"),
    ],
}
PAPER_LABEL = "Schäfer & Turek (1996), bounds"
KEYS = {"cDrag": "CD", "cLift": "CL", "pDiff": "dP"}

# rows of the result tables (key, label) of the cases 2D-2 and 2D-3
PERIODIC_ROWS = [
    (key + statistic, f"{label} ${symbol}$")
    for key, symbol in [("CD", "C_D"), ("CL", "C_L"), ("dP", r"\Delta p")]
    for statistic, label in [
        ("Max", "max"),
        ("Min", "min"),
        ("Mean", "mean"),
        ("Amplitude", "amplitude"),
    ]
] + [
    ("period", r"period $1/f$ [s]"),
    ("St", r"Strouhal number $\mathrm{St} = D f / \bar{U}$"),
    ("dP", r"$\Delta p(t_0 + 1/(2f))$"),
]
TRANSIENT_ROWS = [
    ("CDMax", r"$\max C_D$"),
    ("TCDMax", r"$t(\max C_D)$ [s]"),
    ("CLMax", r"$\max C_L$"),
    ("TCLMax", r"$t(\max C_L)$ [s]"),
    ("dP", r"$\Delta p(t = 8\,\mathrm{s})$"),
]


def readReferences():
    """
    The reference values as {"featflow": ..., "paper": ...}, both indexed by the case
    (reference/featflow.json and reference/1996.json)
    """
    references = {}
    for key, fileName in [("featflow", "featflow.json"), ("paper", "1996.json")]:
        with open(os.path.join(REFERENCE_DIR, fileName), encoding="utf-8") as jsonFile:
            references[key] = json.load(jsonFile)
    return references


def readLog(fileName):
    """
    The number of unknowns and the time series of the indicators of a run without JSON file,
    parsed from the log file. None if the log contains no indicators yet.
    """
    indicators = re.compile(
        r"Benchmark indicators at t = (\S+): cDrag = (\S+), cLift = (\S+), pDiff = (\S+)"
    )
    run = {"time": [], "cDrag": [], "cLift": [], "pDiff": [], "finished": False}
    with open(fileName, encoding="utf-8", errors="replace") as logFile:
        for line in logFile:
            if line.startswith("Total number of dofs"):
                run["numDofs"] = {"total": int(line.split(":")[-1])}
            match = indicators.search(line)
            if match:
                for key, value in zip(["time", "cDrag", "cLift", "pDiff"], match.groups()):
                    run[key].append(float(value))
    return run if run["time"] else None


def readCsv(fileName):
    """The time series of the indicators written at every time step (CSV file), None if empty"""
    data = np.genfromtxt(fileName, delimiter=",", names=True, ndmin=1)
    if data.size == 0:
        return None
    run = {key: data[key].tolist() for key in ["time", "cDrag", "cLift", "pDiff"]}
    run["finished"] = False
    return run


def readRun(folder, name):
    """
    The results of a run: the JSON file written at the end, otherwise (run still going or
    aborted) the CSV file written at every time step, otherwise the log file. Returns the run
    and a note on its source (None for a complete JSON file), or (None, None) if there are no
    results.
    """
    jsonName = os.path.join(folder, f"{name}.json")
    if os.path.exists(jsonName):
        try:
            with open(jsonName, encoding="utf-8") as jsonFile:
                return json.load(jsonFile), None
        except json.JSONDecodeError:
            print(f"{name}: incomplete JSON file, using the CSV or log file")

    logName = os.path.join(folder, f"{name}.log")
    csvName = os.path.join(folder, f"{name}.csv")
    csvRun = readCsv(csvName) if os.path.exists(csvName) else None
    logRun = readLog(logName) if os.path.exists(logName) else None
    run = dict(csvRun or logRun or {})
    if not run:
        return None, None

    numDofs = dict(logRun or {}).get("numDofs")
    if numDofs is not None:
        run.setdefault("numDofs", numDofs)
    source = "CSV" if csvRun else "log"
    return run, f"no JSON file, using the {source} file (up to t = {run['time'][-1]} s)"


def readRuns(case, pattern="*"):
    """
    The results of all runs of a case as {name: result}. Runs without a complete JSON file
    (still going or aborted) are marked with "(running)".
    """
    runs = {}
    for folder in sorted(glob.glob(os.path.join(runsDir(case), pattern))):
        name = os.path.basename(folder)
        run, note = readRun(folder, name)
        if run is None:
            continue
        if note is None:
            runs[name] = run
        else:
            print(f"{name}: {note}")
            runs[f"{name} (running)"] = run
    return runs


def saveFigure(figure, resultsDir, fileName):
    """Save a figure into the results folder"""
    os.makedirs(resultsDir, exist_ok=True)
    path = os.path.join(resultsDir, fileName)
    figure.savefig(path, dpi=150)
    plt.close(figure)
    print(f"Wrote {path}")


###############################################################################
# Case 2D-1: grid convergence
###############################################################################


def plotReferenceValues(axes, references):
    """FeatFlow reference values (lines) and bounds of the paper (shaded) of case 2D-1"""
    for axis, short in zip(axes, KEYS.values()):
        axis.axhspan(*references["paper"]["1"][short], **BOUNDS_STYLE)
        axis.axhline(references["featflow"]["1"][short], **REFERENCE_STYLE)


def plotSchemes(axes):
    """
    Indicators of all schemes over the number of unknowns: lines for the refined coarse meshes,
    unfilled markers (same color) for the default meshes of the regression tests
    """
    for scheme in SCHEMES:
        color = COLORS[scheme]
        runs = sorted(readRuns(1, f"{scheme}_r*").values(), key=lambda run: run["numDofs"]["total"])
        default = readRuns(1, f"{scheme}_default").get(f"{scheme}_default")
        for axis, key in zip(axes, SYMBOLS):
            if runs:
                dofs = [run["numDofs"]["total"] for run in runs]
                values = [run[key] for run in runs]
                axis.semilogx(
                    dofs, values, color=color, marker=MARKERS[scheme], label=LABELS[scheme]
                )
            if default is not None:
                axis.semilogx(
                    default["numDofs"]["total"],
                    default[key],
                    linestyle="none",
                    marker=MARKERS[scheme],
                    markersize=10,
                    markerfacecolor="none",
                    markeredgewidth=1.5,
                    color=color,
                )


def plotConvergence(resultsDir):
    """Plot the indicators of case 2D-1 over the number of unknowns"""
    figure, axes = plt.subplots(1, len(SYMBOLS), figsize=(13, 4.6))
    plotReferenceValues(axes, readReferences())
    plotSchemes(axes)
    for axis, symbol in zip(axes, SYMBOLS.values()):
        axis.set_xlabel(r"number of unknowns $N$")
        axis.set_ylabel(symbol)
        axis.grid(True, which="both", linewidth=0.3)
        axis.xaxis.set_minor_formatter(NullFormatter())
    handles, labels = axes[0].get_legend_handles_labels()
    defaultMarker = tuple(
        plt.Line2D([], [], linestyle="none", marker=m, markersize=10, mfc="none", color="k")
        for m in CELL_MARKERS.values()
    )
    handles += [defaultMarker, plt.Line2D([], [], **REFERENCE_STYLE), Patch(**BOUNDS_STYLE)]
    labels += ["default mesh (regression tests)", "FeatFlow reference value", PAPER_LABEL]
    figure.legend(
        handles,
        labels,
        loc="lower center",
        ncol=3,
        fontsize=9,
        handler_map={tuple: HandlerTuple(ndivide=None)},
    )
    figure.suptitle(r"Case 2D-1 ($\mathrm{Re} = 20$): grid convergence")
    figure.tight_layout(rect=(0.0, 0.12, 1.0, 1.0))
    saveFigure(figure, resultsDir, "convergence_2d1.png")


###############################################################################
# Evaluation of the time-dependent cases
###############################################################################


def refinedMinimum(time, signal, i):
    """Time of the minimum of the parabola through the samples i-1, i, i+1"""
    t0, t1, t2 = time[i - 1 : i + 2]
    y0, y1, y2 = signal[i - 1 : i + 2]
    denominator = (t0 - t1) * (t0 - t2) * (t1 - t2)
    a = (t2 * (y1 - y0) + t1 * (y0 - y2) + t0 * (y2 - y1)) / denominator
    b = (t2 * t2 * (y0 - y1) + t1 * t1 * (y2 - y0) + t0 * t0 * (y1 - y2)) / denominator
    return -b / (2.0 * a) if a > 0.0 else t1


def localMinima(time, signal, minSeparation):
    """
    Times of the local minima of a signal. Of minima closer than minSeparation, only the
    smallest is kept, such that noise near a minimum does not produce spurious periods.
    """
    selected = []
    for i in range(1, len(signal) - 1):
        if not signal[i] < signal[i - 1] or not signal[i] <= signal[i + 1]:
            continue
        if selected and time[i] - time[selected[-1]] < minSeparation:
            if signal[i] < signal[selected[-1]]:
                selected[-1] = i
        else:
            selected.append(i)
    return np.array([refinedMinimum(time, signal, i) for i in selected])


def valueAt(time, signal, t):
    """
    Linear interpolation of a signal at time t, linear extrapolation from the last two samples
    if t lies after the end of the series (e.g. FeatFlow series given at the half steps)
    """
    if t <= time[-1]:
        return np.interp(t, time, signal)
    slope = (signal[-1] - signal[-2]) / (time[-1] - time[-2])
    return signal[-1] + slope * (t - time[-1])


def evaluatePeriodic(series, tStart):
    """
    Evaluate the quantities of case 2D-2 over the last full period [ts, ts + 1/f] of the lift
    coefficient after tStart, where ts is a minimum of the lift coefficient (FeatFlow convention,
    mean = (max + min)/2, amplitude = max - min). The paper defines t0 as the time of the maximum
    lift and evaluates Delta p(t0 + 1/(2f)). If t0 + 1/(2f) lies after the end of the series, the
    periodic value Delta p(t0 - 1/(2f)) is used. Returns None if there is no full period.
    """
    mask = series["time"] >= tStart
    minima = localMinima(series["time"][mask], series["cLift"][mask], MIN_PERIOD)
    if len(minima) < 2:
        return None

    periodStart, periodEnd = minima[-2], minima[-1]
    period = periodEnd - periodStart
    cycle = (series["time"] >= periodStart) & (series["time"] <= periodEnd)
    result = {
        "periodStart": periodStart,
        "period": period,
        "St": DIAMETER / (period * MEAN_VELOCITY),
    }
    for name, key in KEYS.items():
        values = series[name][cycle]
        result[key + "Max"] = values.max()
        result[key + "Min"] = values.min()
        result[key + "Mean"] = 0.5 * (values.max() + values.min())
        result[key + "Amplitude"] = values.max() - values.min()

    result["dP"] = pressureDifferenceAfterMaxLift(series, cycle, period)
    return result


def pressureDifferenceAfterMaxLift(series, cycle, period):
    """
    Delta p(t0 + 1/(2f)) with t0 the (refined) time of the maximum lift in the cycle, or the
    periodic value Delta p(t0 - 1/(2f)) if t0 + 1/(2f) lies after the end of the series
    """
    cycleTime, cycleLift = series["time"][cycle], series["cLift"][cycle]
    iMax = int(np.argmax(cycleLift))
    t0 = cycleTime[iMax]
    if 0 < iMax < len(cycleLift) - 1:
        t0 = refinedMinimum(cycleTime, -cycleLift, iMax)
    target = t0 + 0.5 * period
    if target > series["time"][-1]:
        target = t0 - 0.5 * period
    return valueAt(series["time"], series["pDiff"], target)


def evaluateTransient(series):
    """
    Evaluate the quantities of case 2D-3: maxima of the drag and lift coefficients and their
    times, and the pressure difference at t = 8 s (FeatFlow series end at the last half step,
    their value is extrapolated to t = 8 s). Returns None if the run did not reach t = 8 s.
    """
    if series["time"][-1] < 8.0 - 0.01:
        return None
    iDrag = np.argmax(series["cDrag"])
    iLift = np.argmax(series["cLift"])
    return {
        "CDMax": series["cDrag"][iDrag],
        "TCDMax": series["time"][iDrag],
        "CLMax": series["cLift"][iLift],
        "TCLMax": series["time"][iLift],
        "dP": valueAt(series["time"], series["pDiff"], 8.0),
    }


def timeSeries(run):
    """The time series of the indicators of a run as numpy arrays"""
    return {key: np.asarray(run[key]) for key in ["time", "cDrag", "cLift", "pDiff"]}


###############################################################################
# Tables and plots of the time-dependent cases
###############################################################################


def referenceCell(value, reference):
    """Table cell with a reference: a range with a check mark or a value with the relative error"""
    if reference is None:
        return "-"
    if isinstance(reference, list):
        tolerance = 1e-9 * max(abs(reference[0]), abs(reference[1]))
        inside = "✓" if reference[0] - tolerance <= value <= reference[1] + tolerance else "✗"
        return f"[{reference[0]:.6g}, {reference[1]:.6g}] {inside}"
    error = 100.0 * abs(value - reference) / abs(reference)
    return f"{reference:.6g} ({error:.2f} %)"


def finestFeatflow(references, case):
    """The FeatFlow Q2/P1disc run with the smallest time step size and the most unknowns"""
    return list(references["featflow"][str(case)]["q2p1"].values())[-1][-1]


def writeReport(case, results, rows, references, fileName):
    """Markdown table comparing the results of all runs with the references"""
    columns = [
        (PAPER_LABEL, references["paper"][str(case)]),
        (FEATFLOW_LABEL, finestFeatflow(references, case)),
    ]
    lines = [f"# Case 2D-{case}", ""]
    for name, result in results.items():
        if result is None:
            lines += [f"## {name}", "", "No result (the run is not finished or not periodic).", ""]
            continue
        headers = ["Quantity", "DuMux"] + [header for header, _ in columns]
        lines += [f"## {name}", "", "| " + " | ".join(headers) + " |"]
        lines.append("| :--- |" + " :---: |" * (len(headers) - 1))
        for key, label in rows:
            cells = [label, f"{result[key]:.6g}"]
            cells += [referenceCell(result[key], values.get(key)) for _, values in columns]
            lines.append("| " + " | ".join(cells) + " |")
        lines.append("")
    with open(fileName, "w", encoding="utf-8") as reportFile:
        reportFile.write("\n".join(lines))
    print(f"Wrote {fileName}")


def fetchFeatflow():
    """The module reference/fetch_featflow.py (reading the downloaded FeatFlow data)"""
    sys.path.insert(0, REFERENCE_DIR)
    import fetch_featflow  # pylint: disable=import-outside-toplevel,import-error

    return fetch_featflow


def featflowCurve(scheme, series, dataset):
    """Plotted FeatFlow time series as (label, series, line style), None if not downloaded"""
    if series is None:
        return None
    steps = int(round(1.0 / np.median(np.diff(series["time"][1:]))))
    details = [dataset, "level 6", f"$\\Delta t = 1/{steps}$ s"]
    label = f"FeatFlow {FEATFLOW_SCHEMES[scheme]} ({', '.join(d for d in details if d)})"
    return label, series, FEATFLOW_SERIES_STYLES[scheme]


def featflowCurves(case):
    """The finest Q2/P1disc FeatFlow run (2D-2: new data), if downloaded"""
    dataset = FEATFLOW_DATASETS[case]["q2p1"]
    curve = featflowCurve("q2p1", fetchFeatflow().finestSeries(case, "q2p1"), dataset)
    return [curve] if curve else []


def oldFeatflowCurves():
    """
    Old FeatFlow data of case 2D-2 for both discretizations on the same level and time step
    size (level 6, dt = 1/800 s), if downloaded
    """
    curves = [
        featflowCurve(scheme, fetchFeatflow().oldSeries(scheme), "old data")
        for scheme in FEATFLOW_SCHEMES
    ]
    return [curve for curve in curves if curve]


def decorateAxes(figure, axes):
    """Axis labels and grid of the time series plots, legend below the plots"""
    for axis, key in zip(axes, KEYS):
        axis.set_ylabel(SYMBOLS[key])
        axis.grid(True, linewidth=0.3)
    axes[-1].set_xlabel(r"$t$ [s]")
    handles, labels = axes[0].get_legend_handles_labels()
    figure.legend(handles, labels, loc="lower center", ncol=2, fontsize=9)
    return (len(labels) + 1) // 2


def timeWindow(case, seriesOfRuns):
    """
    Plotted time interval: for case 2D-2 the interval of the new FeatFlow data (t in [9, 10] s)
    if it has been downloaded, else everything (None). Everything is also shown as long as no
    DuMux run has reached this interval yet, such that the progress of running runs is visible.
    """
    curves = featflowCurves(case) if case == 2 else []
    if not curves:
        return None
    time = curves[0][1]["time"]
    if seriesOfRuns and all(series["time"][-1] < time[0] for series in seriesOfRuns.values()):
        return None
    return time[0], time[-1]


def trimmed(series, window):
    """The part of a time series inside the time window (None: everything)"""
    if window is None:
        return series
    inside = (series["time"] >= window[0]) & (series["time"] <= window[1])
    return {key: np.asarray(values)[inside] for key, values in series.items()}


def checkTimePoints(plotName, namedSeries, tolerance=1e-8):
    """
    Report whether the time series of a plot are given at the same points in time as the first
    one (FeatFlow: half steps of the Crank-Nicolson method, DuMux: the time steps)
    """
    namedSeries = [(name, np.asarray(s["time"])) for name, s in namedSeries if len(s["time"]) > 1]
    if len(namedSeries) < 2:
        return
    reference, referenceTimes = namedSeries[0]
    for name, times in namedSeries[1:]:
        nearest = np.clip(np.searchsorted(referenceTimes, times), 1, len(referenceTimes) - 1)
        distance = np.minimum(
            np.abs(times - referenceTimes[nearest - 1]), np.abs(times - referenceTimes[nearest])
        )
        shared = int(np.count_nonzero(distance < tolerance))
        if shared == len(times) == len(referenceTimes):
            print(f"{plotName}: {name} has the same time points as {reference}")
            continue
        print(
            f"{plotName}: {name} shares {shared} of its {len(times)} time points with the "
            f"{len(referenceTimes)} of {reference} (time step {np.median(np.diff(times)):.6g} s"
            f" vs. {np.median(np.diff(referenceTimes)):.6g} s,"
            f" largest offset {distance.max():.3g} s)"
        )


def plotTimeSeries(curves, seriesOfRuns, window, title, fileName):
    """
    Time series of the indicators of the FeatFlow curves (label, series, style) and the DuMux
    runs inside the time window, without markers since the time points of the series differ
    """
    curves = [(label, trimmed(series, window), style) for label, series, style in curves]
    seriesOfRuns = {name: trimmed(series, window) for name, series in seriesOfRuns.items()}
    checkTimePoints(
        os.path.basename(fileName),
        [(label, series) for label, series, _ in curves] + list(seriesOfRuns.items()),
    )

    figure, axes = plt.subplots(len(KEYS), 1, figsize=(8, 9), sharex=True)
    for label, series, style in curves:
        for axis, key in zip(axes, KEYS):
            axis.plot(series["time"], series[key], linewidth=1.2, zorder=3, label=label, **style)
    for (name, series), run in zip(seriesOfRuns.items(), runIndicesPerScheme(seriesOfRuns)):
        for axis, key in zip(axes, KEYS):
            axis.plot(
                series["time"],
                series[key],
                color=COLORS[schemeOfRun(name)],
                linestyle=RUN_LINESTYLES[run % len(RUN_LINESTYLES)],
                label=name,
            )

    legendRows = decorateAxes(figure, axes)
    figure.suptitle(title)
    figure.tight_layout(rect=(0.0, 0.015 + 0.025 * legendRows, 1.0, 1.0))
    saveFigure(figure, *os.path.split(fileName))


def plotAllTimeSeries(case, seriesOfRuns, resultsDir):
    """
    Time series of all runs with the finest FeatFlow run (case_2d_<case>.png), and for case 2D-2
    with the old FeatFlow data of both discretizations (case_2d_2_olddata.png)
    """
    window = timeWindow(case, seriesOfRuns)
    name = r"$\mathrm{Re} = 100$, periodic" if case == 2 else r"$\mathrm{Re} = 0 \to 100 \to 0$"
    fileName = os.path.join(resultsDir, f"case_2d_{case}.png")
    plotTimeSeries(featflowCurves(case), seriesOfRuns, window, f"Case 2D-{case} ({name})", fileName)
    oldCurves = oldFeatflowCurves() if case == 2 else []
    if oldCurves:
        title = f"Case 2D-2 ({name}): old FeatFlow data"
        fileName = os.path.join(resultsDir, "case_2d_2_olddata.png")
        plotTimeSeries(oldCurves, seriesOfRuns, window, title, fileName)


def timeStepShades(featflow):
    """Gray value per FeatFlow time step size (black for the smallest)"""
    steps = sorted({dt for runsOfDt in featflow.values() for dt in runsOfDt}, key=fractionValue)
    return {dt: str(0.75 * i / max(1, len(steps) - 1)) for i, dt in enumerate(steps)}


def fractionValue(dt):
    """Value of a time step size given as 1/<steps per second>"""
    return 1.0 / int(dt.split("/")[1])


def plotFeatflowIndicators(axes, case, references):
    """
    Indicators of all FeatFlow runs over the number of unknowns: one line per discretization
    and time step size (darker for smaller time steps). Returns the legend entries (one per
    discretization with the data set it is taken from, one per time step size).
    """
    featflow = references["featflow"][str(case)]
    shades = timeStepShades(featflow)
    for scheme, runsOfDt in featflow.items():
        for dt, runs in runsOfDt.items():
            for axis, (key, _) in zip(axes, INDICATORS[case]):
                axis.semilogx(
                    [run["numDofs"] for run in runs],
                    [run[key] for run in runs],
                    color=shades[dt],
                    **FEATFLOW_STYLES[scheme],
                )

    handles, labels = [], []
    for scheme in featflow:
        dataset = FEATFLOW_DATASETS[case][scheme]
        handles.append(plt.Line2D([], [], color="k", **FEATFLOW_STYLES[scheme]))
        labels.append(f"FeatFlow {FEATFLOW_SCHEMES[scheme]}" + (f" ({dataset})" if dataset else ""))
    for dt in sorted(shades, key=fractionValue, reverse=True):
        handles.append(plt.Line2D([], [], color=shades[dt], linewidth=3.0))
        labels.append(f"FeatFlow, $\\Delta t = {dt}$ s")
    return handles, labels


def schemeOfRun(name):
    """The scheme of a run, the longest scheme name the run name starts with"""
    return max((scheme for scheme in SCHEMES if name.startswith(scheme)), key=len)


def runIndicesPerScheme(names):
    """Index of each run among the runs of its scheme (0 for the first run of a scheme)"""
    counts = {}
    indices = []
    for name in names:
        scheme = schemeOfRun(name)
        indices.append(counts.get(scheme, 0))
        counts[scheme] = indices[-1] + 1
    return indices


def plotDumuxIndicators(axes, case, runs, results):
    """Indicators of the evaluated DuMux runs at their number of unknowns"""
    evaluated = [name for name, result in results.items() if result is not None]
    for name, run in zip(evaluated, runIndicesPerScheme(evaluated)):
        scheme = schemeOfRun(name)
        for axis, (key, _) in zip(axes, INDICATORS[case]):
            axis.semilogx(
                runs[name]["numDofs"]["total"],
                results[name][key],
                linestyle="none",
                marker=MARKERS[scheme],
                fillstyle=RUN_FILLSTYLES[run % len(RUN_FILLSTYLES)],
                markersize=8,
                color=COLORS[scheme],
                label=f"DuMux {name}",
            )


def plotIndicators(case, runs, results, references, resultsDir):
    """
    Benchmark indicators of all runs over the number of unknowns, together with the FeatFlow
    results of all levels and time step sizes and the bounds of the paper
    """
    figure, axes = plt.subplots(1, len(INDICATORS[case]), figsize=(13, 4.8))
    for axis, (key, label) in zip(axes, INDICATORS[case]):
        axis.axhspan(*references["paper"][str(case)][key], **BOUNDS_STYLE)
        axis.set_xlabel(r"number of unknowns $N$")
        axis.set_ylabel(label)
        axis.grid(True, which="both", linewidth=0.3)
        axis.xaxis.set_minor_formatter(NullFormatter())
    handles, labels = plotFeatflowIndicators(axes, case, references)
    plotDumuxIndicators(axes, case, runs, results)

    dumuxHandles, dumuxLabels = axes[0].get_legend_handles_labels()
    handles += dumuxHandles + [Patch(**BOUNDS_STYLE)]
    labels += dumuxLabels + [PAPER_LABEL]
    figure.legend(handles, labels, loc="lower center", ncol=3, fontsize=9)
    title = r"$\mathrm{Re} = 100$, periodic" if case == 2 else r"$\mathrm{Re} = 0 \to 100 \to 0$"
    figure.suptitle(f"Case 2D-{case} ({title}): benchmark indicators")
    figure.tight_layout(rect=(0.0, 0.03 + 0.045 * ((len(labels) + 2) // 3), 1.0, 1.0))
    saveFigure(figure, resultsDir, f"indicators_2d_{case}.png")


def timeSchemeOfRun(name):
    """The time stepping scheme of a run from the suffix of its name, None for older runs"""
    suffix = name.split(" ")[0].rsplit("_", 1)[-1]
    return next((scheme for scheme in TIME_SCHEMES if scheme.lower() == suffix), None)


def evaluateCase(case, args):
    """
    Evaluate the runs of case 2D-2 or 2D-3 with the selected time stepping schemes, write the
    report and plot the time series
    """
    references = readReferences()
    runs = {
        name: run
        for name, run in readRuns(case).items()
        if "all" in args.time_schemes or timeSchemeOfRun(name) in args.time_schemes
    }
    seriesOfRuns = {name: timeSeries(run) for name, run in runs.items()}
    if not seriesOfRuns:
        print(f"No runs of case 2D-{case} found in {runsDir(case)}/, plotting the references only")

    if case == 2:
        results = {name: evaluatePeriodic(s, args.tstart) for name, s in seriesOfRuns.items()}
    else:
        results = {name: evaluateTransient(s) for name, s in seriesOfRuns.items()}
    rows = PERIODIC_ROWS if case == 2 else TRANSIENT_ROWS

    os.makedirs(args.results_dir, exist_ok=True)
    if results:
        reportFile = os.path.join(args.results_dir, f"case_2d_{case}.md")
        writeReport(case, results, rows, references, reportFile)
    plotAllTimeSeries(case, seriesOfRuns, args.results_dir)
    plotIndicators(case, runs, results, references, args.results_dir)


def main():
    """Evaluate and plot the results of a benchmark case"""
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("--case", type=int, choices=[1, 2, 3], default=1, help="benchmark case")
    parser.add_argument(
        "--tstart",
        type=float,
        default=9.0,
        help="2D-2: the periodic regime is evaluated for t >= tstart (FeatFlow window [9, 10] s)",
    )
    parser.add_argument(
        "--time-schemes",
        nargs="+",
        choices=TIME_SCHEMES + ["all"],
        default=["DIRK2"],
        help="2D-2, 2D-3: evaluate the runs with these time stepping schemes (all: also older runs"
        " without the scheme in the name)",
    )
    parser.add_argument("--results-dir", default=RESULTS_DIR, help="output folder")
    args = parser.parse_args()

    if args.case == 1:
        plotConvergence(args.results_dir)
    else:
        evaluateCase(args.case, args)


if __name__ == "__main__":
    main()
