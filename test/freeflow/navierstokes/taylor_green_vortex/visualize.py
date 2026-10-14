#!/usr/bin/env python3
# SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later
"""
Plot the results of the Taylor-Green vortex benchmark (2D)
run by benchmark.py, see README.md:

  * spatial:   errors over the mesh sizes h^v and h^p with the expected slopes
               (taylorgreen_spatial_<dim>d.png; hybrid PQ2 with default and higher-order
               quadrature rules in taylorgreen_spatial_<dim>d_quadrature.png)
  * temporal:  errors over the time step size with the expected slopes
               (taylorgreen_temporal_<dim>d.png)
  * energy:    kinetic energy E and enstrophy eps over time and their relative deviations from
               the analytical values (taylorgreen_energy_<dim>d.png)

The script reads the manifests runs/taylorgreen_<study>_<dim>d.json written by benchmark.py
and has to be executed in the same directory. The figures are written to the directory
results/ (option --output-dir).
"""

import argparse
import csv
import io
import json
import math
import re
import shutil
import sys
from fractions import Fraction
from pathlib import Path

ERROR_KEYS = ("velocityL2", "velocityH1", "pressureL2", "pressureH1")
# notation as in README.md
ERROR_LABELS = {
    "velocityL2": r"$\|\mathbf{u} - \mathbf{u}_h\|_{L^2(\Omega)}$",
    "velocityH1": r"$\|\mathbf{u} - \mathbf{u}_h\|_{H^1(\Omega)}$",
    "pressureL2": r"$\|p - p_h\|_{L^2(\Omega)}$",
    "pressureH1": r"$\|p - p_h\|_{H^1(\Omega)}$",
}
ENERGY_LABEL = r"$E = \frac{1}{2}\int_\Omega \rho \|\mathbf{u}\|^2 \, \mathrm{d}x$"
ENSTROPHY_LABEL = (
    r"$\boldsymbol{\varepsilon} = \frac{1}{2}\int_\Omega \rho"
    r" \|\boldsymbol{\omega}\|^2 \, \mathrm{d}x$"
)
TITLES = {2: "2D Taylor-Green vortex"}
# index of the color in the color cycle of each momentum scheme
SCHEME_COLORS = {
    "pq1bubble": 0,
    "pq1bubblehybrid": 1,
    "pq2hybrid": 2,
    "pq2hybridhoquad": 4,
    "pq3hybrid": 3,
}
SCHEME_STYLES = {"pq2hybridhoquad": {"marker": "s", "markerfacecolor": "none", "linestyle": "-"}}
DEFAULT_STYLE = {"marker": "o", "linestyle": "-"}
# hybrid PQ2 with higher-order quadrature rules is compared to hybrid PQ2 in a separate figure
QUADRATURE_SCHEMES = ("pq2hybrid", "pq2hybridhoquad")
# directory of the manifests and the runs written by benchmark.py
RUN_DIRECTORY = Path("runs")


def latexRenderingWorks(plt):
    """Enable LaTeX text rendering if a sufficiently complete LaTeX installation is found"""
    if not (shutil.which("latex") and shutil.which("dvipng")):
        return False

    plt.rcParams.update(
        {"text.usetex": True, "text.latex.preamble": r"\usepackage{amsmath}\usepackage{bm}"}
    )
    try:
        fig, ax = plt.subplots()
        ax.set_xlabel(ENSTROPHY_LABEL, fontsize=12)
        fig.savefig(io.BytesIO(), format="png")
        plt.close(fig)
        return True
    except (RuntimeError, OSError):
        plt.close("all")
        print("LaTeX rendering failed, using matplotlib's mathtext instead")
        return False


def importPyplot():
    """Import matplotlib with a LaTeX-like style, returns pyplot and whether LaTeX is used"""
    import matplotlib  # pylint: disable=import-outside-toplevel

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt  # pylint: disable=import-outside-toplevel

    plt.rcParams.update(
        {
            "font.family": "serif",
            "font.size": 11,
            "axes.labelsize": 12,
            "legend.fontsize": 9,
            "lines.markersize": 5,
        }
    )
    usetex = latexRenderingWorks(plt)
    if not usetex:
        # LaTeX-like math rendering without (a complete) LaTeX installation
        plt.rcParams.update(
            {
                "text.usetex": False,
                "font.serif": ["Latin Modern Roman", "DejaVu Serif"],
                "mathtext.fontset": "cm",
                "axes.formatter.use_mathtext": True,
            }
        )
    return plt, usetex


class Context:
    """matplotlib with the LaTeX-like style and the output directory of the figures"""

    def __init__(self, outputDir):
        self.plt, self.usetex = importPyplot()
        self.outputDir = Path(outputDir)

    def tex(self, text):
        """
        Math text for the labels: bold Greek symbols (\\boldsymbol, notation as in README.md)
        require LaTeX, matplotlib's mathtext shows the plain symbols instead
        """
        return text if self.usetex else re.sub(r"\\boldsymbol\{(\\?[A-Za-z]+)\}", r"\1", text)

    def save(self, fig, name):
        """Save and close a figure"""
        self.outputDir.mkdir(parents=True, exist_ok=True)
        fileName = self.outputDir / name
        fig.savefig(fileName, dpi=200)
        self.plt.close(fig)
        print(f"Wrote {fileName}")


def loadManifest(study, dim):
    """The manifest of a study written by benchmark.py (or None)"""
    fileName = RUN_DIRECTORY / f"taylorgreen_{study}_{dim}d.json"
    if not fileName.is_file():
        print(f"{fileName} not found, run 'benchmark.py --study {study} --dim {dim}' first")
        return None
    with open(fileName) as manifestFile:
        return json.load(manifestFile)


def readRows(errorFile):
    """The rows of an error file"""
    with open(errorFile, newline="") as errors:
        reader = csv.DictReader(errors)
        return [{key: float(value) for key, value in row.items()} for row in reader]


def formatOrder(order):
    """Format a (rational) convergence order for LaTeX, e.g. 3/2"""
    fraction = Fraction(order).limit_denominator(12)
    return f"{fraction.numerator}" if fraction.denominator == 1 else f"{fraction}"


def addReferenceSlope(ax, data, order, color, anchor=-1):
    """
    Add a dashed line error ~ x^order through the data point (xs, errors) of the given run
    (default: the finest run; the coarsest run, anchor=0, if the finest errors are dominated by
    other errors) and return the end of the line opposite to the anchor (for its label)
    """
    xs, errors = data
    xRange = [min(xs), max(xs)]
    reference = [errors[anchor] * (x / xs[anchor]) ** order for x in xRange]
    ax.loglog(xRange, reference, "--", color=color, linewidth=1.0, alpha=0.8)
    labelEnd = 0 if xs[anchor] == xRange[1] else 1
    return labelEnd, (xRange[labelEnd], reference[labelEnd])


def labelReferenceSlope(ax, order, end, color):
    """Label a reference slope with its order at the given end of the line"""
    labelEnd, position = end
    ax.annotate(
        rf"${formatOrder(order)}$",
        xy=position,
        xytext=(4, 0) if labelEnd == 1 else (-4, 0),
        textcoords="offset points",
        ha="left" if labelEnd == 1 else "right",
        va="center",
        fontsize=10,
        color=color,
    )


def addLegend(fig, handles, labels):
    """Legend below all panels"""
    ncol = min(len(labels), 3)
    try:
        fig.legend(handles, labels, loc="outside lower center", ncol=ncol)
    except ValueError:
        # matplotlib < 3.7: reserve space for the legend below the panels manually
        fig.set_layout_engine("none")
        numRows = math.ceil(len(labels) / ncol)
        legendHeight = (0.2 + 0.22 * numRows) / fig.get_figheight()
        titleHeight = 0.55 / fig.get_figheight()
        fig.tight_layout(rect=(0, legendHeight, 1, 1 - titleHeight))
        fig.legend(handles, labels, loc="lower center", ncol=ncol, frameon=False)


def groupColors(context, groups):
    """
    Colors of the groups: each momentum scheme has the same color in all plots; if a scheme
    appears in several groups (e.g. with different time schemes), the groups get their own colors
    """
    colorCycle = context.plt.rcParams["axes.prop_cycle"].by_key()["color"]
    schemes = [group["scheme"] for group in groups]
    if len(set(schemes)) == len(schemes):
        return [colorCycle[SCHEME_COLORS[scheme]] for scheme in schemes]
    return [colorCycle[i % len(colorCycle)] for i in range(len(groups))]


def groupStyle(group):
    """Marker and line style of a group"""
    return SCHEME_STYLES.get(group["scheme"], DEFAULT_STYLE)


def lastRows(group):
    """The last row (final time) of the error file of every run of a group"""
    return [readRows(run["errors"])[-1] for run in group["runs"]]


def meshSize(row, key):
    """h^v for velocity and h^p for pressure errors (see Schneider & Koch, 2024)"""
    return row["hVelocity"] if key.startswith("velocity") else row["hPressure"]


def meshSizeLabel(key):
    """Axis label of the mesh size of error norm key"""
    return r"$h^v$ [m]" if key.startswith("velocity") else r"$h^p$ [m]"


def timeStepSize(row, _key):
    """The time step size of a run"""
    return row["dt"]


def plotErrorNorm(ax, key, groups, x, anchor):
    """
    Plot an error norm of all groups (with their colors) together with the expected slopes;
    x(row, key) is the plotted abscissa (a step size, such that the expected order is the slope),
    the expected slopes go through the run with index anchor
    """
    labeledOrders = set()  # label each expected order only once per panel
    for group, color in groups:
        rows = lastRows(group)
        data = ([x(row, key) for row in rows], [row[key] for row in rows])
        ax.loglog(*data, color=color, label=group["label"], **groupStyle(group))
        if key in group["expected"]:
            order = group["expected"][key]
            end = addReferenceSlope(ax, data, order, color, anchor)
            if order not in labeledOrders:
                labelReferenceSlope(ax, order, end, color)
            labeledOrders.add(order)
    ax.set_ylabel(ERROR_LABELS[key])
    ax.grid(True, which="both", alpha=0.3)
    ax.margins(x=0.15)


def plotConvergence(context, groups, x, xLabel, anchor=-1):
    """
    Plot all error norms of the groups together with the expected slopes (see plotErrorNorm),
    xLabel(key) is the axis label of the abscissa of error norm key
    """
    fig, axes = context.plt.subplots(2, 2, figsize=(10, 9), constrained_layout=True)
    coloredGroups = list(zip(groups, groupColors(context, groups)))
    for ax, key in zip(axes.flat, ERROR_KEYS):
        plotErrorNorm(ax, key, coloredGroups, x, anchor)
        ax.set_xlabel(xLabel(key))

    # one legend for all panels, with a proxy entry for the reference slopes
    handles, labels = axes.flat[0].get_legend_handles_labels()
    handles.append(context.plt.Line2D([], [], color="gray", linestyle="--", linewidth=1.0))
    labels.append("expected order")
    addLegend(fig, handles, labels)
    return fig, axes


def plotSpatial(context, dim):
    """
    Errors of the stationary problem under grid refinement, hybrid PQ2 with higher-order
    quadrature rules is compared to hybrid PQ2 in a separate figure
    """
    manifest = loadManifest("spatial", dim)
    if manifest is None:
        return False
    groups = manifest["groups"]
    figures = {
        f"taylorgreen_spatial_{dim}d.png": (
            [g for g in groups if g["scheme"] != "pq2hybridhoquad"],
            "spatial convergence",
        ),
        f"taylorgreen_spatial_{dim}d_quadrature.png": (
            [g for g in groups if g["scheme"] in QUADRATURE_SCHEMES],
            "hybrid PQ2 with default and higher-order quadrature rules",
        ),
    }
    for fileName, (figureGroups, title) in figures.items():
        if not figureGroups:
            continue
        fig, _ = plotConvergence(context, figureGroups, meshSize, meshSizeLabel)
        fig.suptitle(f"{TITLES[dim]}, stationary problem, {title}")
        context.save(fig, fileName)
    return True


def plotTemporal(context, dim):
    """Errors of the instationary problem under refinement of the time step size"""
    manifest = loadManifest("temporal", dim)
    if manifest is None:
        return False
    parameters = manifest["parameters"]
    # on fine time steps, the spatial error dominates: draw the slopes from the coarsest run
    fig, _ = plotConvergence(
        context, manifest["groups"], timeStepSize, lambda key: r"time step size $\Delta t$ [s]", 0
    )
    fig.suptitle(
        f"{TITLES[dim]}, temporal convergence at $t = {parameters['tEnd']:g}$ s "
        f"(${parameters['cells']}^{dim}$ cells)"
    )
    context.save(fig, f"taylorgreen_temporal_{dim}d.png")
    return True


def plotIndicator(context, axes, groups, indicator):
    """
    Plot an indicator (key, name, symbol, definition) over time, normalized with its analytical
    value at t = 0 (axes[0]), and its relative deviation from the analytical value (axes[1])
    """
    key, name, symbol, definition = indicator
    exact = readRows(groups[0]["runs"][0]["errors"])
    initial = exact[0][f"{key}Exact"]
    axes[0].plot(
        [r["t"] for r in exact], [r[f"{key}Exact"] / initial for r in exact], "k--", label="exact"
    )
    for group, color in zip(groups, groupColors(context, groups)):
        rows = readRows(group["runs"][0]["errors"])
        t = [r["t"] for r in rows]
        axes[0].plot(t, [r[key] / initial for r in rows], color=color, label=group["label"])
        axes[1].plot(t, [100 * (r[key] / r[f"{key}Exact"] - 1) for r in rows], color=color)
    axes[0].set_title(context.tex(f"{name} {definition}"), fontsize=11)
    axes[0].set_ylabel(context.tex(rf"${symbol}_h(t) / {symbol}(0)$"))
    axes[1].set_title(f"relative deviation of the {name}", fontsize=11)
    axes[1].set_ylabel(context.tex(rf"$({symbol}_h - {symbol}) / {symbol}$ [%]"))


def plotEnergy(context, dim):
    """Kinetic energy (left) and enstrophy (right) over time and their relative deviations"""
    manifest = loadManifest("energy", dim)
    if manifest is None:
        return False
    # hybrid PQ2 with higher-order quadrature coincides with hybrid PQ2 (see README.md)
    groups = [g for g in manifest["groups"] if g["scheme"] != "pq2hybridhoquad"]

    fig, axes = context.plt.subplots(2, 2, figsize=(10, 8.5), constrained_layout=True)
    indicators = (
        ("kineticEnergy", "kinetic energy", "E", ENERGY_LABEL),
        ("enstrophy", "enstrophy", r"\boldsymbol{\varepsilon}", ENSTROPHY_LABEL),
    )
    for column, indicator in enumerate(indicators):
        plotIndicator(context, axes[:, column], groups, indicator)
    for ax in axes.flat:
        ax.set_xlabel(r"$t$ [s]")
        ax.grid(True, alpha=0.3)
    addLegend(fig, *axes[0][0].get_legend_handles_labels())

    dt = groups[0]["runs"][0]["dt"]
    fig.suptitle(
        f"{TITLES[dim]}, instationary problem, {manifest['parameters']['timeScheme']}, "
        rf"$\Delta t = {dt:g}$ s, $\nu = {manifest['viscosity']:g}$ m$^2$/s"
    )
    context.save(fig, f"taylorgreen_energy_{dim}d.png")
    return True


def main():
    """Create the plots of the selected studies"""
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("--dim", type=int, nargs="+", choices=(2,), default=[2])
    parser.add_argument(
        "--study",
        nargs="+",
        choices=("spatial", "temporal", "energy", "all"),
        default=["all"],
        help="plots to create (only for the studies run by benchmark.py)",
    )
    parser.add_argument(
        "--output-dir", default="results", help="directory of the figures (default: results)"
    )
    args = parser.parse_args()

    context = Context(args.output_dir)
    plots = {
        "spatial": lambda dim: plotSpatial(context, dim),
        "temporal": lambda dim: plotTemporal(context, dim),
        "energy": lambda dim: plotEnergy(context, dim),
    }
    explicit = "all" not in args.study
    studies = args.study if explicit else plots.keys()
    success = True
    for dim in args.dim:
        for study in studies:
            # missing studies are only an error if the plot was requested explicitly
            success &= plots[study](dim) or not explicit
    sys.exit(0 if success else 1)


if __name__ == "__main__":
    main()
