#!/usr/bin/env python3
# SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later


"""
# Script to reproduce Figure 5 of
# Koch et al (2020) JCP
# https://doi.org/10.1016/j.jcp.2020.109370
#
# The script also prints a table with errors and
# convergence rates and fails if the mean rates differ
# from the rates of the paper.
#
# Usage: python3 convergence.py [output directory of the figure]
"""

import numpy as np
import subprocess
import sys
import os

def make_program_call(executable, args):
    call = ['./' + executable]
    for k, v in args.items():
        call += ['-' + k, v]
    return call

# remove the old log file and run
def run_simulation(params, filter=[]):
    print("Running ", params["exec"], params["Vessel.Grid.Cells"], "...")
    if params["exec"] not in filter:
        with open(os.devnull, 'w') as devnull:
            subprocess.run(['rm', '-f', params["Problem.OutputFilename"]], stdout=devnull)
            subprocess.run(make_program_call(params["exec"], params))

executables = ["test_md_embedded_1d3d_1p1p_tpfatpfa_average", "test_md_embedded_1d3d_1p1p_tpfatpfa_surface", "test_md_embedded_1d3d_1p1p_tpfatpfa_kernel"]
method = {"test_md_embedded_1d3d_1p1p_tpfatpfa_average":1, "test_md_embedded_1d3d_1p1p_tpfatpfa_surface":3, "test_md_embedded_1d3d_1p1p_tpfatpfa_kernel":5}
label = {"test_md_embedded_1d3d_1p1p_tpfatpfa_average":"line source (LS)", "test_md_embedded_1d3d_1p1p_tpfatpfa_surface":"surface source (CSS)", "test_md_embedded_1d3d_1p1p_tpfatpfa_kernel":"kernel (DS)"}
cells = [10, 20, 40, 80]
radius = 0.1

runs = []
res = {}
for e in executables:
    for c in cells:
        runs.append({"Vessel.Grid.Cells":str(c//2),
                     "Tissue.Grid.Cells":str(c)+" "+str(c)+" "+str(c//2),
                     "SpatialParams.Radius":str(radius),
                     "Problem.OutputFilename":e+str(c)+".log",
                     "Vtk.EnableVtkOutput":"false",
                     "exec":e})
        res[e] = {}

for run in runs:
    run_simulation(run)
    name, numcells, result = (run["exec"], run["Vessel.Grid.Cells"], np.genfromtxt(run["Problem.OutputFilename"]))
    print("Read file: ", run["Problem.OutputFilename"])
    res[name][int(numcells)] = result

# produce output suitable for latex table
# Table 1
table1 = [["", "", "", "", "", "", ""] for i in range(len(cells)+1)]
table2 = [["", "", "", "", "", "", ""] for i in range(len(cells)+1)]
table3 = [["", "", "", "", "", "", ""] for i in range(len(cells)+1)]

def get_errors(result):
    p3d, p1d, q, h = [], [], [], []
    for cells, norms in sorted(result.items()):
        h.append(norms[1])
        p3d.append(norms[2])
        p1d.append(norms[3])
        q.append(norms[4])
    return (p3d, p1d, q, h)

for exec, result in res.items():
    p3d, p1d, q, h = get_errors(result=result)

    rates3d = [(np.log10(p3d[i+1])-np.log10(p3d[i]))/(np.log10(h[i+1])-np.log10(h[i])) for i in range(len(p3d)-1)]
    rates1d = [(np.log10(p1d[i+1])-np.log10(p1d[i]))/(np.log10(h[i+1])-np.log10(h[i])) for i in range(len(p1d)-1)]
    ratesq = [(np.log10(q[i+1])-np.log10(q[i]))/(np.log10(h[i+1])-np.log10(h[i])) for i in range(len(q)-1)]

    ofs = method[exec]
    for i in range(len(h)):
        table1[i][0] = "{:.4f}".format(h[i])
        table2[i][0] = "{:.4f}".format(h[i])
        table3[i][0] = "{:.4f}".format(h[i])
        table1[i][ofs] = "{:.4e}".format(p3d[i])
        table2[i][ofs] = "{:.4e}".format(p1d[i])
        table3[i][ofs] = "{:.4e}".format(q[i])
        if i > 0:
            table1[i][ofs+1] = "{:.4e}".format(rates3d[i-1])
            table2[i][ofs+1] = "{:.4e}".format(rates1d[i-1])
            table3[i][ofs+1] = "{:.4e}".format(ratesq[i-1])
    table1[len(h)][ofs] = "mean rate"
    table2[len(h)][ofs] = "mean rate"
    table3[len(h)][ofs] = "mean rate"
    table1[len(h)][ofs+1] = "{:.4e}".format(np.mean(rates3d[1:]))
    table2[len(h)][ofs+1] = "{:.4e}".format(np.mean(rates1d[1:]))
    table3[len(h)][ofs+1] = "{:.4e}".format(np.mean(ratesq[1:]))

def print_table(table):
    for row in table:
        print(" & ".join(row) + "\\\\")

print("p3d")
print_table(table1)
print("p1d")
print_table(table2)
print("q")
print_table(table3)

#############################################
# verify mean rates against reference
#############################################
# the reference rates from the paper
reference = {"test_md_embedded_1d3d_1p1p_tpfatpfa_average":[1.3483, 1.3545, 1.8449],
             "test_md_embedded_1d3d_1p1p_tpfatpfa_surface":[1.4693, 1.7800, 2.2715],
             "test_md_embedded_1d3d_1p1p_tpfatpfa_kernel":[2.0583, 2.0927, 2.5822]}

for exec, result in res.items():
    samples = len(result)
    ofs = method[exec]
    rate3d = float(table1[len(h)][ofs+1])
    rate1d = float(table2[len(h)][ofs+1])
    rateq = float(table3[len(h)][ofs+1])
    error_to_reference = np.linalg.norm(np.array(reference[exec])-np.array([rate3d, rate1d, rateq]), ord=2)
    if error_to_reference > 0.1:
        print("\nWrong convergence rates for case {}".format(exec))
        print("Expected {}, obtained {}".format(reference[exec], [rate3d, rate1d, rateq]))
        sys.exit(1)

#############################################
# create plot from paper (Figure 5)
#############################################
if len(sys.argv) > 1:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({"font.size": 13, "axes.titlesize": 13, "axes.labelsize": 13,
                         "xtick.labelsize": 12, "ytick.labelsize": 12, "legend.fontsize": 12})

    # color, marker and line style differ per method, so the figure also reads in black and white
    style = {"test_md_embedded_1d3d_1p1p_tpfatpfa_average": dict(color="#2a78d6", marker="o", linestyle="-"),
             "test_md_embedded_1d3d_1p1p_tpfatpfa_surface": dict(color="#eb6834", marker="^", linestyle="--"),
             "test_md_embedded_1d3d_1p1p_tpfatpfa_kernel": dict(color="#1baf7a", marker="s", linestyle="-.")}
    titles = [r"tissue pressure $p_t$", r"vessel pressure $p_v$", r"exchange rate $q$"]
    orders = [(1.5, 2.0), (1.5, 2.0), (2.0, 2.5)]

    fig, axes = plt.subplots(1, 3, figsize=(12, 4.6), layout="constrained")
    for exec, result in res.items():
        p3d, p1d, q, h = get_errors(result=result)
        hR = np.array(h)/radius
        for ax, errors in zip(axes, [p3d, p1d, q]):
            ax.loglog(hR, errors, markersize=7, linewidth=1.8, label=label[exec], **style[exec])

    def slope_triangle(ax, x_fine, y_low, order, width=1.5):
        x_coarse = width*x_fine
        y_high = y_low*width**order
        ax.plot([x_coarse, x_fine, x_coarse, x_coarse], [y_high, y_low, y_low, y_high], color="0.2", linewidth=1.2)
        ax.annotate(f"{order:g}", (x_coarse, np.sqrt(y_low*y_high)), textcoords="offset points", xytext=(-5, 0),
                    ha="right", va="center", color="0.2")
        ax.annotate("1", (np.sqrt(x_coarse*x_fine), y_low), textcoords="offset points", xytext=(0, -3),
                    ha="center", va="top", color="0.2")

    for k, ax in enumerate(axes):
        data = [line.get_ydata() for line in ax.get_lines()]
        lows = []
        for level, order in zip([-1, -2], sorted(orders[k], reverse=True)):
            lows.append(0.3*min(y[level] for y in data))
            slope_triangle(ax, hR[level], lows[-1], order)
        ax.set_ylim(bottom=0.5*min(lows))
        ax.invert_xaxis()
        ax.set_xticks(hR, [f"{v:g}" for v in hR])
        ax.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
        ax.set_xlabel(r"$h / R$")
        ax.set_title(titles[k])
        ax.grid(True, which="major", color="0.9", linewidth=0.6)
    axes[0].set_ylabel("error")
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="outside lower center", ncol=len(labels), frameon=False)

    fig.savefig(os.path.join(sys.argv[1], "embedded1d3d_vessel_convergence.png"), dpi=150)
