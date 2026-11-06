#!/usr/bin/env python3
# SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

"""
# Script to reproduce Figure 5 of
# Koch et al (2020) JCP
# https://doi.org/10.1016/j.jcp.2020.109370
#
# The script also prints a table with errors and
# convergence rates
"""

import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.realpath(__file__)), "../../../../../bin/testing"))
from dumux_convergencetest import convergenceRates, runLevels

executables = ["test_md_embedded_1d3d_1p1p_tpfatpfa_average", "test_md_embedded_1d3d_1p1p_tpfatpfa_surface", "test_md_embedded_1d3d_1p1p_tpfatpfa_kernel"]
method = {"test_md_embedded_1d3d_1p1p_tpfatpfa_average":1, "test_md_embedded_1d3d_1p1p_tpfatpfa_surface":3, "test_md_embedded_1d3d_1p1p_tpfatpfa_kernel":5}
label = {"test_md_embedded_1d3d_1p1p_tpfatpfa_average":"LS", "test_md_embedded_1d3d_1p1p_tpfatpfa_surface":"CSS", "test_md_embedded_1d3d_1p1p_tpfatpfa_kernel":"DS"}
marker = {"test_md_embedded_1d3d_1p1p_tpfatpfa_average":"o", "test_md_embedded_1d3d_1p1p_tpfatpfa_surface":"^", "test_md_embedded_1d3d_1p1p_tpfatpfa_kernel":"d"}
cells = [10, 20, 40, 80]
radius = 0.1

# length of the vessel elements, the vessel has unit length and half as many cells as the tissue per direction
h = [1.0/(c//2) for c in cells]

def get_errors(executable):
    name = executable + "_figure"
    rows = runLevels(f"./{executable} params.input -SpatialParams.Radius {radius} -Vtk.EnableVtkOutput false -Problem.Name {name}",
                     [f"{c//2} {c}" for c in cells], "-Vessel.Grid.Cells {0} -Tissue.Grid.Cells '{1} {1} {0}'",
                     f"{name}-errors.csv")
    return ([r["L2Rel(p_3d)"] for r in rows], [r["L2Rel(p_1d)"] for r in rows], [r["L2Rel(q)"] for r in rows])

res = {e: get_errors(e) for e in executables}

# produce output suitable for latex table
# Table 1
table1 = [["", "", "", "", "", "", ""] for i in range(len(cells)+1)]
table2 = [["", "", "", "", "", "", ""] for i in range(len(cells)+1)]
table3 = [["", "", "", "", "", "", ""] for i in range(len(cells)+1)]

for exec, (p3d, p1d, q) in res.items():
    rates3d = convergenceRates(h, p3d)
    rates1d = convergenceRates(h, p1d)
    ratesq = convergenceRates(h, q)

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
# create plot from paper (Figure 5)
#############################################
try:
    import matplotlib
    import matplotlib.pyplot as plt

    plt.style.use('ggplot')
    font = {'family': 'sans-serif', 'weight': 'normal', 'size': 8}
    matplotlib.rc('font', **font)

    dpi = 300.0
    fig, axes = plt.subplots(1, 3, dpi=dpi, figsize=(8, 4))

    for exec, (p3d, p1d, q) in res.items():
        hR = np.array(h)/radius
        axes[0].plot(hR, p3d, "--" + marker[exec], label=label[exec])
        axes[0].set_title(r"$|| p^\mathbb{M}_{t,e} -p^\mathbb{M}_t ||_2$")
        axes[1].plot(hR, p1d, "--" + marker[exec], label=label[exec])
        axes[1].set_title(r"$|| p_{v,e} - p_v ||_2$")
        axes[2].plot(hR, q, "--" + marker[exec], label=label[exec])
        axes[2].set_title(r"$|| q_e - q ||_2$")

    x = np.linspace(np.min(hR), np.max(hR), 10)
    axes[0].plot(x, np.power(x*radius*0.4, 1.5), "--k", label=r"$\Delta$ 1.5")
    axes[0].plot(x, np.power(x*radius*0.3, 2), "-.k", label=r"$\Delta$ 2")
    axes[1].plot(x, np.power(x*radius*0.4, 1.5), "--k", label=r"$\Delta$ 1.5")
    axes[1].plot(x, np.power(x*radius*0.3, 2), "-.k", label=r"$\Delta$ 2")
    axes[2].plot(x, np.power(x*radius*2.5, 2), "--k", label=r"$\Delta$ 2")
    axes[2].plot(x, np.power(x*radius*1.2, 2.5), "-.k", label=r"$\Delta$ 2.5")

    for ax in axes:
        ax.set_xscale("log")
        ax.set_xlabel("$h/r_v$")
        ax.set_xlim([np.max(hR), np.min(hR)])
        ax.xaxis.set_minor_formatter(plt.NullFormatter())
        ax.set_yscale("log")
        ax.legend()

    fig.tight_layout(rect=[0.03, 0.07, 1, 0.93], pad=0.4, w_pad=2.0, h_pad=1.0)
    fig.savefig("rates_r01.pdf")

except ImportError:
    print("Skipping plot: matplotlib has not been found.")
