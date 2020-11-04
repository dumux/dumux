#!/usr/bin/env python3

"""
# Script to reproduce Figure 5 of
# Koch et al (2020) JCP
# https://doi.org/10.1016/j.jcp.2020.109370
#
# The scipt also prints a table with errors and
# convergence rates
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
label = {"test_md_embedded_1d3d_1p1p_tpfatpfa_average":"LS", "test_md_embedded_1d3d_1p1p_tpfatpfa_surface":"CSS", "test_md_embedded_1d3d_1p1p_tpfatpfa_kernel":"DS"}
marker = {"test_md_embedded_1d3d_1p1p_tpfatpfa_average":"o", "test_md_embedded_1d3d_1p1p_tpfatpfa_surface":"^", "test_md_embedded_1d3d_1p1p_tpfatpfa_kernel":"d"}
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
try:
    import matplotlib
    import matplotlib.pyplot as plt

    plt.style.use('ggplot')
    font = {'family': 'sans-serif', 'weight': 'normal', 'size': 8}
    matplotlib.rc('font', **font)

    dpi = 300.0
    fig, axes = plt.subplots(1, 3, dpi=dpi, figsize=(8, 4))

    for exec, result in res.items():
        p3d, p1d, q, h = get_errors(result=result)
        hR = np.array(h)/radius
        axes[0].plot(hR, p3d, "--" + marker[exec], label=label[exec])
        axes[0].set_title(r"$\vert\vert p^\mathbb{M}_{t,e} -p^\mathbb{M}_t \vert\vert_2$")
        axes[1].plot(hR, p1d, "--" + marker[exec], label=label[exec])
        axes[1].set_title(r"$\vert\vert p_{v,e} - p_v \vert\vert_2$")
        axes[2].plot(hR, q, "--" + marker[exec], label=label[exec])
        axes[2].set_title(r"$\vert\vert q_e - q \vert\vert_2$")

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
