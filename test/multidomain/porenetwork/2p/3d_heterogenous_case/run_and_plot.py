"""
This is the script to validate FI-R and FI-Theta of PNM
by comparing the capillary pressure equlibirum state with
the static PNM.
"""

import subprocess

subprocess.call(["./pnm_3d_reg"]
                + ["-Problem.Name", "pcScurve_FIR"]
                + ["-Problem.RelShiftThreshold", str(1e-5)]
                + ["-Problem.RegularizationDelta", str(0.01)]
                + ["-Newton.UseLineSearch", str(True)])

subprocess.call(["./pnm_3d_md"]
                + ["-Pnm.Problem.Name", "pcScurve_md"]
                + ["-Constraint.Problem.Name", "pcScurve_md_constraint"])

import numpy as np
import matplotlib
import matplotlib.pyplot as plt


sw_static, pc_static = np.genfromtxt("static_pc_sw_pc-s-curve.txt", skip_header=1, usecols=(0,1)).T

sw_dynamic, pc_global, pc_dynamic = np.genfromtxt("logfile_pcScurve_FIR.txt", skip_header=2, usecols=(2, 1, 5)).T
sw_eq, pc_dynamic_eq = np.genfromtxt("eqPoints_pcScurve_FIR.txt", skip_header=1, usecols= (2,1)).T

sw_md_dynamic, pc_md_dynamic = np.genfromtxt("logfile_pcScurve_md.txt", skip_header=2, usecols=(2,5)).T
sw_md_eq, pc_md_dynamic_eq = np.genfromtxt("eqPoints_pcScurve_md.txt", skip_header=1, usecols= (2,1)).T

plt.plot(sw_dynamic, pc_global, ls = "-.", label="$p_{n,inlet} - p_{w, outlet}$")
plt.plot(sw_static, pc_static, marker = "x", ls = "", label = "Quasi-static solution", linewidth = 1, markersize = 4)
plt.plot(sw_eq, pc_dynamic_eq, marker = "+", ls = "", label = "Equlibirum state, $\Theta = 0.01$", linewidth = 1, markersize = 4)
plt.plot(sw_dynamic, pc_dynamic, lw = 1, alpha = 0.5, label="Dynamic $p_c$, FI-R, $\Theta = 0.01$")

plt.plot(sw_md_eq, pc_md_dynamic_eq, marker = ">", ls = "", label = "Equlibirum state, FI-$\Theta$", linewidth = 1, markersize = 4)
plt.plot(sw_md_dynamic, pc_md_dynamic, lw = 1, alpha = 0.5, label="Dynamic $p_c$, FI-$\Theta$")

plt.xticks(fontsize = 14)
plt.yticks(fontsize = 14)
plt.xlabel("Saturation $S_w$ [-]", fontsize = 14)
plt.ylabel("Capillary pressure $p_c$ [Pa]", fontsize = 14)
plt.legend(frameon=False, fontsize = 14)
plt.tight_layout(rect=[0.0, 0.0, 1.0, 0.95], pad=0.4, w_pad=2.0, h_pad=2.0)
plt.show()
