"""
This is the script for reproducing Fig. 10
"""

import subprocess

static_test_name  = "pc_sw_static"
dynamic_test_name = "pc_sw_dynamic"

regIntervals = [1e-2, 1e-3, 1e-4]
subprocess.call(["./" + static_test_name])

for regInterval in regIntervals:
    subprocess.call(["./" + dynamic_test_name]
                + ["-Problem.Name", "pcScurve_reg-" + str(regInterval)]
                + ["-Problem.RelShiftThreshold", str(1e-5)]
                + ["-Regularization.RegPercentage", str(regInterval)]
                + ["-Newton.NewtonOutputFilename", "NewtonLog_reg-"+str(regInterval)+".txt"]
                + ["-Newton.NewtonOverview", "NewtonOverview_reg-"+str(regInterval)+".txt"])

import numpy as np
import matplotlib
import matplotlib.pyplot as plt


sw_static, pc_static = np.genfromtxt("static_pc_sw_pc-s-curve.txt", skip_header=1, usecols=(0,1)).T

sw2, pc_global2, pc_dynamic2 = np.genfromtxt("logfile_pcScurve_reg-0.01.txt", skip_header=2, usecols=(2, 1, 5)).T
sw2_eq, pc_dynamic2_eq = np.genfromtxt("eqPoints_pcScurve_reg-0.01.txt", skip_header=1, usecols= (2,1)).T

sw3, pc_global3, pc_dynamic3 = np.genfromtxt("logfile_pcScurve_reg-0.001.txt", skip_header=2, usecols=(2, 1, 5)).T

sw4, pc_global4, pc_dynamic4 = np.genfromtxt("logfile_pcScurve_reg-0.0001.txt", skip_header=2, usecols=(2, 1, 5)).T


plt.plot(sw4, pc_global4, ls = "-.", label="$p_{n,inlet} - p_{w, outlet}$")

plt.plot(sw_static, pc_static, marker = "x", ls = "", label = "Quasi-static solution", linewidth = 1, markersize = 4)
plt.plot(sw2_eq, pc_dynamic2_eq, marker = "+", ls = "", label = "Equlibirum state, $\epsilon = 0.01$", linewidth = 1, markersize = 4)
plt.plot(sw4, pc_dynamic4, lw = 1, alpha = 0.5, label="Dynamic $p_c$, EFI-R, $\epsilon = 0.0001$")
plt.plot(sw3, pc_dynamic3, lw = 1, alpha = 1, label="Dynamic $p_c$, EFI-R, $\epsilon = 0.001$")
plt.plot(sw2, pc_dynamic2, lw = 1, alpha = 0.5, label="Dynamic $p_c$, EFI-R, $\epsilon = 0.01$")

plt.xticks(fontsize = 14)
plt.yticks(fontsize = 14)
plt.xlabel("Saturation $S_w$ [-]", fontsize = 14)
plt.ylabel("Capillary pressure $p_c$ [Pa]", fontsize = 14)
plt.legend(frameon=False, fontsize = 14)
plt.tight_layout(rect=[0.0, 0.0, 1.0, 0.95], pad=0.4, w_pad=2.0, h_pad=2.0)
plt.savefig("PcS.pdf", dpi = 900)
