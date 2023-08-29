"""
This is the script for reproducing Fig. 11
"""

import numpy as np
import matplotlib
import matplotlib.pyplot as plt


pc_gloabl, num_inv_static = np.genfromtxt("static_pc_sw_pc-s-curve.txt", skip_header=1, usecols=(1,2)).T

pc_global2, num_inv_dynamic2 = np.genfromtxt("eqPoints_pcScurve_reg-0.01.txt", skip_header=1, usecols= (1,6)).T

pc_global3, num_inv_dynamic3 = np.genfromtxt("eqPoints_pcScurve_reg-0.001.txt", skip_header=1, usecols= (1,6)).T

pc_global4, num_inv_dynamic4 = np.genfromtxt("eqPoints_pcScurve_reg-0.0001.txt", skip_header=1, usecols= (1,6)).T


plt.plot(pc_gloabl, num_inv_static, marker = "x", ls = "", label = "Quasi-static solution", linewidth = 1, markersize = 4)
plt.plot(pc_global4, num_inv_dynamic4, marker = "+", ls = "", label="EFI-R, $\epsilon = 0.0001$")
plt.plot(pc_global3, num_inv_dynamic3, marker = "1", ls = "" , label="EFI-R, $\epsilon = 0.001$")
plt.plot(pc_global2, num_inv_dynamic2, marker = "2", ls = "", label="EFI-R, $\epsilon = 0.01$")

plt.xticks(fontsize = 14)
plt.yticks(fontsize = 14)
plt.xlabel("Global capillary pressure $p_c$ [Pa]", fontsize = 14)
plt.ylabel("Number of invaded throats [-]", fontsize = 14)
plt.legend(frameon=False, fontsize = 14, loc="lower right")
plt.tight_layout(rect=[0.0, 0.0, 1.0, 0.95], pad=0.4, w_pad=2.0, h_pad=2.0)
plt.savefig("Number_invaded_throats.pdf", dpi = 900)
