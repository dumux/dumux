#!/usr/bin/env python3
"""
This is a script to plot constitutive relations (Pci-Swi and Kij-Pcij)
for Pore-Network model combined with the excutable plot_constitutive_relations
complied with source code implememnted in DuMux project
! To get nice plot, set REGULARIZATIONWITHPRESSURE=1
"""
import numpy as np
import matplotlib.pyplot as plt
import subprocess

print("Now we are going to plot the constitutive laws")

swDra, pcDra = np.genfromtxt("Pc_Sw_Dra.log").T
kwinv, kninv = np.genfromtxt("K-Sw_Dra.log", usecols=(0,1)).T
swImb, pcImb = np.genfromtxt("Pc_Sw_Imb.log").T
kwImb, knImb = np.genfromtxt("K-Sw_Imb.log", usecols=(2, 3)).T

def add_arrow(line, position=None, direction='right', size=15, color=None):
    """
    add an arrow to a line.

    line:       Line2D object
    position:   x-position of the arrow. If None, mean of xdata is taken
    direction:  'left' or 'right'
    size:       size of the arrow in fontsize points
    color:      if None, line color is taken.
    """
    if color is None:
        color = line[0].get_color()

    xdata = line[0].get_xdata()
    ydata = line[0].get_ydata()

    if position is None:
        position = xdata.mean()
    # find closest index
    start_ind = np.argmin(np.absolute(xdata - position))
    if direction == 'right':
        end_ind = start_ind + 1
    else:
        end_ind = start_ind - 1

    line[0].axes.annotate('',
        xytext=(xdata[start_ind], ydata[start_ind]),
        xy=(xdata[end_ind], ydata[end_ind]),
        arrowprops=dict(arrowstyle="->", color=color),
        size=size
    )

print("First the local pc-Sw with hysteresis")
line1 = plt.plot(swDra, pcDra, color = "green", ls = "-", label="$p^{\mathrm{dr}}_{\mathrm{c,i}}$")
line2 = plt.plot(swImb, pcImb, color = "black", ls = "--", label="$p^{\mathrm{im}}_{\mathrm{c,i}}$")
plt.xlabel(r"$S_\mathrm{w,i}$")
plt.ylabel(r"$p_\mathrm{c,i}$ in Pa")
plt.xticks((0.0, 0.5, 1.0), labels= ["0", "0.5", "1"])
plt.ylim(0)
add_arrow(line1, direction='left')
add_arrow(line2, direction='right')
plt.legend(loc='upper right', frameon=False)
#plt.show()
plt.tight_layout(rect=[0.0, 0.0, 1.0, 0.95], pad=0.4, w_pad=2.0, h_pad=2.0)
plt.savefig("pore-local-pc-S.pdf", dpi = 900)
plt.clf


fig, ax = plt.subplots(dpi=900, ncols=2, nrows=1, figsize=(6, 3), sharey=True)
dynamicViscosityWetting = 1e-3
dynamicViscosityNonWetting = 1e-3
upwindTermWetting = 1/dynamicViscosityWetting
upwindTermNonWetting = 1/dynamicViscosityNonWetting

line3 = ax[0].plot(pcDra, kwinv*upwindTermWetting, color="mediumblue",  ls = "-", linewidth = 1, label= "$g^{\mathrm{dr}}_{\mathrm{w,ij}}$")
line4 = ax[0].plot(pcDra, kninv*upwindTermNonWetting, color="orange", ls = "-",  linewidth = 1, label = "$g^{\mathrm{dr}}_{\mathrm{n,ij}}$")
ax[0].set_xlabel(r"$p_\mathrm{c,ij}$ in Pa")
ax[0].set_xticks((0, 20000), labels = ["0", "20000"])
ax[0].set_ylabel(r"$g_{\alpha, ij}$ in m$^2$/(Pas)")
ax[0].legend(loc='center right', frameon=False)
add_arrow(line3, direction='left')
add_arrow(line4, direction='left')

line5 = ax[1].plot(pcImb, kwImb*upwindTermWetting, color="blue",  ls = "-",  linewidth = 1, label= "$g^{\mathrm{im}}_{\mathrm{w,ij}}$")
line6 = ax[1].plot(pcImb, knImb*upwindTermNonWetting, color="darkorange", ls = "-",  linewidth = 1, label = "$g^{\mathrm{im}}_{\mathrm{n,ij}}$")
ax[1].set_xlabel(r"$p_\mathrm{c,ij}$ in Pa")
ax[1].set_xticks((0, 20000), labels = ["0", "20000"])
ax[1].set_ylabel(r"$g_{\alpha, ij}$ in m$^2$/(Pas)")
ax[1].legend(loc='center right', frameon=False)
add_arrow(line5, direction='right')
add_arrow(line6, direction='right')

fig.tight_layout(rect=[0.0, 0.0, 1.0, 0.95], pad=0.4, w_pad=2.0, h_pad=2.0)
#plt.show()
fig.savefig("k-pc.pdf", dpi=900)
