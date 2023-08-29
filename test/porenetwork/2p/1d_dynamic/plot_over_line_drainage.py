"""
This is script to run the the cases and reproduce Fig. 7
"""
import os
import subprocess

tests        = ['test_pnm_2p_1d_oilwater_drainage', 'test_pnm_2p_1d_oilwater_drainage_other_implicit']
for test in tests:
    subprocess.run(['./' + test]
                   + ['-TimeLoop.DtInitial', str(0.1)]
                   + ['-TimeLoop.TEnd', str(5)]
                   + ['-TimeLoop.MaxTimeStepSize', str(2)]
                   + ['-Grid.ThroatCrossSectionShape', 'Circle']
                   + ['-Problem.Name', str(test)]
                   + ['-Problem.NonWettingMassFlux', str(5e-8)]
                   + ['-Newton.NewtonOutputFilename', str("NewtonLog.txt")]
                   + ['-Component.LiquidKinematicViscosity', str(1e-7)]
                   + ['-Component.LiquidDensity', str(1000)]
                   + ['-Newton.MaxSteps', str(10)]
                   + ['-Newton.TargetSteps', str(4)]
                   + ['-Newton.UseLineSearch', 'false'])


import numpy as np
import matplotlib.pyplot as plt

poreid = np.linspace(1, 10, 10)
swEntry = 0.101485677973638
singlePoreVolume = 6.4e-11
massFlux = 5e-8
density = 1000
tEnd = 5
swAtPore5 = 1 -  (massFlux/density*tEnd - singlePoreVolume*4*(1-swEntry))/singlePoreVolume
swAnalytical = np.array([swEntry, swEntry, swEntry, swEntry, swAtPore5, 1, 1, 1, 1, 1])


fi_sw = [0.117609, 0.229476, 0.239035, 1, 1, 0.2119, 0.22573, 1, 1, 1]
# data can be find in output file "test_pnm_2p_1d_oilwater_drainage_other_implicit.pvd"
efi_sw = [0.101468, 0.101472, 0.101479, 0.101478, 0.687853, 1, 1, 1, 1, 1]
# data can be found in output file "test_pnm_2p_1d_oilwater_drainage.pvd"

plt.plot(poreid, swAnalytical, label = "$S_{w, ref}$", markersize = 10, linewidth = 4, color = "black", marker = "+")
plt.plot(poreid, fi_sw, label = "FI", markersize = 10, linewidth = 2, color = "blue", marker = ">", alpha = 0.8)
plt.plot(poreid, efi_sw, label = "EFI-R", markersize = 10, linewidth = 2, color ="orange", marker = "o", alpha = 0.8)
plt.xticks(fontsize = 14)
plt.yticks(fontsize = 14)
plt.xlabel("n-th pore", fontsize = 14)
plt.ylabel("$S_{w,n}$", fontsize = 14)
plt.legend(fontsize = 14, loc="lower right")
plt.show()
plt.savefig("Pore-by-pore-drainge.pdf", dpi=900)
