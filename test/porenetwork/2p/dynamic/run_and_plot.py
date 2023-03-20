#!/usr/bin/env python3

"""
We compare averaged newton iterations
and L2 error at the end of time loop. (50s)
"""
import os
import subprocess
import numpy as np
import matplotlib.pyplot as plt
from collections import defaultdict

test        = 'test_pnm_2p_1d'
pvdFileName = 'test_pnm_2p_1d.pvd'
times       = {}

intervals   = np.linspace(0.00001, 0.1, 50)

totalNewtonIterations   = []
averageNewtonIterations = []
l2Error = []
swEntry = 0.101485677973638
singlePoreVolume = 6.4e-11
massFlux = 5e-8
density = 1000
tEnd = 5
swAtPore5 = 1 -  (massFlux/density*tEnd - singlePoreVolume*4*(1-swEntry))/singlePoreVolume

extractedVtpFile = "test_pnm_2p_1d-00001.vtp"
extractedResults = "test_pnm_2p_1d-00001.csv"
homeDir = os.path.expanduser('~')
pvpythonPath = homeDir + '/ParaView-5.6.0-osmesa-MPI-Linux-64bit/bin/pvpython'
extractScript = 'extractresults.py'


swAnalytical = np.array([swEntry, swEntry, swEntry, swEntry, swAtPore5, 1, 1, 1, 1, 1])


for idx, interval in enumerate(intervals):
    subprocess.run(['./' + test]
                   + ['-TimeLoop.DtInitial', str(0.01)]
                   + ['-TimeLoop.TEnd', str(tEnd)]
                   + ['-Regularization.RegPercentage', str(interval)]
                   + ['-Problem.NonWettingMassFlux', str(massFlux)])
    iterations = np.genfromtxt('NewtonLog.txt').T
    totalNewtonIterations.append(np.sum(iterations))
    averageNewtonIterations.append(np.average(iterations))
    subprocess.run(['rm', 'NewtonLog.txt'])
    subprocess.run([pvpythonPath, extractScript, '-f', extractedVtpFile, '-p1', '0.0', '0.0', '0.0', '-p2', '4.5e-3', '0', '0', "-r", '9'])
    swNumerical = np.genfromtxt(extractedResults, skip_header=1, usecols=0, delimiter=",").T
    l2Error.append((np.square(swNumerical - swAnalytical)).mean(axis=0))


fig, ax1 = plt.subplots()

color = 'tab:red'
ax1.plot(intervals, totalNewtonIterations, label = "total newton steps", color = color)
ax1.tick_params(axis ='y', labelcolor = color)
ax1.set_xlabel("regularization interval width $\epsilon$")
ax1.set_ylabel("total newton steps [-]")


ax2 = ax1.twinx()
color = 'tab:green'
ax2.set_ylabel("$\Delta_{S_\mathrm{n}, L_2}$ [-]", color = color)
ax2.plot(intervals, l2Error, color = color)
ax2.tick_params(axis ='y', labelcolor = color)

np.savetxt("output.txt", (intervals, totalNewtonIterations, l2Error))


plt.tight_layout(rect=[0.03, 0.07, 1, 0.93], pad=0.4, w_pad=2.0, h_pad=1.0)

plt.show()
plt.savefig("Test2_Regularization.pdf", dpi = 300)
