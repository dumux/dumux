#!/usr/bin/env python3

"""
We compare averaged newton iterations
and L2 error at the end of time loop. (5s)
"""

import os
import subprocess
import numpy as np
import matplotlib.pyplot as plt
from collections import defaultdict

test        = 'test_pnm_2p_1d_oilwater_drainage'
pvdFileName = 'test_pnm_2p_1d.pvd'
times       = {}

intervals   = np.linspace(1e-4, 1e-1, 5)

swEntry = 0.101485677973638
singlePoreVolume = 6.4e-11
massFlux = 5e-8
density = 1000
tEnd = 5
swAtPore5 = 1 -  (massFlux/density*tEnd - singlePoreVolume*4*(1-swEntry))/singlePoreVolume
Ms = np.array([10, 1, 0.1])
kinematicViscosity = np.array([10, 1, 0.1])*1e-6

extractedVtpFile = "test_pnm_2p_1d_oilwater_drainage-00001.vtp"
extractedResults = "test_pnm_2p_1d_oilwater_drainage-00001.csv"
homeDir = os.path.expanduser('~')
pvpythonPath = homeDir + '/ParaView-5.6.0-osmesa-MPI-Linux-64bit/bin/pvpython'
extractScript = 'extractresults.py'


swAnalytical = np.array([swEntry, swEntry, swEntry, swEntry, swAtPore5, 1, 1, 1, 1, 1])

fig, ax1 = plt.subplots(dpi=300, ncols=3, nrows=1, figsize=(12, 4))

for i, m, mu in zip (range(3), Ms, kinematicViscosity):
    totalNewtonIterations   = []
    averageNewtonIterations = []
    l2Error = []
    for idx, interval in enumerate(intervals):
        subprocess.run(['./' + test]
                       + ['-TimeLoop.DtInitial', str(0.1)]
                       + ['-TimeLoop.TEnd', str(tEnd)]
                       + ['-TimeLoop.MaxTimeStepSize', str(tEnd)] # no limit for max dt
                       + ['-Regularization.RegPercentage', str(interval)]
                       + ['-Problem.NonWettingMassFlux', str(massFlux)]
                       + ['-Newton.NewtonOutputFilename', str("NewtonLog.txt")]
                       + ['-Grid.ThroatCrossSectionShape', 'Circle']
                       + ['-Newton.MaxSteps', str(20)]
                       + ['-Newton.TargetSteps', str(10)]
                       + ['-Component.LiquidKinematicViscosity', str(mu)]
                       + ['-Component.LiquidDensity', str(1000)])

        iterations = np.genfromtxt('NewtonLog.txt').T
        totalNewtonIterations.append(np.sum(iterations))
        averageNewtonIterations.append(np.average(iterations))
        subprocess.run(['rm', 'NewtonLog.txt'])
        subprocess.run([pvpythonPath, extractScript, '-f', extractedVtpFile, '-p1', '0.0', '0.0', '0.0', '-p2', '4.5e-3', '0', '0', "-r", '9'])
        swNumerical = np.genfromtxt(extractedResults, skip_header=1, usecols=0, delimiter=",").T
        l2errorSw =  np.sum(np.power((swNumerical-swAnalytical),2))/10
        l2errorSw = np.sqrt(l2errorSw)
        l2Error.append(l2errorSw)



    color = 'tab:red'
    ax1[i].plot(intervals, totalNewtonIterations,  marker = '^', markersize = 6, color = color, linewidth = 2)
    ax1[i].tick_params(axis ='y', labelcolor = color, labelsize = 14)
    ax1[i].tick_params(axis ='x', labelsize = 14)
    ax1[i].set_xlabel("$\epsilon$", fontsize = 14)
    ax1[i].set_ylabel("iterations [-]", fontsize = 14, color = color)
    ax2 = ax1[i].twinx()
    color = 'tab:green'
    ax2.plot(intervals, l2Error,  marker = 'o', markersize = 6, color = color, label = "$E_{s_w}$", linewidth = 2)
    ax2.tick_params(axis ='y', labelcolor = color, labelsize = 14)
    ax2.tick_params(axis ='x', labelsize = 14)
    y2label = ax2.set_ylabel("$E_{S_w}$ [-]", color = color, fontsize = 14)
    props = dict(boxstyle='round',  alpha=1, facecolor = "white")
    ax1[i].text(0.35, 0.95, "M = " + str(m), transform=ax1[i].transAxes, fontsize = 14 ,
        verticalalignment='top', bbox=props)
    np.savetxt("1D_drainage_output_"+str(mu)+".txt", (intervals, totalNewtonIterations, l2Error))

plt.subplots_adjust(wspace=0, hspace=0)
plt.tight_layout(rect=[0.0, 0.0, 1, 1], pad=0.4, w_pad=2.0, h_pad=1.0)
# plt.show()
plt.savefig("1D_Regularization_Drainage.pdf", dpi = 300)
