#!/usr/bin/env python3

"""
We compare our regularized MFI-PNM with FI PNM.

We simulate 1d-drainage case and compare accuracy and
efficency in terms of average time step size.

We use global L2 error to measure accuracy and
total newton iterations to represent the efficiency.
"""

import os
import subprocess
import numpy as np
import matplotlib.pyplot as plt

tests        = ['test_pnm_2p_1d_oilwater_drainage', 'test_pnm_2p_1d_oilwater_drainage_other_implicit']
pvdFileNames = ['test_pnm_2p_1d_oilwater_drainage.pvd', 'test_pnm_2p_1d_oilwater_drainage_other_implicit.pvd']

swEntry = 0.101485677973638
singlePoreVolume = 6.4e-11
massFlux = 5e-8
density = 1000
tEnd = 5
swAtPore5 = 1 -  (massFlux/density*tEnd - singlePoreVolume*4*(1-swEntry))/singlePoreVolume
swAnalytical = np.array([swEntry, swEntry, swEntry, swEntry, swAtPore5, 1, 1, 1, 1, 1])
kinematicViscosity = np.array([0.1])*1e-6

# We limit the maximum time step size from 0.05 to tEnd
maxTimeStep = np.linspace(0.02, 2.5, 4)
intervals = [0.05, 0.1]

extractedVtpFiles = ["test_pnm_2p_1d_oilwater_drainage-00001.vtp", "test_pnm_2p_1d_oilwater_drainage_other_implicit-00001.vtp"]
extractedResultss =  ["test_pnm_2p_1d_oilwater_drainage-00001.csv", "test_pnm_2p_1d_oilwater_drainage_other_implicit-00001.csv"]
homeDir = os.path.expanduser('~')
pvpythonPath = homeDir + '/ParaView-5.6.0-osmesa-MPI-Linux-64bit/bin/pvpython'
extractScript = 'extractresults.py'

for mu in kinematicViscosity:
    fig, ax = plt.subplots(dpi=300, ncols=2, nrows=1, figsize=(6, 3))
    l2Error = []
    averageTimeStep = []
    totalNewtonIterations =  []
    l2ErrorReg = [[], []]
    averageTimeStepReg = [[], []]
    totalNewtonIterationsReg =  [[], []]

    for maxdt in maxTimeStep:
        subprocess.run(['./' + tests[1]]
                       + ['-TimeLoop.DtInitial', str(maxdt)]
                       + ['-TimeLoop.TEnd', str(tEnd)]
                       + ['-TimeLoop.MaxTimeStepSize', str(maxdt)]
                       + ['-Grid.ThroatCrossSectionShape', 'Circle']
                       + ['-Problem.Name', str(tests[1])]
                       + ['-Problem.NonWettingMassFlux', str(massFlux)]
                       + ['-Newton.NewtonOutputFilename', str("NewtonLog.txt")]
                       + ['-Component.LiquidKinematicViscosity', str(mu)]
                       + ['-Component.LiquidDensity', str(1000)])
        iterations = np.genfromtxt('NewtonLog.txt').T
        totalNewtonIterations.append(np.sum(iterations))
        subprocess.run(['rm', 'NewtonLog.txt'])
        subprocess.run([pvpythonPath, extractScript, '-f', extractedVtpFiles[1], '-p1', '0.0', '0.0', '0.0', '-p2', '4.5e-3', '0', '0', "-r", '9'])
        swNumerical = np.genfromtxt(extractedResultss[1], skip_header=1, usecols=0, delimiter=",").T
        l2errorSw =  np.sum(np.power((swNumerical-swAnalytical),2))/10
        l2errorSw = np.sqrt(l2errorSw)
        l2Error.append(l2errorSw)
        timestepFile = 'time_steps_' + str(tests[1]) +".txt"
        times = np.genfromtxt(timestepFile)
        timesteps = np.diff(times)
        timestep_avg = np.average(timesteps)
        averageTimeStep.append(timestep_avg)

        for idx, interval in enumerate(intervals):
            subprocess.run(['./' + tests[0]]
                           + ['-TimeLoop.DtInitial', str(maxdt)]
                           + ['-TimeLoop.TEnd', str(tEnd)]
                           + ['-TimeLoop.MaxTimeStepSize', str(maxdt)]
                           + ['-Regularization.RegPercentage', str(interval)]
                           + ['-Grid.ThroatCrossSectionShape', 'Circle']
                           + ['-Problem.Name', str(tests[0])]
                           + ['-Problem.NonWettingMassFlux', str(massFlux)]
                           + ['-Newton.NewtonOutputFilename', str("NewtonLog.txt")]
                           + ['-Component.LiquidKinematicViscosity', str(mu)]
                           + ['-Component.LiquidDensity', str(1000)])

            iterations = np.genfromtxt('NewtonLog.txt').T
            totalNewtonIterationsReg[idx].append(np.sum(iterations))
            subprocess.run(['rm', 'NewtonLog.txt'])
            subprocess.run([pvpythonPath, extractScript, '-f', extractedVtpFiles[0], '-p1', '0.0', '0.0', '0.0', '-p2', '4.5e-3', '0', '0', "-r", '9'])
            swNumerical = np.genfromtxt(extractedResultss[0], skip_header=1, usecols=0, delimiter=",").T
            l2errorSw =  np.sum(np.power((swNumerical-swAnalytical),2))/10
            l2errorSw = np.sqrt(l2errorSw)
            l2ErrorReg[idx].append(l2errorSw)
            timestepFile = 'time_steps_' + str(tests[0]) +".txt"
            times = np.genfromtxt(timestepFile)
            timesteps = np.diff(times)
            timestep_avg = np.average(timesteps)
            averageTimeStepReg[idx].append(timestep_avg)

    ax[0].plot(averageTimeStep, l2Error, marker="^", label = 'FI', linewidth = 2, markersize = 6)
    for idx, interval in enumerate(intervals):
        ax[0].plot(averageTimeStepReg[idx], l2ErrorReg[idx], marker="o", label = 'MFI, $\epsilon = $' + str(interval), linewidth = 2, markersize = 6)
    ax[0].set_xlabel("Average time step size")
    ax[0].set_ylabel("$E_{S_{w}}$ [-]")
    ax[0].set_yscale('log')
    ax[0].set_xticks([0.02, 1, 2, 2.5], label = ["0.02", "1", "2", "2.5"])
    ax[0].legend()

    ax[1].plot(averageTimeStep, totalNewtonIterations, marker="^", label = 'FI', linewidth = 2, markersize = 6)
    for idx, interval in enumerate(intervals):
        ax[1].plot(averageTimeStepReg[idx], totalNewtonIterationsReg[idx], marker="o", label = 'MFI, $\epsilon = $' + str(interval), linewidth = 2, markersize = 6)
    ax[1].set_xlabel("Average time step size")
    ax[1].set_ylabel("total newton iterations [-]")
    ax[1].set_xticks([0.02, 1, 2, 2.5], label = ["0.02", "1", "2", "2.5"])
    ax[1].legend()
    plt.subplots_adjust(wspace=0, hspace=0)
    fig.tight_layout(rect=[0.03, 0.07, 1, 0.9], pad=0.4, w_pad=2.0, h_pad=1.0)
    plt.savefig("1D_drainage_accuracy_efficency.pdf", dpi=900)

    np.savetxt("1D_comparison_time_accuracy_efficency_FI.txt", (averageTimeStep, totalNewtonIterations, l2Error))
    for idx in range(2):
        np.savetxt("1D_comparison_time_accuracy_efficency_MFI.txt" + str(idx), (averageTimeStepReg[idx], totalNewtonIterationsReg[idx], l2ErrorReg[idx]))
