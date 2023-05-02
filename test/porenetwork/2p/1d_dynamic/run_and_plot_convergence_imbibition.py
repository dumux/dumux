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

test        = 'test_pnm_2p_1d_oilwater_imbibition'
pvdFileName = 'test_pnm_2p_1d.pvd'
times       = {}

intervals   = np.linspace(0.00001, 0.01, 50)

massFlux = 5e-8
density = 1000
tEnd = 0.02
kinematicViscosity = np.array([10])*1e-6

extractedVtpFile = "test_pnm_2p_1d-00001.vtp"
extractedResults = "test_pnm_2p_1d-00001.csv"
homeDir = os.path.expanduser('~')
pvpythonPath = homeDir + '/ParaView-5.6.0-osmesa-MPI-Linux-64bit/bin/pvpython'
extractScript = 'extractresults.py'

for mu in kinematicViscosity:
    totalNewtonIterations   = []
    averageNewtonIterations = []
    l2Error = []
    subprocess.run(['./' + test]
                   + ['-TimeLoop.DtInitial', str(1e-4)]
                   + ['-TimeLoop.TEnd', str(tEnd)]
                   + ['-TimeLoop.MaxTimeStepSize', str(1e-7)] # use small maxdt to get Ref solution
                   + ['-Regularization.RegPercentage', str(1e-10)]
                   + ['-Problem.NonWettingMassFlux', str(massFlux)]
                   + ['-Newton.NewtonOutputFilename', str("NewtonLog.txt")]
                   + ['-Component.LiquidKinematicViscosity', str(mu)]
                   + ['-Component.LiquidDensity', str(1000)])

    subprocess.run(['rm', 'NewtonLog.txt'])
    subprocess.run([pvpythonPath, extractScript, '-f', extractedVtpFile, '-p1', '0.0', '0.0', '0.0', '-p2', '4.5e-3', '0', '0', "-r", '9'])
    swAnalytical = np.genfromtxt(extractedResults, skip_header=1, usecols=0, delimiter=",").T

    for idx, interval in enumerate(intervals):
        subprocess.run(['./' + test]
                       + ['-TimeLoop.DtInitial', str(1e-3)]
                       + ['-TimeLoop.TEnd', str(tEnd)]
                       + ['-TimeLoop.MaxTimeStepSize', str(tEnd)] # no limit for max dt
                       + ['-Regularization.RegPercentage', str(interval)]
                       + ['-Problem.NonWettingMassFlux', str(massFlux)]
                       + ['-Newton.NewtonOutputFilename', str("NewtonLog.txt")]
                       + ['-Component.LiquidKinematicViscosity', str(mu)]
                       + ['-Component.LiquidDensity', str(1000)])

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

    np.savetxt("1D_imbibition_output_"+str(mu)+".txt", (intervals, totalNewtonIterations, l2Error))


    plt.tight_layout(rect=[0.03, 0.07, 1, 0.93], pad=0.4, w_pad=2.0, h_pad=1.0)
    plt.savefig("1D_Regularization_Imbibition_"+str(mu)+".pdf", dpi = 300)
