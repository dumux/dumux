#!/usr/bin/env python3

"""
Multi-domain & Regularization.

1d drainage case, compare accuracy & efficency.

We use global L2 error to measure accuracy,
iteration numbers to represent efficiency.
"""

import os
import subprocess
import numpy as np
import matplotlib.pyplot as plt


#############################################
##########     Given parameters  ############
#############################################
swEntry = 0.101485677973638 # saturation below which invaison happens
snEntry = 1 - swEntry # nonwetting saturation to invade next pore
singlePoreVolume = 6.4e-11 # single pore volume in m³
massFlux = 5e-10 # injection rate in kg/s
density = 1000 # density of injection fluid in kg/m³
kinematicViscosity = 1e-7 # kinemaic viscosity in m²/s
volumeFlux = massFlux/density # volume flux of injected flux in m³/kg
tEnd = 500 # total injection time in s
numTotalPores = 10 # total pore number
maxTimeStep = [10, 20, 50] # We limit the maximum time step size from 5s to 10s
regularizationDelta = [0.4, 0.6] # we test few different regularization Delta


##############################################
####    Calculate analytical Solution  #######
##############################################
timeToInvadeOnePore = snEntry * singlePoreVolume / volumeFlux
numberInvPores = int (tEnd/timeToInvadeOnePore) # how many pores are totaly invaded
totalInjection = volumeFlux * tEnd # total volume of injected fluid
previousTotalInvPoreVolume = singlePoreVolume * numberInvPores * snEntry
snLastPore = (totalInjection - previousTotalInvPoreVolume)/singlePoreVolume
swLastPore = 1 - snLastPore # saturation in the partially invaded pore
swAnalytical = [swEntry] * numberInvPores # previous pores with critical values
swAnalytical.append(swLastPore) # last pore get partially invaded
for i in range(numTotalPores - numberInvPores - 1):
    swAnalytical.append(1.0)     # the rest pores remain fully saturated
swAnalytical = np.array(swAnalytical)


##############################################
#######    Simulation Preparation  ###########
##############################################
testName         = ['test_pnm_2p_1d_drainage_reg', 'test_pnm_2p_1d_drainage_md']
pvdFileName      = ['test_pnm_2p_1d_drainage_reg.pvd', 'test_pnm_2p_1d_drainage_md.pvd']
extractedVtpFileName = ["test_pnm_2p_1d_drainage_reg-00001.vtp", "test_pnm_2p_1d_drainage_md-00001.vtp"]
extractedResultsName =  ["test_pnm_2p_1d_drainage_reg-00001.csv", "test_pnm_2p_1d_drainage_md-00001.csv"]
homeDir = os.path.expanduser('~')
pvpythonPath = homeDir + '/ParaView-5.6.0-osmesa-MPI-Linux-64bit/bin/pvpython' # path of pvpython
extractScript = 'extractresults.py'

l2Error = []
averageTimeStep = []
totalNewtonIterations =  []
l2ErrorReg = [[], []]
averageTimeStepReg = [[], []]
totalNewtonIterationsReg =  [[], []]


##############################################
############   Run & Get Data  ###############
##############################################

# Test case using theta-regularization
for maxdt in maxTimeStep:
    print("First we run the drainage case using theta-regularization")
    for idx, interval in enumerate(regularizationDelta):
        print("The regularization Delta is: ", interval)
        subprocess.run(['./' + testName[0]]
                       + ['-TimeLoop.DtInitial', str(maxdt)]
                       + ['-TimeLoop.TEnd', str(tEnd)]
                       + ['-TimeLoop.MaxTimeStepSize', str(maxdt)]
                       + ['-TimeLoop.CheckPoints', str(tEnd)]
                       + ['-Grid.ThroatCrossSectionShape', 'Circle']
                       + ['-Pnm.Problem.Name', str(testName[0])]
                       + ['-Problem.RegularizationDelta', str(interval)]
                       + ['-Problem.NonWettingMassFlux', str(massFlux)]
                       + ['-Newton.UseLineSearch', str(True)]
                       + ['-Newton.MaxSteps', str(500)]
                       + ['-Newton.TargetSteps', str(500)]
                       + ['-Newton.MaxTimeStepDivisions', str(1)])
        iterations = np.genfromtxt('NewtonLog.txt')
        totalNewtonIterationsReg[idx].append(iterations)
        subprocess.run(['rm', 'NewtonLog.txt'])
        subprocess.run([pvpythonPath, extractScript, '-f', extractedVtpFileName[0], '-p1', '0.0', '0.0', '0.0', '-p2', '4.5e-3', '0', '0', "-r", '9'])
        swNumerical = np.genfromtxt(extractedResultsName[0], skip_header=1, usecols=0, delimiter=",").T
        l2errorSw =  np.sum(np.power((swNumerical-swAnalytical),2))/10
        l2errorSw = np.sqrt(l2errorSw)
        l2ErrorReg[idx].append(l2errorSw)
        timestepFile = 'time_steps_' + str(testName[0]) +".txt"
        times = np.genfromtxt(timestepFile)
        timesteps = np.diff(times)
        timestep_avg = np.average(timesteps)
        averageTimeStepReg[idx].append(timestep_avg)

    # Test case using multi-domain constraint
    print("Next, we run the drainage case using multi-domain method :")
    subprocess.run(['./' + testName[1]] + ['params_md.input']
                   + ['-TimeLoop.DtInitial', str(maxdt)]
                   + ['-TimeLoop.TEnd', str(tEnd)]
                   + ['-TimeLoop.MaxTimeStepSize', str(maxdt)]
                   + ['-TimeLoop.CheckPoints', str(tEnd)]
                   + ['-Grid.ThroatCrossSectionShape', 'Circle']
                   + ['-Pnm.Problem.Name', str(testName[1])]
                   + ['-Problem.NonWettingMassFlux', str(massFlux)]
                   + ['-Newton.UseLineSearch', str(False)]
                   + ['-Newton.MaxSteps', str(500)]
                   + ['-Newton.TargetSteps', str(500)]
                   + ['-Newton.MaxTimeStepDivisions', str(1)])
    iterations = np.genfromtxt('NewtonLog.txt')
    totalNewtonIterations.append(iterations)
    subprocess.run(['rm', 'NewtonLog.txt'])
    subprocess.run([pvpythonPath, extractScript, '-f', extractedVtpFileName[1], '-p1', '0.0', '0.0', '0.0', '-p2', '4.5e-3', '0', '0', "-r", '9'])
    swNumerical = np.genfromtxt(extractedResultsName[1], skip_header=1, usecols=0, delimiter=",").T
    l2errorSw =  np.sum(np.power((swNumerical-swAnalytical),2))/10
    l2errorSw = np.sqrt(l2errorSw)
    l2Error.append(l2errorSw)
    timestepFile = 'time_steps_' + str(testName[1]) +".txt"
    times = np.genfromtxt(timestepFile)
    timesteps = np.diff(times)
    timestep_avg = np.average(timesteps)
    averageTimeStep.append(timestep_avg)


##############################################
############        Plot        ##############
##############################################

fig, ax = plt.subplots(dpi=300, ncols=2, nrows=1, figsize=(6, 3)) # two subplots, acc & eff

# 1st subplot about accuracy
ax[0].plot(averageTimeStep, l2Error, marker="^", label = 'MD-Constraint', linewidth = 1, markersize = 6)
for idx, interval in enumerate(regularizationDelta):
    ax[0].plot(averageTimeStepReg[idx], l2ErrorReg[idx], marker="o", label = 'R, $\epsilon = $' + str(interval), linewidth = 1, markersize = 6)
ax[0].set_xlabel("Average time step size")
ax[0].set_ylabel("$E_{S_{w}}$ [-]")
ax[0].set_yscale('log')
ax[0].set_xticks(maxTimeStep, label = map(str, maxTimeStep))
ax[0].legend()

# 2nd subplot about efficiency
ax[1].plot(averageTimeStep, totalNewtonIterations, marker="^", label = 'MD-Constraint', linewidth = 1, markersize = 6)
for idx, interval in enumerate(regularizationDelta):
    ax[1].plot(averageTimeStepReg[idx], totalNewtonIterationsReg[idx], marker="o", label = 'R, $\epsilon = $' + str(interval), linewidth = 1, markersize = 6)
ax[1].set_xlabel("Average time step size")
ax[1].set_ylabel("total newton iterations [-]")
ax[1].set_yscale('log')
ax[1].set_xticks(maxTimeStep, label = map(str, maxTimeStep))
ax[1].legend()
plt.subplots_adjust(wspace=0, hspace=0)
fig.tight_layout(rect=[0.03, 0.07, 1, 0.9], pad=0.4, w_pad=2.0, h_pad=1.0)
plt.show()
plt.savefig("1D_drainage_accuracy_efficency.pdf", dpi=900)
