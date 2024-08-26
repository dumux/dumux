#!/usr/bin/env python3

"""
1d drainage case.

This scrip compares accuracy & efficiency under large time step size

for FI-N, FI-R and FI-Theta.

We use global L2 error to measure accuracy and iteration numbers to measure efficiency.
"""

import os
import subprocess
import numpy as np
import matplotlib.pyplot as plt
from collections import defaultdict

#############################################
##########     Given parameters  ############
#############################################
swEntry = 0.101485677973638 # saturation threshold for invasion
snEntry = 1 - swEntry # Sn to invade next pore
singlePoreVolume = 6.4e-11 # single pore volume in m³
massFlux = 5e-10 # injection rate in kg/s
density = 1000 # density of injection fluid in kg/m³
kinematicViscosity = 1e-7 # kinemaic viscosity in m²/s
volumeFlux = massFlux/density # volume flux of injected flux in m³/kg
tEnd = 500 # total injection time in s
numTotalPores = 10 # total pore number
maxTimeStep = [10, 20, 50] # We limit the time step size
regularizationDelta = [0.4, 0.6] # we test different regularization Delta


##############################################
####    Calculate analytical Solution  #######
##############################################
timeToInvadeOnePore = snEntry * singlePoreVolume / volumeFlux
numberInvPores = int (tEnd/timeToInvadeOnePore) # number of totally invaded pores
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
testName         = ['test_pnm_2p_1d_drainage', 'test_pnm_2p_1d_drainage_reg', 'test_pnm_2p_1d_drainage_md']
pvdFileName      = ['test_pnm_2p_1d_drainage.pvd', 'test_pnm_2p_1d_drainage_reg.pvd', 'test_pnm_2p_1d_drainage_md.pvd']
extractedVtpFileName = ['test_pnm_2p_1d_drainage-00001.vtp', 'test_pnm_2p_1d_drainage_reg-00001.vtp', 'test_pnm_2p_1d_drainage_md-00001.vtp']
extractedResultsName =  ['test_pnm_2p_1d_drainage-00001.csv', 'test_pnm_2p_1d_drainage_reg-00001.csv', 'test_pnm_2p_1d_drainage_md-00001.csv']
homeDir = os.path.expanduser('~')
pvpythonPath = homeDir + '/ParaView-5.6.0-osmesa-MPI-Linux-64bit/bin/pvpython' # path of pvpython
extractScript = 'extractresults.py'

l2Error = defaultdict(list)
averageTimeStep = defaultdict(list)
totalNewtonIterations = defaultdict(list)

##############################################
############   Run & Get Data  ###############
##############################################

for maxdt in maxTimeStep: # we loop over different maxDt
    # FI-N, here we don't set accuracy check, because for large interval (samll accruacy check value), the results are same
    print("We run the drainage case using FI-N method")
    subprocess.run(['./' + testName[0]]
                   + ['-TimeLoop.DtInitial', str(maxdt)]
                   + ['-TimeLoop.TEnd', str(tEnd)]
                   + ['-TimeLoop.MaxTimeStepSize', str(maxdt)]
                   + ['-TimeLoop.CheckPoints', str(tEnd)]
                   + ['-Grid.ThroatCrossSectionShape', 'Circle']
                   + ['-Pnm.Problem.Name', str(testName[0])]
                   + ['-Problem.NonWettingMassFlux', str(massFlux)]
                   + ['-Newton.UseLineSearch', str(True)]
                   + ['-Newton.MaxSteps', str(500)]
                   + ['-Newton.TargetSteps', str(500)]
                   + ['-Newton.MaxTimeStepDivisions', str(1)])
    # get total newton iterations
    iterations = np.genfromtxt('NewtonLog.txt')
    totalNewtonIterations['FI-N'].append(iterations)
    subprocess.run(['rm', 'NewtonLog.txt'])

    # get L2 eror
    subprocess.run([pvpythonPath, extractScript, '-f', extractedVtpFileName[0], '-p1', '0.0', '0.0', '0.0', '-p2', '4.5e-3', '0', '0', "-r", '9'])
    swNumerical = np.genfromtxt(extractedResultsName[0], skip_header=1, usecols=0, delimiter=",").T
    l2errorSw =  np.sum(np.power((swNumerical-swAnalytical),2))/10 # calculate L2 error between anayltical solution and numerical solution
    l2errorSw = np.sqrt(l2errorSw)
    l2Error['FI-N'].append(l2errorSw)

    # get average time step size
    timestepFile = 'time_steps_' + str(testName[0]) +".txt"
    times = np.genfromtxt(timestepFile)
    timesteps = np.diff(times)
    timestep_avg = np.average(timesteps)
    averageTimeStep['FI-N'].append(timestep_avg)

    # FI-R
    print("We run the drainage case using regularization")
    for idx, interval in enumerate(regularizationDelta): # idx represents index for different Regularization Delta
        print("The regularization Delta is: ", interval)
        subprocess.run(['./' + testName[1]]
                       + ['-TimeLoop.DtInitial', str(maxdt)]
                       + ['-TimeLoop.TEnd', str(tEnd)]
                       + ['-TimeLoop.MaxTimeStepSize', str(maxdt)]
                       + ['-TimeLoop.CheckPoints', str(tEnd)]
                       + ['-Grid.ThroatCrossSectionShape', 'Circle']
                       + ['-Pnm.Problem.Name', str(testName[1])]
                       + ['-Problem.RegularizationDelta', str(interval)]
                       + ['-Problem.NonWettingMassFlux', str(massFlux)]
                       + ['-Newton.UseLineSearch', str(True)]
                       + ['-Newton.MaxSteps', str(500)]
                       + ['-Newton.TargetSteps', str(500)]
                       + ['-Newton.MaxTimeStepDivisions', str(1)])

        # get total iterations
        iterations = np.genfromtxt('NewtonLog.txt')
        totalNewtonIterations['FI-R_'+ str(idx)].append(iterations)
        subprocess.run(['rm', 'NewtonLog.txt'])

        # get L2 error
        subprocess.run([pvpythonPath, extractScript, '-f', extractedVtpFileName[1], '-p1', '0.0', '0.0', '0.0', '-p2', '4.5e-3', '0', '0', "-r", '9'])
        swNumerical = np.genfromtxt(extractedResultsName[1], skip_header=1, usecols=0, delimiter=",").T
        l2errorSw =  np.sum(np.power((swNumerical-swAnalytical),2))/10 # calculate L2 error between anayltical solution and numerical solution
        l2errorSw = np.sqrt(l2errorSw)
        l2Error['FI-R_'+ str(idx)].append(l2errorSw)

        # get time step size
        timestepFile = 'time_steps_' + str(testName[1]) +".txt"
        times = np.genfromtxt(timestepFile)
        timesteps = np.diff(times)
        timestep_avg = np.average(timesteps)
        averageTimeStep['FI-R_'+ str(idx)].append(timestep_avg)

    # FI-Theta
    print("We run the drainage case using multi-domain method :")
    subprocess.run(['./' + testName[2]] + ['params_md.input']
                   + ['-TimeLoop.DtInitial', str(maxdt)]
                   + ['-TimeLoop.TEnd', str(tEnd)]
                   + ['-TimeLoop.MaxTimeStepSize', str(maxdt)]
                   + ['-TimeLoop.CheckPoints', str(tEnd)]
                   + ['-Grid.ThroatCrossSectionShape', 'Circle']
                   + ['-Pnm.Problem.Name', str(testName[2])]
                   + ['-Problem.NonWettingMassFlux', str(massFlux)]
                   + ['-Newton.UseLineSearch', str(False)]
                   + ['-Newton.MaxSteps', str(500)]
                   + ['-Newton.TargetSteps', str(500)]
                   + ['-Newton.MaxTimeStepDivisions', str(1)])

    # get total iterations
    iterations = np.genfromtxt('NewtonLog.txt')
    totalNewtonIterations['FI-Theta'].append(iterations)
    subprocess.run(['rm', 'NewtonLog.txt'])
    subprocess.run([pvpythonPath, extractScript, '-f', extractedVtpFileName[2], '-p1', '0.0', '0.0', '0.0', '-p2', '4.5e-3', '0', '0', "-r", '9'])

    # get L2 error
    swNumerical = np.genfromtxt(extractedResultsName[2], skip_header=1, usecols=0, delimiter=",").T
    l2errorSw =  np.sum(np.power((swNumerical-swAnalytical),2))/10
    l2errorSw = np.sqrt(l2errorSw)
    l2Error['FI-Theta'].append(l2errorSw)

    # get time step size
    timestepFile = 'time_steps_' + str(testName[2]) +".txt"
    times = np.genfromtxt(timestepFile)
    timesteps = np.diff(times)
    timestep_avg = np.average(timesteps)
    averageTimeStep['FI-Theta'].append(timestep_avg)


##############################################
############        Plot        ##############
##############################################

fig, ax = plt.subplots(dpi=300, ncols=2, nrows=1, figsize=(6, 3)) # two subplots, acc & eff

# 1st subplot about accuracy

ax[0].plot(averageTimeStep['FI-N'], l2Error['FI-N'], marker="^", label = 'FI-N', linewidth = 1, markersize = 6)

for idx, interval in enumerate(regularizationDelta):
    ax[0].plot(averageTimeStep['FI-R_'+ str(idx)], l2Error['FI-R_'+ str(idx)], marker="o", label = 'FI-R, $ \delta = $' + str(interval), linewidth = 1, markersize = 6)

ax[0].plot(averageTimeStep['FI-Theta'], l2Error['FI-Theta'], marker="^", label = 'FI-$\Theta$', linewidth = 1, markersize = 6)

ax[0].set_xlabel("Average time step size")
ax[0].set_ylabel("$E_{S_{w}}$ [-]")
ax[0].set_yscale('log')
ax[0].set_xticks(maxTimeStep, label = map(str, maxTimeStep))
ax[0].legend()

# 2nd subplot about efficiency

ax[1].plot(averageTimeStep['FI-N'], totalNewtonIterations['FI-N'], marker="^", label = 'FI-N', linewidth = 1, markersize = 6)

for idx, interval in enumerate(regularizationDelta):
    ax[1].plot(averageTimeStep['FI-R_'+ str(idx)], totalNewtonIterations['FI-R_'+ str(idx)], marker="o", label = 'FI-R, $ \delta = $' + str(interval), linewidth = 1, markersize = 6)

ax[1].plot(averageTimeStep['FI-Theta'], totalNewtonIterations['FI-Theta'], marker="^", label = 'FI-$\Theta$', linewidth = 1, markersize = 6)

ax[1].set_xlabel("Average time step size")
ax[1].set_ylabel("total newton iterations [-]")
ax[1].set_yscale('log')
ax[1].set_xticks(maxTimeStep, label = map(str, maxTimeStep))
ax[1].legend()
plt.subplots_adjust(wspace=0, hspace=0)
fig.tight_layout(rect=[0.03, 0.07, 1, 0.9], pad=0.4, w_pad=2.0, h_pad=1.0)
plt.show()
plt.savefig("1D_drainage_accuracy_efficency.pdf", dpi=900)
