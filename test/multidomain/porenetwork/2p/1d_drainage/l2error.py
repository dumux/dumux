#!/usr/bin/env python3

"""
Multi-domain & Regularization.

1d drainage case, compare accuracy & efficency.

We use global L2 error to measure accuracy,
iteration numbers to represent efficiency.
"""

import glob, os
import re
import os
import subprocess
import numpy as np
import matplotlib.pyplot as plt
import xml.etree.ElementTree as ET

##################################################
##  functions to process output vtp files      ###
##################################################
def find_all_vtp_files_withPrefix(prefix):
    vtp_file_lists = []
    for file in glob.glob(prefix + "*.vtp"):
        vtp_file_lists.append(file)
    return vtp_file_lists

def extract_number(f):
    s = re.findall("(\d+).vtp",f)
    return (int(s[0]) if s else -1,f)

def find_last_vtp_file(vtp_file_lists):
    return max(vtp_file_lists,key=extract_number)

def delete_all_vtp_files():
    for item in os.listdir():
        if item.endswith(".vtp"):
            os.remove(item)

def get_time_steps_from_pvd(pvd_file):
    time_steps = []
    tree = ET.parse(pvd_file)
    root = tree.getroot()

    for collection in root.findall('.//DataSet'):
        timestep_attr = collection.attrib.get('timestep')
        if timestep_attr is not None:
            time_step = float(timestep_attr)
            time_steps.append(time_step)

    return time_steps

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
tEnd = 400 # total injection time in s
numTotalPores = 10 # total pore number
maxTimeStep = [1, 0.5, 0.25] # We limit the maximum time step size from 5s to 10s
regularizationDelta = [0.1, 0.2] # we test few different regularization Delta


##############################################
####    Calculate analytical Solution  #######
##############################################
def analyticalSolution(time):
    timeToInvadeOnePore = snEntry * singlePoreVolume / volumeFlux
    numberInvPores =  np.floor(time/timeToInvadeOnePore).astype(int) # how many pores are totaly invaded
    totalInjection = volumeFlux * time # total volume of injected fluid
    previousTotalInvPoreVolume = singlePoreVolume * numberInvPores * snEntry
    snLastPore = (totalInjection - previousTotalInvPoreVolume)/singlePoreVolume
    swLastPore = 1 - snLastPore # saturation in the partially invaded pore
    swAnalytical = [swEntry] * numberInvPores # previous pores with critical values
    swAnalytical.append(swLastPore) # last pore get partially invaded
    for i in range(numTotalPores - numberInvPores - 1):
        swAnalytical.append(1.0)     # the rest pores remain fully saturated
    swAnalytical = np.array(swAnalytical)

    return swAnalytical


##############################################
#######    Simulation Preparation  ###########
##############################################
testName         = ['test_pnm_2p_1d_drainage_reg', 'test_pnm_2p_1d_drainage_md']
pvdFileName      = ['test_pnm_2p_1d_drainage_reg.pvd', 'test_pnm_2p_1d_drainage_md.pvd']
#extractedVtpFileName = ["test_pnm_2p_1d_drainage_reg-00001.vtp", "test_pnm_2p_1d_drainage_md-00001.vtp"]
#extractedResultsName =  ["test_pnm_2p_1d_drainage_reg-00001.csv", "test_pnm_2p_1d_drainage_md-00001.csv"]
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
for itime, maxdt in enumerate(maxTimeStep):
    print("First we run the drainage case using theta-regularization")
    for idx, interval in enumerate(regularizationDelta):
        print("The regularization Delta is: ", interval)
        delete_all_vtp_files()
        subprocess.run(['./' + testName[0]]
                       + ['-TimeLoop.DtInitial', str(maxdt)]
                       + ['-TimeLoop.TEnd', str(tEnd)]
                       + ['-TimeLoop.MaxTimeStepSize', str(maxdt)]
                       + ['-Grid.ThroatCrossSectionShape', 'Circle']
                       + ['-Pnm.Problem.Name', str(testName[0])]
                       + ['-Problem.RegularizationDelta', str(interval/(2**itime))]
                       + ['-Problem.NonWettingMassFlux', str(massFlux)]
                       + ['-Newton.UseLineSearch', str(True)]
                       + ['-Newton.MaxSteps', str(500)]
                       + ['-Newton.TargetSteps', str(500)]
                       + ['-Newton.MaxTimeStepDivisions', str(1)])
        iterations = np.genfromtxt('NewtonLog.txt')
        totalNewtonIterationsReg[idx].append(iterations)
        subprocess.run(['rm', 'NewtonLog.txt'])

        pvd_file = pvdFileName[0]
        times = get_time_steps_from_pvd(pvd_file)
        print("times are :     ", times)
        timesteps = np.diff(times)
        vtpFiles = find_all_vtp_files_withPrefix(testName[0])
        print("unsorted vtp files: ", vtpFiles)
        vtpFiles.sort(key=lambda x: extract_number(x))
        print("sorted vtp files:", vtpFiles)

        assert len(vtpFiles) == len(times), "Number of files must match"

        l2errorSw = 0
        for i in range(1, len(vtpFiles)):
            subprocess.run([pvpythonPath, extractScript, '-f', vtpFiles[i], '-p1', '0.0', '0.0', '0.0', '-p2', '4.5e-3', '0', '0', "-r", '9', "-of", "output", "-v", '0'])
            swNumerical = np.genfromtxt("output.csv", skip_header=1, usecols=0, delimiter=",").T
            swAnalytical = analyticalSolution(times[i])
            l2errorSw += timesteps[i-1]*(np.sum(np.power((swNumerical-swAnalytical),2))/10)

        l2errorSw = np.sqrt(l2errorSw)
        l2ErrorReg[idx].append(l2errorSw)
        timestep_avg = np.average(timesteps)
        averageTimeStepReg[idx].append(timestep_avg)

    # Test case using multi-domain constraint
    print("Next, we run the drainage case using multi-domain method :")
    delete_all_vtp_files()
    subprocess.run(['./' + testName[1]] + ['params_md.input']
                   + ['-TimeLoop.DtInitial', str(maxdt)]
                   + ['-TimeLoop.TEnd', str(tEnd)]
                   + ['-TimeLoop.MaxTimeStepSize', str(maxdt)]
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

    pvd_file = pvdFileName[1]
    times = get_time_steps_from_pvd(pvd_file)
    timesteps = np.diff(times)
    vtpFiles = find_all_vtp_files_withPrefix(testName[1])
    vtpFiles.sort(key=lambda x: extract_number(x))

    assert len(vtpFiles) == len(times), "Number of files must match"

    l2errorSw = 0
    for i in range(1, len(vtpFiles)):
        subprocess.run([pvpythonPath, extractScript, '-f', vtpFiles[i], '-p1', '0.0', '0.0', '0.0', '-p2', '4.5e-3', '0', '0', "-r", '9', "-of", "output", "-v", '0'])
        swNumerical = np.genfromtxt("output.csv", skip_header=1, usecols=0, delimiter=",").T
        swAnalytical = analyticalSolution(times[i])
        l2errorSw += timesteps[i-1]*(np.sum(np.power((swNumerical-swAnalytical),2))/10)

    l2errorSw = np.sqrt(l2errorSw)
    l2Error.append(l2errorSw)
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
