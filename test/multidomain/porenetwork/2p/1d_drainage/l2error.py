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
maxTimeStep =  [4, 2, 1, 0.5, 0.25] # We limit the maximum time step size from 5s to 10s
regularizationDelta = [0.1, 0.01] # we test few different regularization Delta


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
testName         = ['test_pnm_2p_1d_drainage', 'test_pnm_2p_1d_drainage_reg', 'test_pnm_2p_1d_drainage_md']
pvdFileName      = ['test_pnm_2p_1d_drainage.pvd', 'test_pnm_2p_1d_drainage_reg.pvd', 'test_pnm_2p_1d_drainage_md.pvd']
homeDir = os.path.expanduser('~')
pvpythonPath = homeDir + '/ParaView-5.6.0-osmesa-MPI-Linux-64bit/bin/pvpython' # path of pvpython
extractScript = 'extractresults.py'

l2ErrorFI = [[], []]
averageTimeStepFI = [[], []]
totalNewtonIterationsFI =  [[], []]

l2ErrorRegNonReduce = [[], []]
averageTimeStepRegNonReduce = [[], []]
totalNewtonIterationsRegNonReduce =  [[], []]

l2ErrorRegReduce = [[], []]
averageTimeStepRegReduce = [[], []]
totalNewtonIterationsRegReduce =  [[], []]

l2ErrorMd = [ [] ]
averageTimeStepMd =  [ [] ]
totalNewtonIterationsMd =   [ [] ]

##############################################
############   Run & Get Data  ###############
##############################################

for itime, maxdt in enumerate(maxTimeStep):
    # Test case using classical FI
    for idx, interval in enumerate(regularizationDelta):
        print("Delta is: ", interval)
        delete_all_vtp_files()
        subprocess.run(['./' + testName[0]]
                       + ['-TimeLoop.DtInitial', str(maxdt)]
                       + ['-TimeLoop.TEnd', str(tEnd)]
                       + ['-TimeLoop.MaxTimeStepSize', str(maxdt)]
                       + ['-Grid.ThroatCrossSectionShape', 'Circle']
                       + ['-Pnm.Problem.Name', str(testName[0])]
                       + ['-Problem.NonWettingMassFlux', str(massFlux)]
                       + ['-Newton.UseLineSearch', str(True)]
                       + ['-Newton.MaxSteps', str(500)]
                       + ['-Newton.TargetSteps', str(500)]
                       + ['-InvasionState.AccuracyCriterion', str(1-interval)])
        iterations = np.genfromtxt('NewtonLog.txt')
        totalNewtonIterationsFI[idx].append(iterations)
        subprocess.run(['rm', 'NewtonLog.txt'])
        pvd_file = pvdFileName[0]
        times = get_time_steps_from_pvd(pvd_file)
        timesteps = np.diff(times)
        l2error = np.genfromtxt( str(testName[0]) + '.log' )
        l2ErrorFI[idx].append(l2error)
        timestep_avg = np.average(timesteps)
        averageTimeStepFI[idx].append(timestep_avg)

        # run regularization test case without reducing delta
        delete_all_vtp_files()
        subprocess.run(['./' + testName[1]]
                       + ['-TimeLoop.DtInitial', str(maxdt)]
                       + ['-TimeLoop.TEnd', str(tEnd)]
                       + ['-TimeLoop.MaxTimeStepSize', str(maxdt)]
                       + ['-Grid.ThroatCrossSectionShape', 'Circle']
                       + ['-Pnm.Problem.Name', str(testName[1])]
                       + ['-Problem.RegularizationDelta', str(interval)]
                       + ['-Problem.NonWettingMassFlux', str(massFlux)]
                       + ['-Newton.UseLineSearch', str(True)]
                       + ['-Newton.MaxSteps', str(500)]
                       + ['-Newton.TargetSteps', str(500)])
        iterations = np.genfromtxt('NewtonLog.txt')
        totalNewtonIterationsRegNonReduce[idx].append(iterations)
        subprocess.run(['rm', 'NewtonLog.txt'])
        pvd_file = pvdFileName[1]
        times = get_time_steps_from_pvd(pvd_file)
        timesteps = np.diff(times)
        l2error = np.genfromtxt( str(testName[1]) + '.log' )
        l2ErrorRegNonReduce[idx].append(l2error)
        timestep_avg = np.average(timesteps)
        averageTimeStepRegNonReduce[idx].append(timestep_avg)

        # run regularization test case with reduced delta
        delete_all_vtp_files()
        subprocess.run(['./' + testName[1]]
                       + ['-TimeLoop.DtInitial', str(maxdt)]
                       + ['-TimeLoop.TEnd', str(tEnd)]
                       + ['-TimeLoop.MaxTimeStepSize', str(maxdt)]
                       + ['-Grid.ThroatCrossSectionShape', 'Circle']
                       + ['-Pnm.Problem.Name', str(testName[1])]
                       + ['-Problem.RegularizationDelta', str(interval/(2**itime))]
                       + ['-Problem.NonWettingMassFlux', str(massFlux)]
                       + ['-Newton.UseLineSearch', str(True)]
                       + ['-Newton.MaxSteps', str(500)]
                       + ['-Newton.TargetSteps', str(500)])
        iterations = np.genfromtxt('NewtonLog.txt')
        totalNewtonIterationsRegReduce[idx].append(iterations)
        subprocess.run(['rm', 'NewtonLog.txt'])
        pvd_file = pvdFileName[1]
        times = get_time_steps_from_pvd(pvd_file)
        timesteps = np.diff(times)
        l2error = np.genfromtxt( str(testName[1]) + '.log' )
        l2ErrorRegReduce[idx].append(l2error)
        timestep_avg = np.average(timesteps)
        averageTimeStepRegReduce[idx].append(timestep_avg)


    # run regularization test case with reduced delta
    delete_all_vtp_files()
    subprocess.run(['./' + testName[2]]
                    + ['params_md.input ']
                    + ['-TimeLoop.DtInitial', str(maxdt)]
                    + ['-TimeLoop.TEnd', str(tEnd)]
                    + ['-TimeLoop.MaxTimeStepSize', str(maxdt)]
                    + ['-Grid.ThroatCrossSectionShape', 'Circle']
                    + ['-Pnm.Problem.Name', str(testName[2])]
                    + ['-Problem.NonWettingMassFlux', str(massFlux)]
                    + ['-Newton.UseLineSearch', str(False)]
                    + ['-Newton.MaxSteps', str(500)]
                    + ['-Newton.TargetSteps', str(500)])
    iterations = np.genfromtxt('NewtonLog.txt')
    totalNewtonIterationsMd[0].append(iterations)
    subprocess.run(['rm', 'NewtonLog.txt'])
    pvd_file = pvdFileName[1]
    times = get_time_steps_from_pvd(pvd_file)
    timesteps = np.diff(times)
    l2error = np.genfromtxt( str(testName[2]) + '.log' )
    l2ErrorMd[0].append(l2error)
    timestep_avg = np.average(timesteps)
    averageTimeStepMd[0].append(timestep_avg)

##############################################
############        Plot        ##############
##############################################

fig, ax = plt.subplots(dpi=300, ncols=2, nrows=1, figsize=(8, 4)) # two subplots, acc & eff

# 1st subplot about accuracy
for idx, interval in enumerate(regularizationDelta):
    ax[0].plot(averageTimeStepFI[idx], l2ErrorFI[idx], marker="o", label = 'FI-N, $\\tilde{\delta} = $' + str(interval), linewidth = 1, markersize = 6)
    ax[0].plot(averageTimeStepRegNonReduce[idx], l2ErrorRegNonReduce[idx], marker="^", label = 'FI-R, $\delta = $' + str(interval), linewidth = 1, markersize = 6)
    ax[0].plot(averageTimeStepRegReduce[idx], l2ErrorRegReduce[idx], marker=">", label = 'FI-R, decreased $\delta$', linewidth = 1, markersize = 6)

ax[0].plot(averageTimeStepMd[0], l2ErrorMd[0], marker="*", label = 'FI-$\Theta$', linewidth = 1, markersize = 6)

ax[0].set_xlabel("Average time step size")
ax[0].set_ylabel("$E_{S_{w}}$ [-]")
ax[0].set_yscale('log')
ax[0].set_xticks([0.5, 1, 2, 4])
ax[0].legend(loc='upper center', bbox_to_anchor=(0.47, 1.32), ncol=2, prop={'size': 8})

# 2nd subplot about efficiency
for idx, interval in enumerate(regularizationDelta):
    ax[1].plot(averageTimeStepFI[idx], totalNewtonIterationsFI[idx], marker="o", label = 'FI-N, $\\tilde{\delta} = $' + str(interval), linewidth = 1, markersize = 6)
    ax[1].plot(averageTimeStepRegNonReduce[idx], totalNewtonIterationsRegNonReduce[idx], marker="^", label = 'FI-R, $\delta = $' + str(interval), linewidth = 1, markersize = 6)
    ax[1].plot(averageTimeStepRegReduce[idx], totalNewtonIterationsRegReduce[idx], marker=">", label = 'FI-R, decreased $\delta$', linewidth = 1, markersize = 6)

ax[1].plot(averageTimeStepMd[0], totalNewtonIterationsMd[0], marker="*", label = 'FI-$\Theta$', linewidth = 1, markersize = 6)

ax[1].set_xlabel("Average time step size")
ax[1].set_ylabel("total newton iterations [-]")
ax[1].set_yscale('log')
ax[1].set_xticks([0.5, 1, 2, 4])
ax[1].legend(loc='upper center', bbox_to_anchor=(0.5, 1.32), ncol=2, prop={'size': 8})


plt.subplots_adjust(wspace=0, hspace=0)
fig.tight_layout(rect=[0.03, 0.07, 1, 0.9], pad=0.4, w_pad=2.0, h_pad=1.0)
plt.savefig("1D_drainage_accuracy_efficency.pdf", dpi=900)
