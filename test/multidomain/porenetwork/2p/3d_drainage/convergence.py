"""
This script does convergence study for general 3d case
"""

import glob, os
import re
import subprocess
import numpy as np
import matplotlib.pyplot as plt
import vtk
from vtk.util.numpy_support import vtk_to_numpy

##################################################
##  functions to process output vtp files      ###
##################################################
def find_all_vtp_files():
    vtp_file_lists = []
    for file in glob.glob("*.vtp"):
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

def get_Sw_from_vtp(vtp_file_name):
    reader = vtk.vtkXMLPolyDataReader()
    reader.SetFileName(vtp_file_name)
    reader.Update()
    polydata = reader.GetOutput()
    sw_vtk = polydata.GetPointData().GetArray('S_aq')
    sw = vtk_to_numpy(sw_vtk)
    return sw

if __name__ == "__main__":
    # generate grid
    if os.path.exists('lattice_network.dgf'):
        print("use the exisiting file in your folder")
    else:
        subprocess.call(['python3', 'generate_lattice_grid.py'])
    # input parameters
    tEnd = 20
    # We first try few Detla and Dt to see if we can observe convergence behavior
    swNumerical = {}
    L2Error = []
    testName = 'pnm_3d_drainage_reg'
    startingDelta = 1e-3
    startingMaxDt = 1
    Deltas = []
    MaxDts = []
    sRef = np.array([])
    L2Error = 10
    ### Step 1 : Get Referecne solution numerically ##
    Thresholdbias = 1e-4
    minRefineTime = 4
    maxRefineTime = 10
    i = 0
    while ((L2Error > Thresholdbias and i < maxRefineTime) or i < minRefineTime):
        delta = startingDelta/(2**i)
        maxDt = startingMaxDt/(2**i)
        delete_all_vtp_files()
        subprocess.run(['./' + testName]
                       + ['-TimeLoop.DtInitial', str(maxDt)]
                       + ['-TimeLoop.TEnd', str(tEnd)]
                       + ['-TimeLoop.MaxTimeStepSize', str(maxDt)]
                       + ['-Grid.File', 'lattice_network.dgf']
                       + ['-Grid.ThroatCrossSectionShape', 'Square']
                       + ['-Pnm.Problem.Name', str(testName)]
                       + ['-Problem.RegularizationDelta', str(delta)]
                       + ['-Problem.NonWettingMassFlux', str(5e-8)]
                       + ['-Newton.UseLineSearch', str(True)]
                       + ['-Newton.MaxSteps', str(10)]
                       + ['-Newton.TargetSteps', str(4)]
                       + ['-Newton.MaxTimeStepDivisions', str(20)])
        all_vtp = find_all_vtp_files()
        extractedVtpFile = find_last_vtp_file(all_vtp)
        swNumerical[delta] = get_Sw_from_vtp(extractedVtpFile)
        Deltas.append(delta)
        MaxDts.append(maxDt)
        if (i > 0):
            L2Error = np.sqrt(np.sum(np.power((swNumerical[delta] - swNumerical[delta*2]),2))/50)
        sRef = swNumerical[delta]
        i += 1

    print("After ", (i-1), "-th Refinement, we get SRef: ", sRef, "Error with last refinement is: ", L2Error)

    np.savetxt('reference.txt', sRef)
    # ask user if he wants to see the change of bias?
    user_input = input("Do you want to plot the change of bias? (yes/no): ")
    if user_input.lower() == "yes":
        Delta1s = Deltas[:-1] # list without last element
        Delta2s = Deltas[1:]  # list without 1st element
        L2ErrorChange = []
        for delta1, delta2 in zip(Delta1s, Delta2s):
            l2errorSw =  np.sum(np.power((swNumerical[delta1] - swNumerical[delta2]),2))/50
            l2errorSw = np.sqrt(l2errorSw)
            L2ErrorChange.append(l2errorSw)


        plt.plot( np.arange( len(L2ErrorChange) ), L2ErrorChange)
        plt.xlabel("Refinement Time")
        plt.ylabel("L2 bias compared to last refinement")
        plt.show()

    # ## Step 2: Plot the convergence ##
    # sRef = np.genfromtxt('reference.txt')  ## in case the ref solution has been obtained
    DeltaRange = [0.2/(2**i) for i in range(10)]
    print(DeltaRange)

    ConvergenceL2Error = []
    for idx, delta in enumerate(DeltaRange):
        delete_all_vtp_files()
        subprocess.run(['./' + testName]
                       + ['-TimeLoop.DtInitial', str(4.0/(2**idx))]
                       + ['-TimeLoop.TEnd', str(tEnd)]
                       + ['-TimeLoop.MaxTimeStepSize', str(4.0/(2**idx))]
                       + ['-Grid.File', 'lattice_network.dgf']
                       + ['-Grid.ThroatCrossSectionShape', 'Square']
                       + ['-Pnm.Problem.Name', str(testName)]
                       + ['-Problem.RegularizationDelta', str(delta)]
                       + ['-Problem.NonWettingMassFlux', str(5e-8)]
                       + ['-Newton.UseLineSearch', str(True)]
                       + ['-Newton.MaxSteps', str(10)]
                       + ['-Newton.TargetSteps', str(4)]
                       + ['-Newton.MaxTimeStepDivisions', str(20)])
        all_vtp = find_all_vtp_files()
        extractedVtpFile = find_last_vtp_file(all_vtp)
        swNumerical[delta] = get_Sw_from_vtp(extractedVtpFile)
        L2ErrorConvergence = np.sqrt(np.sum(np.power((swNumerical[delta] - sRef),2))/50)
        ConvergenceL2Error.append(L2ErrorConvergence)
    plt.plot(np.arange(len(DeltaRange)), ConvergenceL2Error, marker="o", label = 'L2 Error', linewidth = 1, markersize = 6)
    plt.xlabel('Refinement time')
    plt.ylabel('L2 error')
    plt.yscale('log')
    plt.show()
