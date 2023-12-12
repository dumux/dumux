"""
The script is used to:
1. calcualte the analytical solution for 1D drainage case
2. compare the numerical results and the analytical solution
3. plot the result with matplotlib
"""

import glob, os
import re
import subprocess
import numpy as np
import matplotlib.pyplot as plt

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

##################################################
##  functions for constitutive relations      ####
##################################################
defaultPoreRadius = 2e-4
# default is water air system, 0.2mm pore radius, and cubic pores
def Pc(sw, sigma=0.0725, poreRadius=defaultPoreRadius, expFactor=6.83):
    Pc = 2.0*sigma / (poreRadius*(1.0 - np.exp(-expFactor*sw)))
    return Pc

def Sw(pc, sigma=0.0725, poreRadius=defaultPoreRadius, expFactor=6.83):
    Sw =-1.0/expFactor* np.log(1.0 - 2.0*sigma/(poreRadius*pc))
    return Sw

def dPc_dSw(sw, sigma=0.0725, poreRadius=defaultPoreRadius, expFactor=6.83):
    e = np.exp(expFactor*sw)
    dPc_dSw = -(2.0*expFactor*sigma*e) / (poreRadius*(1.0-e)*(1.0-e))
    return dPc_dSw


def PcEntry(theta=0, sigma=0.0725, radius=defaultPoreRadius*0.5):
    cosTheta = np.cos(theta)
    sinTheta = np.sin(theta)
    shapeFactor = 1.0/(4.0*np.pi)       # shape factor for circular thoat
    #shapeFactor = 1.0/16.0              # shape factor for square throat
    factorD = np.pi - 3.0*theta + 3*sinTheta*cosTheta - (cosTheta*cosTheta) / (4.0*shapeFactor)
    factorF = (1 + np.sqrt(1 + 4*shapeFactor*factorD / (cosTheta*cosTheta))) / (1.0 + 2.0*np.sqrt(np.pi*shapeFactor))
    PcEntry = sigma / radius * cosTheta * (1 + 2*np.sqrt(np.pi*shapeFactor)) * factorF
    return PcEntry


if __name__ == "__main__":
    PcValue = PcEntry()
    swEntry = Sw(PcValue)
    singlePoreVolume = 6.4e-11
    massFlux = 5e-10
    density = 1000
    tEnd = 470
    poreNumber = 10
    swAtPore5 = 1 -  (massFlux/density*tEnd - singlePoreVolume*4*(1-swEntry))/singlePoreVolume
    swAnalytical = np.array([swEntry, swEntry, swEntry, swEntry, swAtPore5, 1, 1, 1, 1, 1])
    pcAnalytical = np.array([PcEntry(), PcEntry(), PcEntry(), PcEntry(), Pc(swAtPore5), 0, 0, 0, 0, 0])
    timeToInvadeFirstPore = (1.0 - swEntry) * singlePoreVolume / (massFlux/density)
    print("The throat should be invaded when invading pore has the Saturation: ", swEntry)

    for i in range(poreNumber-1):
        print("The " + str(i+1) + "-th throat is invaded at: ", timeToInvadeFirstPore*(i+1))

    print("At time", tEnd, "s, saturation is :", swAnalytical)
    print("Now check the analytical at time ", tEnd, "s:")

    homeDir = os.path.expanduser('~')
    pvpythonPath = homeDir + '/ParaView-5.6.0-osmesa-MPI-Linux-64bit/bin/pvpython'
    extractScript = 'extractresults.py'

    Deltas = [1e-1, 1e-3, 1e-4]
    swNumerical = {}
    testName = 'test_pnm_2p_reg'
    poreid = np.linspace(1, 10, 10)
    for delta in Deltas:
        delete_all_vtp_files()
        subprocess.run(['./' + testName]
                       + ['-TimeLoop.DtInitial', str(0.1)]
                       + ['-TimeLoop.TEnd', str(tEnd)]
                       + ['-Grid.ThroatCrossSectionShape', 'Circle']
                       + ['-Pnm.Problem.Name', str(testName)]
                       + ['-Problem.NonWettingMassFlux', str(5e-10)]
                       + ['-Newton.NewtonOutputFilename', str("NewtonLog.txt")]
                       + ['-Component.LiquidKinematicViscosity', str(1e-7)]
                       + ['-Component.LiquidDensity', str(1000)]
                       + ['-Newton.MaxSteps', str(10)]
                       + ['-Newton.TargetSteps', str(4)]
                       + ['-Newton.UseLineSearch', 'false']
                       + ['-Pnm.Problem.RegularizationDelta', str(delta)])
        all_vtp = find_all_vtp_files()
        extractedVtpFile = find_last_vtp_file(all_vtp)
        extractedCsvFile = os.path.splitext(extractedVtpFile)[0]+'.csv'
        subprocess.run([pvpythonPath, extractScript, '-f', extractedVtpFile, '-p1', '0.0', '0.0', '0.0', '-p2', '4.5e-3', '0', '0', "-r", '9'])
        swNumerical[delta] = np.genfromtxt(extractedCsvFile, names=True, usecols=("S_aq"), delimiter=",").T
        print(swNumerical[delta])

        plt.plot(poreid, swNumerical[delta], label = str(delta), markersize = 10, linewidth = 1, marker= ">", alpha = 0.8)

    # plt.axhline(y = PcValue, color = 'r', linestyle = '-')
    plt.plot(poreid, swAnalytical, label = "$S_{w, ref}$", markersize = 10, linewidth = 1, color = "black", marker = "+")

    plt.xticks(fontsize = 14)
    plt.yticks(fontsize = 14)
    plt.xlabel("n-th pore", fontsize = 14)
    plt.ylabel("$S_{w,n}$", fontsize = 14)
    plt.legend(fontsize = 14, loc="lower right")
    plt.show()
    plt.savefig("Pore-by-pore-drainge.pdf", dpi=900)
