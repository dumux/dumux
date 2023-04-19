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

test        = 'test_pnm_2p_1d_oilwater'
pvdFileName = 'test_pnm_2p_1d.pvd'

swEntry = 0.101485677973638
singlePoreVolume = 6.4e-11
massFlux = 5e-8
density = 1000
tEnd = 5
swAtPore5 = 1 -  (massFlux/density*tEnd - singlePoreVolume*4*(1-swEntry))/singlePoreVolume
kinematicViscosity = np.array([10, 1, 0.1])*1e-6

extractedVtpFile = "test_pnm_2p_1d-00001.vtp"
extractedResults = "test_pnm_2p_1d-00001.csv"
homeDir = os.path.expanduser('~')
pvpythonPath = homeDir + '/ParaView-5.6.0-osmesa-MPI-Linux-64bit/bin/pvpython'
extractScript = 'extractresults.py'

swAnalytical = np.array([swEntry, swEntry, swEntry, swEntry, swAtPore5, 1, 1, 1, 1, 1])

iterations = np.genfromtxt('NewtonLog.txt').T
print("total iterations: ", np.sum(iterations))
print("\n")
subprocess.run(['rm', 'NewtonLog.txt'])
subprocess.run([pvpythonPath, extractScript, '-f', extractedVtpFile, '-p1', '0.0', '0.0', '0.0', '-p2', '4.5e-3', '0', '0', "-r", '9'])
swNumerical = np.genfromtxt(extractedResults, skip_header=1, usecols=0, delimiter=",").T
print("L2error", np.square(swNumerical - swAnalytical).mean(axis=0))

times = np.genfromtxt("time_steps.txt")
subprocess.run(['rm', 'time_steps.txt'])
timesteps = np.diff(times)
print("Average time step size:", np.average(timesteps))
