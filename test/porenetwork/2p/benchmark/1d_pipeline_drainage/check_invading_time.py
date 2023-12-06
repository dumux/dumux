#!/usr/bin/env python3

"""
This is the script to reproduce Fig 8

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


extractedVtpFiles = ["test_pnm_2p_1d_oilwater_drainage-00001.vtp", "test_pnm_2p_1d_oilwater_drainage_other_implicit-00001.vtp"]
extractedResultss =  ["test_pnm_2p_1d_oilwater_drainage-00001.csv", "test_pnm_2p_1d_oilwater_drainage_other_implicit-00001.csv"]


swEntry = 0.101485677973638
singlePoreVolume = 6.4e-11
massFlux = 5e-10
density = 1000
tEnd = 5


timeToInvadeFirstPore = (1.0 - swEntry) * singlePoreVolume / (massFlux/density)

print("The first throat is invaded at: ", timeToInvadeFirstPore)
