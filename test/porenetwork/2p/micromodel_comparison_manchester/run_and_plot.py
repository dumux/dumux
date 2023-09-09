#!/usr/bin/env python3

"""
We compare total newton iterations, averaged iterations
and averaged time steps
"""
import os
import subprocess
import numpy as np
import matplotlib.pyplot as plt

executables = ["pnm_3d_comparison_manchester", "pnm_3d_comparison_no_regularization_manchester"]
newtonFiles = ["NewtonLog.txt", "NewtonLog_no_reg.txt"]
newtonOverview = ["NewtonOverview.txt", "NewtonOverview_no_reg.txt"]

# remove log file
for nf in newtonFiles:
    subprocess.run(["rm", str(nf)])

# run both cases with output file
for e, nf, no in zip(executables, newtonFiles, newtonOverview):
    subprocess.call(["./" + e,
                     "-Problem.Name", str(e),
                     "-Newton.NewtonOutputFilename", str(nf),
                     "-Newton.NewtonOverview", str(no)])
