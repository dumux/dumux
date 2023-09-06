#!/usr/bin/env python3

"""
We compare total newton iterations, averaged iterations
and averaged time steps
"""
import os
import subprocess
import numpy as np
import matplotlib.pyplot as plt

MFIcase = "pnm_3d_comparison"
FIcase = "pnm_3d_comparison_no_regularization"

maxTimeSteps = [5e-3, 1e-3, 5e-4, 1e-4]
intervals = [0.001, 0.01, 0.05, 0.1]

# 1st comparsion study, influence of maxDT for both schemes

for maxDt in maxTimeSteps:
    subprocess.call(["./" + MFIcase,
                 "-Problem.Name", str(MFIcase) + "_maxDt_" + str(maxDt),
                 "-TimeLoop.MaxTimeStepSize", str(maxDt),
                 "-Newton.NewtonOutputFilename", "NewtonLog_MFI_maxDt_" + str(maxDt) + ".txt",
                 "-Newton.NewtonOverview", "NewtonOverview_MFI_maxDt_" + str(maxDt) + ".txt"])

    subprocess.call(["./" + FIcase,
                 "-Problem.Name", str(FIcase) + "_maxDt_" + str(maxDt),
                 "-TimeLoop.MaxTimeStepSize", str(maxDt),
                 "-Newton.NewtonOutputFilename", "NewtonLog_FI_maxDt_" + str(maxDt) + ".txt",
                 "-Newton.NewtonOverview", "NewtonOverview_FI_maxDt_" + str(maxDt) + ".txt"])


# 2nd comparison study, influence of interval and accuracy

for epsilon in intervals:
    subprocess.call(["./" + MFIcase,
             "-Problem.Name", str(MFIcase) + "_epsilon_" + str(epsilon),
             "-Regularization.RegPercentage", str(epsilon),
             "-Newton.NewtonOutputFilename", "NewtonLog_MFI_epsilon_" + str(epsilon) + ".txt",
             "-Newton.NewtonOverview", "NewtonOverview_MFI_epsilon_" + str(epsilon) + ".txt"])

subprocess.call(["./" + FIcase,
         "-Problem.Name", str(FIcase) + "_no_limit",
         "-Newton.NewtonOutputFilename", "NewtonLog_MFI_epsilon_no_limit.txt",
         "-Newton.NewtonOverview", "NewtonOverview_MFI_epsilon_no_limit.txt"])
