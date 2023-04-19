#!/usr/bin/env python3
# SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later


import sys
import numpy as np
import subprocess

# Run benchmark M2.2 evaporation scenarios
subprocess.run(["./test_richards_benchmark_tpfa", "params_evaporation_sand.input"])
subprocess.run(["./test_richards_benchmark_tpfa", "params_evaporation_loam1.input"])
subprocess.run(["./test_richards_benchmark_tpfa", "params_evaporation_loam2.input"])
subprocess.run(["./test_richards_benchmark_tpfa", "params_evaporation_clay.input"])

# Create equally spaced plots to improve test robustness.
# In the case of slightly different time steps this makes sure
# we always compare the same number of data points
for refDataFile in ["rate_analytical_sand.dat", "rate_analytical_loam1.dat", "rate_analytical_loam2.dat", "rate_analytical_clay.dat"]:
    tRef, rateRef = np.genfromtxt(refDataFile).T
    # sample output data at equally spaced times using linear interpolation
    tEquallySpaced = np.linspace(np.min(tRef), np.max(tRef), 50, endpoint=True)

    sampledAnaRate = np.interp(tEquallySpaced, tRef, rateRef)
    sanitizedOutput = np.vstack((tEquallySpaced, sampledAnaRate)).T
    np.savetxt(refDataFile, sanitizedOutput)

    numericDataFile = refDataFile.replace("analytical", "actual")
    t, rate = np.genfromtxt(numericDataFile).T
    sampledRate = np.interp(tEquallySpaced, t, rate)
    sanitizedOutput = np.vstack((tEquallySpaced, sampledRate)).T
    np.savetxt(numericDataFile, sanitizedOutput)

try:
    import matplotlib.pyplot as plt
except ImportError:
    print("Matplotlib not found so no result plots will be created.")
    sys.exit(0)

# Create Fig. 5abcd of Vanderborght 2005
fig, axes = plt.subplots(nrows=2, ncols=2, figsize=(8, 8))

axes[0,0].set_ylabel(r"$E_\mathrm{act}$ in (cm day${}^{-1})$")
axes[1,0].set_ylabel(r"$E_\mathrm{act}$ in (cm day${}^{-1})$")
axes[1,0].set_xlabel(r"t (days)")
axes[1,1].set_xlabel(r"t (days)")

axes[0,0].set_ylim([0.0, 0.1])
axes[1,0].set_ylim([0.0, 0.3])

plots = {"sand": (0,0), "loam1": (0,1), "loam2": (1,0), "clay": (1,1)}
for soil, index in plots.items():
    t, e = np.genfromtxt("rate_actual_{}.dat".format(soil)).T
    t, e = t, e*0.1
    axes[index].plot(t, e, "k", label="DuMux")
    t, e = np.genfromtxt("rate_analytical_{}.dat".format(soil)).T
    t, e = t, e*0.1
    axes[index].plot(t, e, marker="x", linestyle='none', label="Benchmark")
    axes[index].legend()

plt.savefig("benchmark_evaporation.png")
