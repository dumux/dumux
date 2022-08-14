#!/usr/bin/env python3

import sys
import numpy as np
import subprocess

# Run benchmark M3.1 root scenario
subprocess.run(["./test_1p_rootbenchmark_tpfa", "-Grid.Refinement", "5"])

try:
    import matplotlib.pyplot as plt
except ImportError:
    print("Matplotlib not found so no result plots will be created.")
    sys.exit(0)

# Create Fig. of Schnepf et at 2020 benchmark
fig, axes = plt.subplots(nrows=1, ncols=1, figsize=(5, 5), dpi=72)
ax = axes

# set labels
ax.set_xlabel(r"root water pressure head (cm)")
ax.set_ylabel(r"depth (cm)")

depth = np.genfromtxt("depth.txt")
head = np.genfromtxt("pressure_head.txt")
ax.plot(head, depth)

textOutput = ""
textOutput += ",".join(f"{i:.10e}" for i in depth) + "\n"
textOutput += ",".join(f"{i:.10e}" for i in head) + "\n"

with open("root.txt", "w") as out:
    out.write(textOutput)

fig.tight_layout()
fig.savefig("root.png")
