#!/usr/bin/env python3
# SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later


import numpy as np
import matplotlib.pyplot as plt

# normalize to maximum concentration
normalize = True

input_data = {
    "total": "clearance_tracer_amounts.dat",
}

data = {label: np.genfromtxt(name).T for (label, name) in input_data.items()}
data = {label: (dataset[0]/60.0, dataset[1]) for (label, dataset) in data.items()}

fig, ax = plt.subplots(1, 1, figsize=(4,4))

for label, (t, c) in data.items():
    if normalize:
        c = c/np.max(c)
    ax.plot(t, c, "-", linewidth=3, label=label)

    ax.set_xlabel("Time in minutes")
    if normalize:
        ax.set_ylabel("Remaining fraction of total amount of tracer")
    else:
        ax.set_ylabel("Total amount of tracer in mol")

ax.spines.right.set_visible(False)
ax.spines.top.set_visible(False)

fig.tight_layout()

plt.show()