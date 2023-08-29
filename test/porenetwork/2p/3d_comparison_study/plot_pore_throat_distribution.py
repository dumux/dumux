"""
This is the script to reproduce Fig. 12
"""

import numpy as np
import matplotlib.pyplot as plt

from  collections import Counter
import pandas as pd

coordinationNumber, poreRadii = np.genfromtxt("3d_comparison_pore_info.txt", usecols=(4,5), skip_header=True, delimiter=",").T
throatRadii, throatLength = np.genfromtxt("3d_comparison_throat_info.txt", skip_header=True, delimiter=",").T

fig, ax = plt.subplots(dpi=300, ncols=1, nrows=3, figsize=(4, 8))



ax[0].hist(poreRadii*1e6, bins=40, edgecolor="black")
ax[0].set_xlabel("Pore radius [$\mu m$]")
ax[0].set_ylabel("# of pores")

ax[1].hist(throatRadii*1e6, bins=40, edgecolor="black")
ax[1].set_xlabel("Throat radius [$\mu m$]")
ax[1].set_ylabel("# of throats")

coordinationNumber = (np.rint(coordinationNumber)).astype(int)
sCNumber = sorted(coordinationNumber)
sorted_counted = Counter(sCNumber)
range_length = list(range(max(coordinationNumber)))
data_series = {}
for i in range_length:
    data_series[i] = 0
for key, value in sorted_counted.items():
    data_series[key] = value
data_series = pd.Series(data_series)
x_vlaues = data_series.index
ax[2].bar(x_vlaues, data_series.values)
ax[2].set_xlabel("Coordination number")
ax[2].set_ylabel("# of pores")


#plt.show()
plt.tight_layout(rect=[0.03, 0.07, 1, 0.93], pad=0.4, w_pad=2.0, h_pad=1.0)
plt.savefig("pore-throat-distribution.pdf", dpi = 900)
