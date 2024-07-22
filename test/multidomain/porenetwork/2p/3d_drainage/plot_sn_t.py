import numpy as np
import matplotlib.pyplot as plt

time, sw = np.genfromtxt("logfile_pnm_3d.txt", usecols= (0, 1)).T
plt.plot(time, 1-sw)
plt.xlabel("Time[s]")
plt.ylabel("Sw")
plt.show()
