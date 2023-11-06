import numpy as np
import matplotlib.pyplot as plt

t, sw = np.genfromtxt("logfile_test_pnm_2pnc_reg.txt", usecols= (0, 1)).T
t_square, sw_square = np.genfromtxt("../logfile_test_pnm_2pnc_reg.txt", usecols= (0, 1)).T

plt.plot(t, sw, label = "Circular throat")
plt.plot(t_square, sw_square, label = "Square throat")
plt.xlabel("Time [s]")
plt.ylabel("Saturation [-]")
plt.legend()
plt.show()
