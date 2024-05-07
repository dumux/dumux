import os
import subprocess
import numpy as np
import matplotlib.pyplot as plt

initialRegDealta = 1e-6
timeStepSize = [0.8, 0.4, 0.2, 0.1, 0.05]
l2error = []
for idx, timeStep in enumerate(timeStepSize):
    subprocess.run(['./' + 'test_pnm_2p_1d_drainage_reg']
                   + ['-TimeLoop.MaxTimeStepSize', str(timeStep)]
                   + ['-Grid.ThroatCrossSectionShape', 'Circle']
                   + ['-Problem.RegularizationDelta', str(initialRegDealta/(2**idx))])
    l2error.append(np.genfromtxt("test_pnm_2p.log"))

plt.plot(timeStepSize, l2error)
plt.xlabel("time step size [s]")
plt.ylabel("L2-Error")
plt.show()
