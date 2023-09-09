import numpy as np

times_3d_reg = np.genfromtxt("logfile_pnm_3d_comparison.txt")
times_3d = np.genfromtxt("logfile_pnm_3d_comparison_no_regularization.txt")
timesteps_3d = np.diff(times_3d)
timesteps_3d_reg = np.diff(times_3d_reg)
avg_time = np.average(timesteps_3d)
avg_time_reg = np.average(timesteps_3d_reg)
print("Averaged time step for 3d comparsion case: ", avg_time)
print("Averaged time step for 3d comparison case using regularization: ", avg_time_reg)

iterations = np.genfromtxt("NewtonLog.txt")
