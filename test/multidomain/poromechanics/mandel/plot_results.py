
import os
import subprocess
import sys

try:
    import pyvista as pv
except ImportError:
    subprocess.check_call([sys.executable, "-m", "pip", "install", "pyvista"])
    import pyvista as pv

import numpy as np
import matplotlib.pyplot as plt

# Set the path to the pvd files
root = os.path.dirname(os.path.abspath(__file__))
flow_pvd = "test_mandel_flow.pvd"
mech_pvd = "test_mandel_mech.pvd"

flow_reader = pv.get_reader(os.path.join(root,flow_pvd))
mech_reader = pv.get_reader(os.path.join(root,mech_pvd))

# set parameters for the plot
num_time_step = 5
indices = np.linspace(0, len(flow_reader.time_values) - 1, num_time_step, dtype=int)
selected_times = np.array(flow_reader.time_values)[indices]

# get the position for plot over line function
flow_mesh = flow_reader.read()[0]
start_x = flow_mesh.bounds[0]
end_x = flow_mesh.bounds[1]
start_y = np.mean(flow_mesh.bounds[3:5])
end_y = np.mean(flow_mesh.bounds[3:5])
line_start = [start_x, start_y, 0]
line_end = [end_x, end_y, 0]

# extract the pressure
p = []
p_exact = []
for time in selected_times:
    flow_reader.set_active_time_value(time)
    flow_mesh = flow_reader.read()[0]
    sampled_line = flow_mesh.sample_over_line(line_start, line_end, resolution=100)
    p.append(sampled_line.point_data["p"])
    p_exact.append(sampled_line.point_data["pExact"])
x = sampled_line.points[:,0]
p = np.array(p)
p_exact = np.array(p_exact)

plt.figure()
selected_points = 10
point_distance = len(x) // selected_points
for i in range(num_time_step):
    if i == 0:
        plt.plot(x, p[i], label=f"t={selected_times[i]}",color=f"C{i}")
    else:
        plt.scatter(x[::point_distance], p[i][::point_distance], label=f"t={selected_times[i]}", marker="x", color=f"C{i}")
        plt.plot(x, p_exact[i], label=f"t={selected_times[i]} exact", color=f"C{i}")
plt.legend()
plt.show()
