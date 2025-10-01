#! /usr/bin/env python3
import os
import subprocess
import sys

print("=" * 60)
print("DuMux Mandel Test Results Visualization Script")
print("=" * 60)

print("\n[1/8] Importing required libraries...")
try:
    import pyvista as pv
    print("✓ PyVista already available")
except ImportError:
    print("⚠ PyVista not found, installing...")
    subprocess.check_call([sys.executable, "-m", "pip", "install", "pyvista"])
    import pyvista as pv
    print("✓ PyVista installed successfully")

import numpy as np
import matplotlib.pyplot as plt
print("✓ NumPy and Matplotlib imported")

print("\n[2/8] Setting up file paths...")
# Set the path to the pvd files
root = os.path.dirname(os.path.abspath(__file__))
flow_pvd = "test_mandel_flow.pvd"
mech_pvd = "test_mandel_mech.pvd"

print(f"  Working directory: {root}")
print(f"  Flow PVD file: {flow_pvd}")
print(f"  Mechanics PVD file: {mech_pvd}")

print("\n[3/8] Loading PVD files...")
try:
    flow_reader = pv.get_reader(os.path.join(root,flow_pvd))
    print("✓ Flow PVD file loaded successfully")
except Exception as e:
    print(f"✗ Error loading flow PVD file: {e}")
    sys.exit(1)

try:
    mech_reader = pv.get_reader(os.path.join(root,mech_pvd))
    print("✓ Mechanics PVD file loaded successfully")
except Exception as e:
    print(f"✗ Error loading mechanics PVD file: {e}")
    sys.exit(1)

print("\n[4/8] Configuring plot parameters...")
# set parameters for the plot
num_time_step = 5
total_time_steps = len(flow_reader.time_values)
print(f"  Total time steps available: {total_time_steps}")
print(f"  Number of time steps to plot: {num_time_step}")

indices = np.linspace(0, total_time_steps - 1, num_time_step, dtype=int)
selected_times = np.array(flow_reader.time_values)[indices]
print(f"  Selected time indices: {indices}")
print(f"  Selected time values: {selected_times}")

print("\n[5/8] Setting up sampling line...")
# get the position for plot over line function
flow_mesh = flow_reader.read()[0]
start_x = flow_mesh.bounds[0]
end_x = flow_mesh.bounds[1]
start_y = np.mean(flow_mesh.bounds[3:5])
end_y = np.mean(flow_mesh.bounds[3:5])
line_start = [start_x, start_y, 0]
line_end = [end_x, end_y, 0]

print(f"  Mesh bounds: x=[{start_x:.3f}, {end_x:.3f}], y=[{flow_mesh.bounds[2]:.3f}, {flow_mesh.bounds[3]:.3f}]")
print(f"  Sampling line: from [{line_start[0]:.3f}, {line_start[1]:.3f}, {line_start[2]:.3f}] to [{line_end[0]:.3f}, {line_end[1]:.3f}, {line_end[2]:.3f}]")

print("\n[6/8] Extracting pressure data...")
# extract the pressure
p = []
p_exact = []
for i, time in enumerate(selected_times):
    print(f"  Processing time step {i+1}/{len(selected_times)}: t = {time:.6f}")
    flow_reader.set_active_time_value(time)
    flow_mesh = flow_reader.read()[0]
    sampled_line = flow_mesh.sample_over_line(line_start, line_end, resolution=39)
    p.append(sampled_line.point_data["p"])
    p_exact.append(sampled_line.point_data["pExact"])
    print(f"    ✓ Sampled {len(sampled_line.point_data['p'])} points along line")
x = sampled_line.points[:,0]
p = np.array(p)
p_exact = np.array(p_exact)

print(f"✓ Data extraction complete: {p.shape[0]} time steps, {p.shape[1]} spatial points")

print("\n[7/8] Normalizing pressure data...")
# Normalize pressure by initial maximum pressure
p_initial_max = np.max(p[0])
p_normalized = p / p_initial_max
p_exact_normalized = p_exact / p_initial_max
print(f"  Initial maximum pressure: {p_initial_max:.6f}")
print(f"✓ Pressure data normalized")

print("\n[8/8] Preparing visualization...")
plt.figure(figsize=(12, 6))
selected_points = 10
point_distance = len(x) // selected_points
print(f"  Plotting {num_time_step} time steps")
print(f"  Using {selected_points} scatter points per time step for clarity")

from matplotlib.lines import Line2D
from matplotlib.patches import Patch

# Plot data without labels first
for i in range(num_time_step):
    if i == 0:
        plt.plot(x, p_normalized[i], color=f"C{i}", linewidth=2)
        print(f"    ✓ Plotted line for t={selected_times[i]:.6f}")
    else:
        plt.scatter(x[::point_distance], p_normalized[i][::point_distance], marker="x", color=f"C{i}", s=50)
        plt.plot(x, p_exact_normalized[i], color=f"C{i}", linestyle='--', alpha=0.7)
        print(f"    ✓ Plotted scatter + exact solution for t={selected_times[i]:.6f}")

# Create custom legend
legend_elements = []
# Add time entries with colors
for i in range(num_time_step):
    legend_elements.append(Patch(facecolor=f"C{i}", label=f"t={selected_times[i]:.3f}"))
# Add solution type markers
legend_elements.append(Line2D([0], [0], marker='x', color='gray', linestyle='', markersize=8, label='Numerical'))
legend_elements.append(Line2D([0], [0], color='gray', linestyle='--', label='Analytical'))

plt.legend(handles=legend_elements, bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=9)
plt.xlabel('x-coordinate')
plt.ylabel('Normalized Pressure')
plt.title('Mandel Test: Normalized Pressure Evolution Over Time')
plt.xlim(0, 100)
plt.ylim(bottom=0)
plt.grid(True, alpha=0.3)
plt.tight_layout()

print("\n[9/9] Displaying plot...")
print("✓ Plot window opened - close the window to exit")
plt.show()

print("\n" + "=" * 60)
print("Visualization complete! Check the plot window for results.")
print("=" * 60)
