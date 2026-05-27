#! /usr/bin/env python3
import os
import subprocess
import sys

print("=" * 60)
print("DuMux Mandel Test Results Visualization Script")
print("=" * 60)

print("\n[1/9] Importing required libraries...")
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
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
print("✓ NumPy and Matplotlib imported")

print("\n[2/9] Setting up file paths...")
# Set the path to the pvd files: prefer the directory containing the script,
# but fall back to the current working directory (e.g. when running from a
# build directory where the PVDs live).
flow_pvd = "test_mandel_flow.pvd"
mech_pvd = "test_mandel_mech.pvd"
script_dir = os.path.dirname(os.path.abspath(__file__))
if os.path.isfile(os.path.join(script_dir, flow_pvd)):
    root = script_dir
else:
    root = os.getcwd()

print(f"  Working directory: {root}")
print(f"  Flow PVD file: {flow_pvd}")
print(f"  Mechanics PVD file: {mech_pvd}")

print("\n[3/9] Loading PVD files...")
try:
    flow_reader = pv.get_reader(os.path.join(root, flow_pvd))
    print("✓ Flow PVD file loaded successfully")
except Exception as e:
    print(f"✗ Error loading flow PVD file: {e}")
    sys.exit(1)

try:
    mech_reader = pv.get_reader(os.path.join(root, mech_pvd))
    print("✓ Mechanics PVD file loaded successfully")
except Exception as e:
    print(f"✗ Error loading mechanics PVD file: {e}")
    sys.exit(1)

print("\n[4/9] Configuring plot parameters...")
# set parameters for the plot
num_time_step = 5
total_time_steps = len(flow_reader.time_values)
print(f"  Total time steps available: {total_time_steps}")
print(f"  Number of time steps to plot: {num_time_step}")

indices = np.linspace(0, total_time_steps - 1, num_time_step, dtype=int)
selected_times = np.array(flow_reader.time_values)[indices]
print(f"  Selected time indices: {indices}")
print(f"  Selected time values: {selected_times}")

print("\n[5/9] Setting up sampling line...")
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

print("\n[6/9] Extracting pressure and displacement data...")
# extract the pressure and displacement
p = []
p_exact = []
u = []
u_exact = []
resolution = 39
sampled_flow = None
sampled_mech = None
for i, time in enumerate(selected_times):
    print(f"  Processing time step {i+1}/{len(selected_times)}: t = {time:.6f}")

    flow_reader.set_active_time_value(time)
    flow_mesh = flow_reader.read()[0]
    sampled_flow = flow_mesh.sample_over_line(line_start, line_end, resolution=resolution)
    p.append(sampled_flow.point_data["p"])
    p_exact.append(sampled_flow.point_data["pExact"])

    mech_reader.set_active_time_value(time)
    mech_mesh = mech_reader.read()[0]
    sampled_mech = mech_mesh.sample_over_line(line_start, line_end, resolution=resolution)
    u.append(sampled_mech.point_data["u"])
    u_exact.append(sampled_mech.point_data["uExact"])

    print(f"    ✓ Sampled {len(sampled_flow.point_data['p'])} flow points, "
          f"{len(sampled_mech.point_data['u'])} mech points along line")

assert sampled_flow is not None and sampled_mech is not None, "No time steps were sampled"
x = sampled_flow.points[:, 0]
x_mech = sampled_mech.points[:, 0]
p = np.array(p)
p_exact = np.array(p_exact)
u = np.array(u)
u_exact = np.array(u_exact)

print(f"✓ Data extraction complete: pressure {p.shape}, displacement {u.shape}")

print("\n[7/9] Normalizing pressure data...")
# Normalize pressure by initial maximum pressure
p_initial_max = np.max(p[0])
p_normalized = p / p_initial_max
p_exact_normalized = p_exact / p_initial_max
print(f"  Initial maximum pressure: {p_initial_max:.6f}")
print(f"✓ Pressure data normalized")

print("\n[8/9] Preparing visualization...")
selected_points = 10
point_distance = max(1, len(x) // selected_points)
point_distance_mech = max(1, len(x_mech) // selected_points)
print(f"  Plotting {num_time_step} time steps into 3 separate figures (p, u_x, u_y)")

# Helper that draws one quantity into a fresh figure
def make_figure(x_arr, num_data, exact_data, stride, title, ylabel,
                ylim_bottom_zero=False):
    fig, ax = plt.subplots(figsize=(10, 6))
    for i in range(num_time_step):
        if i == 0:
            # initial condition: single line for clarity
            ax.plot(x_arr, num_data[i], color=f"C{i}", linewidth=2)
        else:
            ax.scatter(x_arr[::stride], num_data[i][::stride],
                       marker="x", color=f"C{i}", s=50)
            ax.plot(x_arr, exact_data[i], color=f"C{i}",
                    linestyle='--', alpha=0.7)

    legend_elements: list = [
        Patch(facecolor=f"C{i}", label=f"t={selected_times[i]:.0f}s")
        for i in range(num_time_step)
    ]
    legend_elements.append(Line2D([0], [0], marker='x', color='gray',
                                  linestyle='', markersize=8, label='Numerical'))
    legend_elements.append(Line2D([0], [0], color='gray',
                                  linestyle='--', label='Analytical'))

    ax.set_xlabel('x-coordinate [m]')
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    ax.set_xlim(x_arr.min(), x_arr.max())
    if ylim_bottom_zero:
        ax.set_ylim(bottom=0)
    ax.grid(True, alpha=0.3)
    ax.legend(handles=legend_elements, bbox_to_anchor=(1.02, 1),
              loc='upper left', fontsize=9)
    fig.tight_layout()
    return fig

# Figure 1: normalized pressure
fig_p = make_figure(x, p_normalized, p_exact_normalized,
                    point_distance,
                    'Mandel Test: Normalized Pressure Evolution Over Time',
                    'p / p$_0^{max}$ [-]',
                    ylim_bottom_zero=True)

# Figure 2: x-displacement (component 0)
fig_ux = make_figure(x_mech, u[..., 0], u_exact[..., 0],
                     point_distance_mech,
                     'Mandel Test: Horizontal Displacement u$_x$ Over Time',
                     'u$_x$ [m]')

# Figure 3: y-displacement (component 1)
fig_uy = make_figure(x_mech, u[..., 1], u_exact[..., 1],
                     point_distance_mech,
                     'Mandel Test: Vertical Displacement u$_y$ Over Time',
                     'u$_y$ [m]')

# Save figures next to the PVD inputs
for fig, name in [(fig_p, "results_pressure.png"),
                  (fig_ux, "results_ux.png"),
                  (fig_uy, "results_uy.png")]:
    out_path = os.path.join(root, name)
    fig.savefig(out_path, dpi=150, bbox_inches='tight')
    print(f"  ✓ Saved {out_path}")

print("\n[9/9] Displaying plots...")
print("✓ Plot windows opened - close all windows to exit")
plt.show()

print("\n" + "=" * 60)
print("Visualization complete! Check the plot window for results.")
print("=" * 60)
