#!/usr/bin/env python3
"""
2D GMSH geometry: Gas channel + pillar box + cone exit
Uses the GMSH Python API for better control over mesh refinement
"""

import gmsh
import math

# Initialize GMSH
gmsh.initialize()
gmsh.model.add("channel_with_pillars")

# ============================================================
# USER PARAMETERS
# ============================================================

channel_L = 1.0       # Total length of gas channel
channel_H = 0.1       # Height of gas channel

box_W = 0.5 * channel_L  # Width of pillar box
box_H = 0.4           # Height of pillar box

pillar_R = 0.025      # Radius of each circular pillar
pillar_nx = 6         # Number of pillars in x-direction
pillar_ny = 4         # Number of pillars in y-direction

exit_W = 0.1 * box_W
y_offset_exit = 0.1 * box_H

# Mesh sizes
lc_channel = 0.015     # Open channel regions
lc_box = 0.015         # Pillar box
lc_pillar = 0.001     # Around pillars
lc_exit = 0.005       # Cone exit

pillar_margin = 0.1

# ============================================================
# DERIVED COORDINATES
# ============================================================

x_chan_left = 0.0
x_chan_right = channel_L
y_chan_bot = 0.0
y_chan_top = channel_H

x_box_left = 0.5 * (channel_L - box_W)
x_box_right = x_box_left + box_W
y_box_top = y_chan_bot
y_box_bot = y_chan_bot - box_H

x_exit_left = 0.5 * (x_box_left + x_box_right) - 0.5 * exit_W
x_exit_right = x_exit_left + exit_W
y_exit = y_box_bot - y_offset_exit

# Pillar grid parameters
x_inner_left = x_box_left + pillar_R
x_inner_right = x_box_right - pillar_R
y_inner_top = y_box_top - pillar_R - pillar_margin
y_inner_bot = y_box_bot + pillar_R

dx_step = 3 * pillar_R
dy_step = math.sqrt(3) / 2 * dx_step

width_needed = (pillar_nx - 1) * dx_step + 0.5 * dx_step
x_offset = (x_inner_right - x_inner_left - width_needed) / 2

height_needed = (pillar_ny - 1) * dy_step
y_offset = (y_inner_top - y_inner_bot - height_needed) / 2

# ============================================================
# CREATE GEOMETRY - Use OCC only
# ============================================================

# Channel outer points
p1 = gmsh.model.occ.addPoint(x_chan_left, y_chan_bot, 0, lc_channel)
p2 = gmsh.model.occ.addPoint(x_box_left, y_chan_bot, 0, lc_box)
p3 = gmsh.model.occ.addPoint(x_box_right, y_chan_bot, 0, lc_box)
p4 = gmsh.model.occ.addPoint(x_chan_right, y_chan_bot, 0, lc_channel)
p5 = gmsh.model.occ.addPoint(x_chan_right, y_chan_top, 0, lc_channel)
p6 = gmsh.model.occ.addPoint(x_chan_left, y_chan_top, 0, lc_channel)

# Box corners
p7 = gmsh.model.occ.addPoint(x_box_left, y_box_bot, 0, lc_box)
p8 = gmsh.model.occ.addPoint(x_box_right, y_box_bot, 0, lc_box)

# Cone exit points
p9 = gmsh.model.occ.addPoint(x_exit_left, y_exit, 0, lc_exit)
p10 = gmsh.model.occ.addPoint(x_exit_right, y_exit, 0, lc_exit)

# Channel lines
l1 = gmsh.model.occ.addLine(p1, p2)
l2 = gmsh.model.occ.addLine(p3, p4)
l3 = gmsh.model.occ.addLine(p4, p5)
l4 = gmsh.model.occ.addLine(p5, p6)
l5 = gmsh.model.occ.addLine(p6, p1)

# Box interface
l6 = gmsh.model.occ.addLine(p2, p7)
l7 = gmsh.model.occ.addLine(p7, p9)
l8 = gmsh.model.occ.addLine(p9, p10)
l9 = gmsh.model.occ.addLine(p10, p8)
l10 = gmsh.model.occ.addLine(p8, p3)
l11 = gmsh.model.occ.addLine(p2, p3)

# Create surfaces
cl1 = gmsh.model.occ.addCurveLoop([l1, l11, l2, l3, l4, l5])
s1 = gmsh.model.occ.addPlaneSurface([cl1])

cl2 = gmsh.model.occ.addCurveLoop([l6, l7, l8, l9, l10, -l11])
s2 = gmsh.model.occ.addPlaneSurface([cl2])

gmsh.model.occ.synchronize()

# Create pillar disks and collect centers using OCC
pillar_centers = []
disk_surfaces = []
for j in range(pillar_ny):
    for i in range(pillar_nx):
        cx = x_inner_left + x_offset + i * dx_step + (j % 2) * 0.5 * dx_step
        cy = y_inner_bot + y_offset + j * dy_step

        pillar_centers.append((cx, cy))

        # Create disk using OpenCASCADE primitive
        disk = gmsh.model.occ.addDisk(cx, cy, 0, pillar_R, pillar_R)
        disk_surfaces.append((2, disk))

gmsh.model.occ.synchronize()

# Cut all disks from box surface at once
if len(disk_surfaces) > 0:
    ov, ov_map = gmsh.model.occ.cut([(2, s2)], disk_surfaces, removeObject=True, removeTool=True)
    gmsh.model.occ.synchronize()

gmsh.model.occ.synchronize()

# ============================================================
# SET MESH SIZES WITH DISTANCE FIELDS
# ============================================================

ramp_distance = 2.5 * pillar_R
field_list = []

# Create distance and threshold fields for each pillar
for idx, (cx, cy) in enumerate(pillar_centers):
    # Add a point at pillar center for distance field
    pt = gmsh.model.occ.addPoint(cx, cy, 0, 1.0)
    gmsh.model.occ.synchronize()

    # Distance field to this pillar
    dist_field = 10 + idx
    gmsh.model.mesh.field.add("Distance", dist_field)
    gmsh.model.mesh.field.setNumbers(dist_field, "PointsList", [pt])
    gmsh.model.mesh.field.setNumber(dist_field, "NumPointsPerCurve", 100)

    # Threshold field for linear ramp
    thresh_field = 100 + idx
    gmsh.model.mesh.field.add("Threshold", thresh_field)
    gmsh.model.mesh.field.setNumber(thresh_field, "InField", dist_field)
    gmsh.model.mesh.field.setNumber(thresh_field, "SizeMin", lc_pillar)
    gmsh.model.mesh.field.setNumber(thresh_field, "SizeMax", lc_box)
    gmsh.model.mesh.field.setNumber(thresh_field, "DistMin", 0)
    gmsh.model.mesh.field.setNumber(thresh_field, "DistMax", ramp_distance)

    field_list.append(thresh_field)

gmsh.model.occ.synchronize()

# Combine all pillar fields with Min operation
if len(field_list) > 0:
    min_field = 200
    gmsh.model.mesh.field.add("Min", min_field)
    gmsh.model.mesh.field.setNumbers(min_field, "FieldsList", field_list)
    gmsh.model.mesh.field.setAsBackgroundMesh(min_field)
# ============================================================
# PHYSICAL GROUPS
# ============================================================

gmsh.model.addPhysicalGroup(1, [l5], tag=10001, name="inlet")
gmsh.model.addPhysicalGroup(1, [l3], tag=10002, name="outlet_chan")
gmsh.model.addPhysicalGroup(1, [l4], tag=10003, name="top_wall")
gmsh.model.addPhysicalGroup(1, [l1], tag=10004, name="chan_bot_left")
gmsh.model.addPhysicalGroup(1, [l2], tag=10005, name="chan_bot_right")

gmsh.model.addPhysicalGroup(2, [s1], tag=20001, name="gas_channel")
gmsh.model.addPhysicalGroup(2, [s2], tag=20002, name="pillar_box")

# ============================================================
# GENERATE MESH
# ============================================================

gmsh.option.setNumber("Mesh.Algorithm", 6)  # Frontal-Delaunay
gmsh.option.setNumber("Mesh.CharacteristicLengthMin", lc_pillar)
gmsh.option.setNumber("Mesh.CharacteristicLengthMax", lc_channel)

gmsh.model.mesh.generate(2)

# ============================================================
# WRITE OUTPUT
# ============================================================

# Set mesh file format to version 2
gmsh.option.setNumber("Mesh.MshFileVersion", 2.2)

# Scale to meter (we had everything in mm for better control over mesh sizes, but we want output in m)
gmsh.option.setNumber("Mesh.ScalingFactor", 1e-3)

output_file = "/Users/timokoch/dumux/dumux/test/freeflow/navierstokes/2p/pores/channel_with_pillars.msh"
gmsh.write(output_file)
print(f"Mesh written to {output_file}")

gmsh.finalize()
