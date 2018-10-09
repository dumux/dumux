#!/usr/bin/env python
"""Simple mesh generation with python + matplotlib

"""

from __future__ import absolute_import, division, print_function
import h5py
import xml.etree.cElementTree as et
from xml.dom import minidom
import matplotlib.pyplot as plt
import matplotlib.tri as tri
import numpy as np
from matplotlib.path import Path
import sys

x0 = 0.0
x1 = 20.0
y0 = 0.0
y1 = 20.0
x_steps = 36
y_steps = 36

xList = np.linspace(x0,x1,x_steps)
yList = np.linspace(y0,y1,y_steps)

x_points = []
y_points = []

for i in range(len(xList)):
    for j in range(len(yList)):
        x_points.append(xList[i])
        y_points.append(yList[j])

triang = tri.Triangulation(x_points, y_points)
triangles = triang.triangles

midpoints_x = []
midpoints_y = []

for i in range(len(triangles)):
    midpoints_x.append((x_points[triangles[i][0]]+x_points[triangles[i][1]]+x_points[triangles[i][2]])/3.0)
    midpoints_y.append((y_points[triangles[i][0]]+y_points[triangles[i][1]]+y_points[triangles[i][2]])/3.0)

################################
#
# Define parameters
#
#
# [z, h, u ,v, ks]
################################

parameters = []
nCells = len(triangles)

in_h = [0.0 for i in range(nCells)]
in_u = [0.0 for i in range(nCells)]
in_v = [0.0 for i in range(len(triangles))]
in_ks = [0.0 for i in range(nCells)]
in_z = [0.0 for i in range(nCells)]
in_fluxId = [0.0 for i in range(nCells)]
in_zoneId = [i for i in range(nCells)]
in_defaultId = [i for i in range(nCells)]

for i in range(nCells):
    x = midpoints_x[i]
    y = midpoints_y[i]

    if x < 10.0:
        in_h[i] = 4.00
    else:
        in_h[i] = 1.0

    in_z[i] = 0.0
    in_u[i] = 0.0
    in_v[i] = 0.0
    in_ks[i] = 0.00
    in_defaultId[i] = i
    in_zoneId[i] = 0.0
    in_fluxId[i] = 0.0

##################################   write a hdf5 file ##########################
xy = []
for i in range(len(x_points)):
    xy.append([x_points[i],y_points[i]])

flat_tri = np.asarray(triangles).flatten()
time = 0
f = h5py.File("init.hdf5","w")
grid = f.create_group("grids/grid_0")
grid["topology"] = flat_tri
grid.attrs['TopologyType'] = 'Triangle'
topo_dimension = "%d" %(nCells)
grid.attrs['dimension'] = topo_dimension

#save the geometry
grid.attrs["GeometryType"] = "vertices"
grid["vertices"] = xy
geo_dimension = "%d %d" %(len(xy),2)
grid.attrs["Dimensions"] = geo_dimension

group_name = "timesteps/start"
grid = f.create_group(group_name)
grid.attrs["time"] = time
grid["h"] = in_h
grid["u"] = in_u
grid["v"] = in_v
grid["ks"] = in_ks
grid["z"] = in_z
grid["zoneId"] = in_zoneId
grid["defaultId"] = in_defaultId
grid["fluxId"] = in_fluxId

f.close()

print(len(in_ks))
#### Plot the triangulation.
plt.figure()
plt.gca().set_aspect('equal')
plt.triplot(triang, 'bo-')
plt.plot(midpoints_x,midpoints_y,"ro")
plt.title('triplot of Delaunay triangulation')
#plt.axes([19.0,24.0,0.0,4.0])
plt.show()
