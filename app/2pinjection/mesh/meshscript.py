#!/usr/bin/env python3
import math
import numpy as np
import gmsh
from geometry import *
import sys
##
# script to generate 2d quadratic mesh where multiple layers and fractures exist


######
# input
######
# add lines
only_sub_mesh = False
if len(sys.argv) > 1:
    only_sub_mesh = True

numLayer = 5
thickness = np.array([800,150,100,150,800])
width = 2000
initial_depth = 500

numFracture = 2
fractureAngle = np.array([80,80])
fractureWidth = np.array([10,10])
fractureMidPoint = np.array([(500,-1500),(1500,-1500)])

## mesh refinement
transfinite = np.array([[20,5,40,5,20],
                       [10,10,20,10,10]])
coefficient = np.array([[-1.2,1,-0.1,1,1.2],
                      [-1.2,1,1,1,1.2]])
progression = np.array([[1,1,0,1,1],
                        [1,1,1,1,1]])
# pos left = 0, right = 1
def fractureIntersection(refx, refy, width, angle, y, pos):
    deltax = (y - refy) / math.tan(angle/180*np.pi)
    if pos == "l":
        return refx + deltax - width/2
    if pos == "r":
        return refx + deltax + width/2
    raise ValueError("error in fracture calculation")

### index
numRow = numLayer + 1
numCol = 2 + 2 * numFracture
rows = np.arange(0, numRow, dtype=int)
cols = np.arange(0, numCol, dtype=int)

# geometry parameters
depth = [-np.sum(thickness[0:row]) - initial_depth for row in rows]


gmsh.initialize()
gmsh.model.add("myMesh")

# add boundary points of domain
for row in rows:
    newPoint(row,0, 0, depth[row])
    newPoint(row, cols[-1], width, depth[row])

# add fracture intersection points
for row in range(numRow):
    for fractureIdx in range(numFracture):
        fracMidx, fracMidy = fractureMidPoint[fractureIdx]
        fracAngle = fractureAngle[fractureIdx]
        fracWidth = fractureWidth[fractureIdx]
        xL = fractureIntersection(fracMidx,fracMidy,fracWidth,fracAngle,depth[row],"l")
        xR = fractureIntersection(fracMidx,fracMidy,fracWidth,fracAngle,depth[row],"r")
        newPoint(row, 2*fractureIdx+1, xL, depth[row] )
        newPoint(row, 2*fractureIdx+2, xR, depth[row])


#####
#change here for the subdomain
##
startRowIdx = 2 if only_sub_mesh else 0
endRowIdx = numRow - 2 if only_sub_mesh else numRow

startColIdx = 0 if only_sub_mesh else 0
endColIdx = 4 if only_sub_mesh else numCol

for row in range(startRowIdx, endRowIdx):
    for col in range(startColIdx, endColIdx-1):
        newLine(0,row,col)

for row in range(startRowIdx,endRowIdx-1):
    for col in range(startColIdx, endColIdx):
        newLine(1,row,col)

# add surfaces
for row in range(startRowIdx,endRowIdx-1):
    for col in range(startColIdx, endColIdx-1):
        newSurface(row,col)

for lineName,lineIdx in lines.items():
    dir = int(lineName[0])
    idx = int(lineName[2]) if dir==0 else int(lineName[1])
    meshType = "Progression" if progression[dir][idx] else "Bump"
    gmsh.model.geo.mesh.setTransfiniteCurve(lineIdx,transfinite[dir][idx],coef=coefficient[dir][idx], meshType = meshType)

for s, sIdx in surfaces.items():
    gmsh.model.geo.mesh.setTransfiniteSurface(sIdx)
    gmsh.model.geo.mesh.setRecombine(2, sIdx)

for row in range(numRow):
    for col in range(numCol):
        if row < startRowIdx or row >= endRowIdx:
            removePoint(row,col)
        elif col < startColIdx or col >= endColIdx:
            removePoint(row,col)

removePoint(0,0)
gmsh.model.geo.synchronize()
gmsh.option.setNumber("Mesh.Smoothing", 1)
gmsh.option.setNumber("Mesh.MshFileVersion", 2.2)
gmsh.model.mesh.generate(2)
gmsh.fltk.run()

filename = "fracture_sub.msh" if only_sub_mesh == True else "2pinjection.msh"
gmsh.write(filename)
gmsh.finalize()
