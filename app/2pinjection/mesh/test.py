import gmsh
from geometry import *
gmsh.initialize()
gmsh.model.add("test")

newPoint(0,0,0,0)
newPoint(0,1,0,1)
newPoint(1,0,1,0)
newPoint(1,1,1,1)

newLine(0,0,0)
newLine(0,1,0)
newLine(1,0,0)
newLine(1,0,1)

newSurface(0,0)
showGeometry()

gmsh.finalize()
