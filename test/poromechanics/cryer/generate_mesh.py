"""
Generate the 1/8-sphere mesh for the Cryer consolidation benchmark.

Surfaces:
  Physical surface 1  =  x=0 symmetry plane  (yz quarter-disk)
  Physical surface 2  =  y=0 symmetry plane  (xz quarter-disk)
  Physical surface 3  =  z=0 symmetry plane  (xy quarter-disk)
  Physical surface 4  =  sphere surface (drained outer boundary)

Usage:
  python3 generate_mesh.py [--lc 0.12] [--output sphere_eighth.msh]
"""
import argparse
import math
import gmsh

parser = argparse.ArgumentParser()
parser.add_argument("--lc", type=float, default=0.12, help="Mesh characteristic length")
parser.add_argument("--output", type=str, default="sphere_eighth.msh")
args = parser.parse_args()

gmsh.initialize()
gmsh.model.add("sphere_eighth")

# Build 1/8 sphere in the positive octant using OCC built-in
# angle1=0, angle2=pi/2 → polar from equator to north pole
# angle3=pi/2           → azimuth 0..90°
R0 = 1.0
gmsh.model.occ.addSphere(0, 0, 0, R0,
    angle1=0.0, angle2=math.pi/2, angle3=math.pi/2)
gmsh.model.occ.synchronize()

# Identify surfaces by bounding box
# After addSphere with these angles, numbering is:
#   Surface 1: sphere surface (xmax=ymax=zmax=R0, all nonzero)
#   Surface 2: z=0 plane     (zmin=zmax=0)
#   Surface 3: y=0 plane     (ymin=ymax=0)
#   Surface 4: x=0 plane     (xmin=xmax=0)
sphere_surf = []
sym_x = []
sym_y = []
sym_z = []

for _, s in gmsh.model.getEntities(dim=2):
    bb = gmsh.model.getBoundingBox(2, s)
    xmin, ymin, zmin, xmax, ymax, zmax = bb
    if abs(xmax) < 1e-6:        # x=0 plane
        sym_x.append(s)
    elif abs(ymax) < 1e-6:      # y=0 plane
        sym_y.append(s)
    elif abs(zmax) < 1e-6:      # z=0 plane
        sym_z.append(s)
    else:                        # sphere surface
        sphere_surf.append(s)

assert sphere_surf, "Could not identify sphere surface"
assert sym_x, "Could not identify x=0 symmetry plane"
assert sym_y, "Could not identify y=0 symmetry plane"
assert sym_z, "Could not identify z=0 symmetry plane"

gmsh.model.addPhysicalGroup(2, sym_x,      tag=1, name="symmetry_x")
gmsh.model.addPhysicalGroup(2, sym_y,      tag=2, name="symmetry_y")
gmsh.model.addPhysicalGroup(2, sym_z,      tag=3, name="symmetry_z")
gmsh.model.addPhysicalGroup(2, sphere_surf,tag=4, name="sphere_surface")
gmsh.model.addPhysicalGroup(3, [1],        tag=1, name="domain")

# Global mesh size
gmsh.option.setNumber("Mesh.CharacteristicLengthMin", args.lc * 0.5)
gmsh.option.setNumber("Mesh.CharacteristicLengthMax", args.lc)

gmsh.option.setNumber("Mesh.Algorithm",   5)   # Delaunay 2D
gmsh.option.setNumber("Mesh.Algorithm3D", 4)   # Frontal 3D
gmsh.option.setNumber("Mesh.Optimize",    1)

gmsh.model.mesh.generate(3)

# Write in msh2 format (DuMux ALUGrid reads this)
gmsh.option.setNumber("Mesh.MshFileVersion", 2.2)
gmsh.write(args.output)
print(f"Mesh written to {args.output}")
print(f"  Nodes:    {gmsh.model.mesh.getNodes()[0].shape[0] // 3}")

gmsh.finalize()
