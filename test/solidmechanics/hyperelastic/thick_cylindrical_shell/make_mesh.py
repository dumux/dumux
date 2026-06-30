# SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later
#
# Generate a tetrahedral mesh of half a thick cylindrical shell (x <= 0).
#
# Reproduces the geometry of the OpenGeoSys F-bar locking benchmark
# "Thick cylindrical shell under pressure":
#   inner radius ri = 0.008 m, outer radius ro = 0.010 m (thickness 0.002 m),
#   height H = 0.015 m, only the x <= 0 half of the annulus is modelled.
#
# A structured (r, theta, z) lattice of hexahedra is built and each hexahedron
# is split into 6 tetrahedra using a consistent (conforming) Kuhn/Freudenthal
# subdivision along the main diagonal. The output is a Gmsh-2.2 ASCII file that
# Dune::ALUGrid<3,3,simplex,conforming> can read.
#
# Usage: python3 make_mesh.py [nr nt nz]   (default 2 20 10, matching OGS)
import sys
import numpy as np
import meshio

ri, ro, H = 0.008, 0.010, 0.015
nr, nt, nz = (int(a) for a in sys.argv[1:4]) if len(sys.argv) >= 4 else (2, 20, 10)


def node_id(i, j, k):
    return i * (nt + 1) * (nz + 1) + j * (nz + 1) + k


# nodes on the (r, theta, z) lattice; theta from 90 deg (top, +y) over 180 deg
# (left, -x) to 270 deg (bottom, -y), so all points have x <= 0
points = np.empty(((nr + 1) * (nt + 1) * (nz + 1), 3))
for i in range(nr + 1):
    r = ri + (ro - ri) * i / nr
    for j in range(nt + 1):
        theta = 0.5 * np.pi + np.pi * j / nt
        for k in range(nz + 1):
            points[node_id(i, j, k)] = (r * np.cos(theta), r * np.sin(theta), H * k / nz)

# 6-tetrahedra Kuhn subdivision of each lattice cube (shares diagonal v0-v7)
KUHN = [(0, 1, 3, 7), (0, 3, 2, 7), (0, 2, 6, 7),
        (0, 6, 4, 7), (0, 4, 5, 7), (0, 5, 1, 7)]

tets = []
for i in range(nr):
    for j in range(nt):
        for k in range(nz):
            v = [node_id(i + (c & 1), j + ((c >> 1) & 1), k + ((c >> 2) & 1)) for c in range(8)]
            for a, b, c, d in KUHN:
                t = [v[a], v[b], v[c], v[d]]
                p = points[t]
                # ensure positive orientation (positive signed volume)
                if np.dot(np.cross(p[1] - p[0], p[2] - p[0]), p[3] - p[0]) < 0:
                    t[2], t[3] = t[3], t[2]
                tets.append(t)

mesh = meshio.Mesh(points=points, cells=[("tetra", np.array(tets, dtype=int))])
meshio.write("grid.msh", mesh, file_format="gmsh22", binary=False)
print(f"wrote grid.msh: {len(points)} nodes, {len(tets)} tets ({nr}x{nt}x{nz} hexes)")
