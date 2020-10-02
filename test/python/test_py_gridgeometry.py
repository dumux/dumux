#!/usr/bin/env python3

from dune.grid import structuredGrid
from dumux.discretization import GridGeometry

gridView = structuredGrid([0,0],[1,1],[5,5])

gridGeometry = GridGeometry(gridView, discMethod="cctpfa")
gridGeometry.update()

print("The total number of scvs is {}".format(gridGeometry.numScv()))
print("The total number of scvfs is {}".format(gridGeometry.numScvf()))

for e in gridView.elements:
    fvGeometry = gridGeometry.localView()
    fvGeometry.bind(e)

    for scv in fvGeometry.scvs():
        print("scv dofIndex: {}".format(scv.dofIndex()))
        print("scv center: {}".format(scv.center()))
        print("scv volume: {}".format(scv.volume()))
