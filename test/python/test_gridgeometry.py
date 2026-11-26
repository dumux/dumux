#!/usr/bin/env python3
# SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later


from dune.grid import structuredGrid
from dumux.discretization import GridGeometry

gridView = structuredGrid([0, 0], [1, 1], [5, 5])

gridGeometry = GridGeometry(gridView, discMethod="cctpfa")

print("The total number of scvs is {}".format(gridGeometry.numScv))
print("The total number of scvfs is {}".format(gridGeometry.numScvf))

for e in gridView.elements:
    fvGeometry = gridGeometry.boundLocalView(e)
    for scv in fvGeometry.scvs:
        print(f"scv dofIndex: {scv.dofIndex}")
        print(f"scv center: {scv.center}")
        print(f"scv volume: {scv.volume}")


def checkLocalCounts(discMethod):
    """The local view has to count what it iterates, and the two counts are not the same one."""
    localGridGeometry = GridGeometry(gridView, discMethod=discMethod)
    differ = False
    for element in gridView.elements:
        fvGeometry = localGridGeometry.boundLocalView(element)
        numScv = sum(1 for _ in fvGeometry.scvs)
        numScvf = sum(1 for _ in fvGeometry.scvfs)
        if fvGeometry.numScv != numScv:
            raise Exception(f"{discMethod}: numScv is {fvGeometry.numScv}, iterated {numScv}")
        if fvGeometry.numScvf != numScvf:
            raise Exception(f"{discMethod}: numScvf is {fvGeometry.numScvf}, iterated {numScvf}")
        differ = differ or numScv != numScvf
    if not differ:
        raise Exception(
            f"{discMethod}: no element has more faces than sub-control volumes, so "
            f"this grid cannot tell the two counts apart"
        )
    print(f"{discMethod}: local numScv/numScvf agree with the iterated entities")


checkLocalCounts("cctpfa")
checkLocalCounts("box")


def checkTemporaryLocalView(discMethod):
    """Iterating the entities of a local view nobody holds on to has to keep the view alive."""
    localGridGeometry = GridGeometry(gridView, discMethod=discMethod)
    numScv, numScvf = localGridGeometry.numScv, localGridGeometry.numScvf
    dofs = set()
    for element in gridView.elements:
        for scv in localGridGeometry.boundLocalView(element).scvs:
            dofs.add(scv.dofIndex)
        for scvf in localGridGeometry.boundLocalView(element).scvfs:
            if not (scvf.index < numScvf and scvf.insideScvIdx < numScv):
                raise Exception(f"{discMethod}: scvfs of a dropped local view read freed memory")
    if dofs != set(range(localGridGeometry.numDofs)):
        raise Exception(f"{discMethod}: the scvs of a dropped local view read freed memory")
    print(f"{discMethod}: the entities of a temporary local view are valid")


checkTemporaryLocalView("cctpfa")
checkTemporaryLocalView("box")
