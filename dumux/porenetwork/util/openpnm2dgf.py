# SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

"""
Creates a DUNE Grid Format file (*.dgf)
from an OpenPNM porenetwork file (*.pnm)
"""

import argparse
import os
import copy
import openpnm as op  # pylint: disable=import-error
import numpy as np
from openpnm.topotools import trim  # pylint: disable=import-error


def load(filename):
    """Create workspace instance"""

    workspace = op.Workspace()
    workspace.clear()
    project = workspace.load_project(filename=filename)
    print(project)
    return project["net"]


def removeBoundaryThroats(net):
    """DOC ME!"""

    pores = net.pores("boundary")
    net["throat.boundary"] = False

    for idx, throat in enumerate(net["throat.conns"]):
        if net["pore.boundary"][throat[0]] or net["pore.boundary"][throat[1]]:
            net["throat.boundary"][idx] = True

    throats = net.throats("boundary")
    cnBoundary = net["throat.conns"][throats]

    net["pore.delete"] = False
    net["pore.delete"][pores] = True

    for label in [
        "pore.xmin",
        "pore.xmax",
        "pore.ymin",
        "pore.ymax",
        "pore.zmin",
        "pore.zmax",
    ]:
        try:
            for idx, throat in enumerate(cnBoundary):
                poreVals = [net[label][throat[0]], net[label][throat[1]]]
                val = np.max(poreVals)  # is True if any of the two pores' label is True
                # make sure both pores get same label
                net[label][throat[0]] = val
                net[label][throat[1]] = val
        except KeyError:
            print(f"{label} not found. Skipping.")

    trim(network=net, throats=throats)
    trim(network=net, pores=net["pore.delete"])


def getDualGridCoordinationNumber(net):
    """DOC ME!"""

    numPoresTotal = len(net["pore.coords"])
    coordNumVoid = np.zeros(numPoresTotal, dtype=np.int64)
    coordNumSolid = np.zeros(numPoresTotal, dtype=np.int64)

    for voidThroat in net["throat.conns"][net["throat.void_void"]]:
        coordNumVoid[voidThroat[0]] += 1
        coordNumVoid[voidThroat[1]] += 1

    for solidThroat in net["throat.conns"][net["throat.solid_solid"]]:
        coordNumSolid[solidThroat[0]] += 1
        coordNumSolid[solidThroat[1]] += 1

    return coordNumVoid, coordNumSolid


def sanitizeDualNetwork(net, splitSinglePores):
    """DOC ME!"""
    # pylint: disable=too-many-locals,too-many-branches,too-many-statements

    print("Num. void throats before dual grid sanitation", np.sum(net["throat.void_void"]))
    print(net)

    # calculate coordination number for void-void and solid-solid connections
    coordNumVoid, coordNumSolid = getDualGridCoordinationNumber(net)

    poreKeys = [key for key in net if key.startswith("pore.")]
    throatKeys = [key for key in net if key.startswith("throat.")]

    # HACK We want to add information for the new fake pores
    # and throats but it is not clear how to change the entries in net.
    # We therefore replace the net object
    # by a dictionary which behaves similarly.
    # This should be improved!
    if splitSinglePores:
        newNet = {}
        for key in net:
            newNet[key] = net[key]
        net = newNet

    singlePores = []

    for domain in ["void", "solid"]:
        if domain == "void":
            coordinationNumber = coordNumVoid
            # label = 22
            # domainType = 0
        else:
            coordinationNumber = coordNumSolid
            # label = 33
            # domainType = 1

        # identify pores only connect to pores of other phase and
        # add an artificial throat or delete those pores
        for vIdx, coordNum in enumerate(coordinationNumber):
            if net["pore." + domain][vIdx] == 1 and coordNum < 1:
                if domain == "void":
                    print(
                        "Void pore",
                        vIdx,
                        "is only connected to solid pores ...",
                        "Splitting" if splitSinglePores else "Deleting",
                    )
                else:
                    print(
                        "Solid pore",
                        vIdx,
                        "is only connected to void pores ...",
                        "Splitting" if splitSinglePores else "Deleting",
                    )
                singlePores.append(vIdx)

        if splitSinglePores:
            op.topotools.add_boundary_pores(
                network=net,
                pores=singlePores,
                apply_label="fake_" + domain,
                offset=[1e-10, 1e-10, 1e-10],
            )
            singlePores = []  # reset
            net["pore." + domain][net["pore.fake_" + domain]] = True
            net["throat." + domain][net["throat.fake_" + domain]] = True

            # add some articial data
            for throat in net.throats("fake_" + domain):
                vIdx = net["throat.conns"][throat][0]
                for key in poreKeys:
                    if "diameter" in key:
                        net[key] = np.append(
                            net[key], net[key][vIdx]
                        )  # use same diameter for both pores
                        # TODO is this the right approach? # pylint: disable=fixme
                    if "volume" in key:
                        net[key] = np.append(
                            net[key], 0.5 * net[key][vIdx]
                        )  # divide volume on both pores
                        net[key][vIdx] *= 0.5

                poreVolume = net["pore.region_volume"][vIdx]
                eqRadius = (poreVolume * 3.0 / 4.0 * 1.0 / np.pi) ** (1.0 / 3.0)  # assume sphere
                eqArea = np.pi * eqRadius**2

                poreCenter0 = net["pore.centroid"][vIdx]
                poreCenter1 = poreCenter0 + np.array([1e-10, 1e-10, 1e-10])
                throatCenter = 0.5 * (poreCenter0 + poreCenter1)

                for key in throatKeys:
                    if "diameter" in key:
                        net[key] = np.append(net[key], 2 * eqRadius)  # use equivalent diameter
                    if "area" in key:
                        net[key] = np.append(net[key], eqArea)  # use equivalent area
                    if "length" in key:
                        net[key] = np.append(net[key], eqRadius)  # use equivalent radius
                    if "perimeter" in key:
                        net[key] = np.append(net[key], eqRadius)  # use equivalent radius
                    if "centroid" in key:
                        net[key] = np.append(net[key], throatCenter.reshape(-1, 3), axis=0)

            print(
                "Num. void throats after dual grid sanitation",
                np.sum(net["throat.void_void"]),
            )

    if not splitSinglePores:
        op.io.project_to_vtk(project=net.project, filename="before_delete")
        singlePores = np.array(singlePores)
        throatsToDelete = []

        # for eIdx, throat in enumerate(net.throats('interconnect')):
        for eIdx, throat in enumerate(net["throat.conns"]):
            if net["throat.solid_void"][eIdx] == 1 or net["throat.void_solid"][eIdx] == 1:
                vIdx0 = throat[0]
                vIdx1 = throat[1]
                if vIdx0 in singlePores or vIdx1 in singlePores:
                    throatsToDelete.append(eIdx)

        trim(network=net, throats=throatsToDelete, pores=singlePores)
        op.io.project_to_vtk(project=net.project, filename="after_delete")

    return net


def _computeVertexElementData(net, policy):
    """Compute element and vertex data for DGF output"""
    # pylint: disable=too-many-locals

    vertices = net["pore.coords"]
    elements = net["throat.conns"]

    poreInscribedRadius = policy["pore.inscribedRadius"](net)
    poreExtendedRadius = policy["pore.extendedRadius"](net)
    poreVolume = policy["pore.volume"](net)
    poreLabel = policy["pore.label"](net)

    throatRadius = policy["throat.radius"](net)
    throatLength = policy["throat.length"](net)
    throatArea = policy["throat.area"](net)
    throatShapeFactor = policy["throat.throatShapeFactor"](net)
    throatLabel = policy["throat.label"](net)
    throatCenter = net["throat.global_peak"]

    # create dgf data
    isDualNetwork = "pore.solid" in net.keys()
    if isDualNetwork:
        # 1D array (one entry for each throat)
        # containing 0 for void, 1 for solid and 2 for interconnect
        throatDomainType = np.zeros(len(elements), dtype=np.int64)
        throatDomainType[net["throat.solid_solid"]] = 1
        throatDomainType[net["throat.void_solid"]] = 2
        throatDomainType[net["throat.solid_void"]] = 2

        # 1D array (one entry for each pore)
        # containing 0 for void, 1 for solid and 2
        poreDomainType = np.zeros(len(vertices), dtype=np.int64)
        poreDomainType[net["pore.solid"]] = 1
        vertexData = np.stack(
            [
                vertices[:, 0],
                vertices[:, 1],
                vertices[:, 2],
                poreInscribedRadius,
                poreExtendedRadius,
                poreVolume,
                poreLabel,
                poreDomainType,
            ],
            axis=1,
        ).reshape(len(poreLabel), 8)
        elementData = np.stack(
            [
                elements[:, 0],
                elements[:, 1],
                throatRadius,
                throatLength,
                throatArea,
                throatShapeFactor,
                throatCenter[:, 0],
                throatCenter[:, 1],
                throatCenter[:, 2],
                throatLabel,
                throatDomainType,
            ],
            axis=1,
        ).reshape(len(throatLabel), 11)
    else:
        vertexData = np.stack(
            [
                vertices[:, 0],
                vertices[:, 1],
                vertices[:, 2],
                poreInscribedRadius,
                poreExtendedRadius,
                poreVolume,
                poreLabel,
            ],
            axis=1,
        ).reshape(len(poreLabel), 7)
        elementData = np.stack(
            [
                elements[:, 0],
                elements[:, 1],
                throatRadius,
                throatLength,
                throatArea,
                throatShapeFactor,
                throatCenter[:, 0],
                throatCenter[:, 1],
                throatCenter[:, 2],
                throatLabel,
            ],
            axis=1,
        ).reshape(len(throatLabel), 10)

    return vertexData, elementData


def writeDGF(filename, net, policy):
    """DOC ME!"""

    vertexData, elementData = _computeVertexElementData(net, policy)
    isDualNetwork = "pore.solid" in net.keys()

    # write DGF file
    with open(filename, "w") as outputfile:
        outputfile.write("DGF\n")
        outputfile.write(
            "% Vertex parameters: PoreInscribedRadius PoreExtendedRadius PoreVolume PoreLabel"
            f"{' PoreDomainType' if isDualNetwork else ''}\n"
        )
        outputfile.write(
            "% Element parameters: ThroatInscribedRadius ThroatLength ThroatCrossSectionalArea "
            "ThroatShapeFactor ThroatCenterX ThroatCenterY ThroatCenterZ ThroatLabel"
            f"{' ThroatDomainType' if isDualNetwork else ''}\n"
        )
        outputfile.write("Vertex\n")
        outputfile.write("parameters " + str(len(vertexData[0]) - 3) + "\n")  # vertex data
        for i in range(len(vertexData)):
            outputfile.write(" ".join([str(v) for v in vertexData[i]]) + "\n")
        outputfile.write("\n#\n")
        outputfile.write("SIMPLEX\n")
        outputfile.write("parameters " + str(len(elementData[0]) - 2) + "\n")  # cell data
        for i in range(len(elementData)):
            outputfile.write(
                " ".join([str(v if idx > 1 else int(v)) for idx, v in enumerate(elementData[i])])
            )
            outputfile.write("\n")
        outputfile.write("\n\n#\n")
        outputfile.write("BOUNDARYDOMAIN\ndefault 1\n#")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Convert .pnm to .dgf")
    parser.add_argument("file", type=str, help="The .pnm file")
    parser.add_argument("-n", "--outputname", type=str, help="The output file name", default="")
    parser.add_argument(
        "-d",
        "--deleteBoundaryThroats",
        action="store_true",
        help="Delete the boundary throats",
    )
    parser.add_argument(
        "-s",
        "--splitSinglePores",
        action="store_true",
        help="Split pores only connected to the other phase. Will be deleted otherwise.",
    )
    args = vars(parser.parse_args())

    # net contains topology and some labels, pore sizes etc
    network = load(args["file"])

    # maybe remove boundary throats
    if args["deleteBoundaryThroats"]:
        print("Removing boundary throats")
        removeBoundaryThroats(network)

    # check if we have a dual network
    if "pore.solid" in network.keys():
        network = sanitizeDualNetwork(network, args["splitSinglePores"])

    # begin policy //////////////////////////////////////////
    # TODO: How can we conveniently load / exchange this? # pylint: disable=fixme
    def throatLengthFunc(net):
        """DOC ME!"""
        elements = net["throat.conns"]
        length = copy.deepcopy(net["throat.direct_length"])
        poreRadius = net["pore.extended_diameter"] / 2.0

        for eIdx, element in enumerate(elements):
            vIdx0, vIdx1 = element[0], element[1]
            length[eIdx] = length[eIdx] - poreRadius[vIdx0] - poreRadius[vIdx1]

        return np.clip(length, 1e-9, np.max(length))

    def throatShapeFactorFunc(net):  # pylint: disable=unused-argument
        """DOC ME!"""
        transmissibilityG = net["throat.cross_sectional_area"] / net["throat.perimeter"] ** 2
        maxG = 1.0 / 16.0  # for circle
        return np.clip(transmissibilityG, 0.0, maxG)

    def poreLabelFunc(net):  # pylint: disable=unused-argument
        """DOC ME!"""
        label = np.full(len(net["pore.coords"]), -1, dtype=int)
        label[net["pore.xmin"]] = 1
        label[net["pore.xmax"]] = 2
        label[net["pore.ymin"]] = 3
        label[net["pore.ymax"]] = 4
        try:
            label[net["pore.zmin"]] = 5
            label[net["pore.zmax"]] = 6
        except (IndexError, KeyError):
            pass

        if "pore.fake_void" in net.keys() and "pore.fake_solid" in net.keys():
            label[net["pore.fake_void"]] = 22
            label[net["pore.fake_solid"]] = 33

        return label

    def throatLabelFunc(net, policy=None):
        """DOC ME!"""
        # TODO add policy # pylint: disable=fixme
        _ = policy
        pLabels = poreLabelFunc(net)
        label = np.full(len(net["throat.conns"]), -1, dtype=int)

        for eIdx, element in enumerate(net["throat.conns"]):
            label[eIdx] = np.max(
                pLabels[element]
            )  # take the label with higher value, -> add policy here

        if "throat.fake_void" in net.keys() and "throat.fake_solid" in net.keys():
            label[net["throat.fake_void"]] = 22
            label[net["throat.fake_solid"]] = 33

        return label

    POLICY = {
        "pore.inscribedRadius": lambda n: n["pore.inscribed_diameter"] / 2.0
        if "pore.inscribed_diameter" in n
        else n["pore.extended_diameter"] / 2.0,
        "pore.extendedRadius": lambda n: n["pore.extended_diameter"] / 2.0,
        "pore.volume": lambda n: n["pore.region_volume"],
        "pore.label": poreLabelFunc,
        "throat.radius": lambda n: n["throat.inscribed_diameter"] / 2.0,
        "throat.length": throatLengthFunc,
        "throat.area": lambda n: n["throat.cross_sectional_area"],
        "throat.throatShapeFactor": throatShapeFactorFunc,
        "throat.label": throatLabelFunc,
    }
    # end policy /////////////////////////////////////////////

    if not args["outputname"]:
        outputFilename = os.path.splitext(args["file"])[0] + ".dgf"
    else:
        outputFilename = os.path.splitext(args["outputname"])[0] + ".dgf"

    print("Writing", outputFilename)
    writeDGF(outputFilename, network, POLICY)
    # export to vtk
    op.io.project_to_vtk(project=network.project, filename=outputFilename)
    print("Writing", outputFilename + ".vtp")
