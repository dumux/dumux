import os
import openpnm as op
import numpy as np
import argparse
from openpnm.topotools import trim
import copy


def load(filename):
    # create workspace instance
    ws = op.Workspace()
    ws.clear()
    project = ws.load_project(filename=args["file"])
    print(project)
    return project["net_01"], project["geo_01"]


def removeBoundaryThroats(net):
    pores = net.pores("boundary")
    net["throat.boundary"] = False

    for idx, throat in enumerate (net["throat.conns"]):
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
        except:
            print(f"{label} not found. Skipping.")

    trim(network=net, throats=throats)
    trim(network=net, pores=net["pore.delete"])


def poreLabel(net):
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


def throatLabel(net, policy=None):
    # TODO add policy
    pLabels = poreLabel(net)
    label = np.full(len(net["throat.conns"]), -1, dtype=int)

    for eIdx, element in enumerate(net["throat.conns"]):
        label[eIdx] = np.max(
            pLabels[element]
        )  # take the label with higher value, TODO add policy here

    if "throat.fake_void" in net.keys() and "throat.fake_solid" in net.keys():
        label[net["throat.fake_void"]] = 22
        label[net["throat.fake_solid"]] = 33

    return label


def getDualGridCoordinationNumber(net):
    numPoresTotal = len(net["pore.coords"])
    coordNumVoid = np.zeros(numPoresTotal, dtype=np.int64)
    coordNumSolid = np.zeros(numPoresTotal, dtype=np.int64)

    for vT in net["throat.conns"][net["throat.void"]]:
        coordNumVoid[vT[0]] += 1
        coordNumVoid[vT[1]] += 1

    for sT in net["throat.conns"][net["throat.solid"]]:
        coordNumSolid[sT[0]] += 1
        coordNumSolid[sT[1]] += 1

    return coordNumVoid, coordNumSolid


def sanitizeDualNetwork(net, splitSinglePores):
    geo = net.project["geo_01"]

    print("Num. void throats before dual grid sanitation", np.sum(net["throat.void"]))
    print(net)

    # calculate coordination number for void-void and solid-solid connections
    coordNumVoid, coordNumSolid = getDualGridCoordinationNumber(net)

    poreKeys = [key for key in geo if key.startswith("pore.")]
    throatKeys = [key for key in geo if key.startswith("throat.")]

    # HACK We want to add geometry information for the new fake pores and throats but it is not clear how to change the entries in geo.
    # We therefore replace the geo object by a dictionary which behaves similarly.
    # This should be improved!
    if splitSinglePores:
        newGeo = {}
        for key in geo:
            newGeo[key] = geo[key]
        geo = newGeo

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

        # identify pores only connect to pores of other phase and add an artificial throat or delete those pores
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

            # add some articial geometry data
            for throat in net.throats("fake_" + domain):
                vIdx = net["throat.conns"][throat][0]
                for key in poreKeys:
                    if "diameter" in key:
                        geo[key] = np.append(
                            geo[key], geo[key][vIdx]
                        )  # use same diamter for both pores # TODO is this the right approach?
                    if "volume" in key:
                        geo[key] = np.append(
                            geo[key], 0.5 * geo[key][vIdx]
                        )  # divide volume on both pores
                        geo[key][vIdx] *= 0.5

                poreVolume = geo["pore.region_volume"][vIdx]
                eqRadius = (poreVolume * 3.0 / 4.0 * 1.0 / np.pi) ** (
                    1.0 / 3.0
                )  # assume sphere
                eqArea = np.pi * eqRadius ** 2

                poreCenter0 = geo["pore.centroid"][vIdx]
                poreCenter1 = poreCenter0 + np.array([1e-10, 1e-10, 1e-10])
                throatCenter = 0.5 * (poreCenter0 + poreCenter1)

                for key in throatKeys:
                    if "diameter" in key:
                        geo[key] = np.append(
                            geo[key], 2 * eqRadius
                        )  # use equivalent diameter
                    if "area" in key:
                        geo[key] = np.append(geo[key], eqArea)  # use equivalent area
                    if "length" in key:
                        geo[key] = np.append(
                            geo[key], eqRadius
                        )  # use equivalent radius
                    if "perimeter" in key:
                        geo[key] = np.append(
                            geo[key], eqRadius
                        )  # use equivalent radius
                    if "centroid" in key:
                        geo[key] = np.append(
                            geo[key], throatCenter.reshape(-1, 3), axis=0
                        )

            print(
                "Num. void throats after dual grid sanitation",
                np.sum(net["throat.void"]),
            )

    if not splitSinglePores:
        op.io.VTK.save(network=net, filename="before_delete")
        singlePores = np.array(singlePores)
        throatsToDelete = []

        # for eIdx, throat in enumerate(net.throats('interconnect')):
        for eIdx, throat in enumerate(net["throat.conns"]):
            if net["throat.interconnect"][eIdx] == 1:
                vIdx0 = throat[0]
                vIdx1 = throat[1]
                if vIdx0 in singlePores or vIdx1 in singlePores:
                    throatsToDelete.append(eIdx)

        trim(network=net, throats=throatsToDelete, pores=singlePores)
        op.io.VTK.save(network=net, filename="after_delete")

    return net, geo


def writeDGF(filename, net, geo, policy):
    isDualNetwork = "pore.solid" in net.keys()

    vertices = net["pore.coords"]
    elements = net["throat.conns"]

    poreRadius = policy["pore.radius"](net, geo)
    poreVolume = policy["pore.volume"](net, geo)
    poreLabel = policy["pore.label"](net, geo)

    throatRadius = policy["throat.radius"](net, geo)
    throatLength = policy["throat.length"](net, geo)
    throatArea = policy["throat.area"](net, geo)
    throatShapeFactor = policy["throat.throatShapeFactor"](net, geo)
    throatLabel = policy["throat.label"](net, geo)
    throatCenter = net["throat.global_peak"]

    # create dgf data
    if isDualNetwork:
        # 1D array (one entry for each throat) containg 0 for void, 1 for solid and 2 for interconnect
        throatDomainType = np.zeros(len(elements), dtype=np.int64)
        throatDomainType[net["throat.solid"]] = 1
        throatDomainType[net["throat.interconnect"]] = 2

        # 1D array (one entry for each pore) containg 0 for void, 1 for solid and 2
        poreDomainType = np.zeros(len(vertices), dtype=np.int64)
        poreDomainType[net["pore.solid"]] = 1
        vdata = np.stack(
            [
                vertices[:, 0],
                vertices[:, 1],
                vertices[:, 2],
                poreRadius,
                poreVolume,
                poreLabel,
                poreDomainType,
            ],
            axis=1,
        ).reshape(len(poreLabel), 7)
        edata = np.stack(
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
        vdata = np.stack(
            [
                vertices[:, 0],
                vertices[:, 1],
                vertices[:, 2],
                poreRadius,
                poreVolume,
                poreLabel,
            ],
            axis=1,
        ).reshape(len(poreLabel), 6)
        edata = np.stack(
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

    # write DGF file
    with open(filename, "w") as outputfile:
        outputfile.write("DGF\n")
        outputfile.write(
            "% Vertex parameters: PoreInscribedRadius PoreVolume PoreLabel"
            + (" PoreDomainType" if isDualNetwork else "")
            + "\n"
        )
        outputfile.write(
            "% Element parameters: ThroatInscribedRadius ThroatLength ThroatCrossSectionalArea ThroatShapeFactor ThroatCenterX ThroatCenterY ThroatCenterZ ThroatLabel"
            + (" ThroatDomainType" if isDualNetwork else "")
            + "\n"
        )
        outputfile.write("Vertex\n")
        outputfile.write("parameters " + str(len(vdata[0]) - 3) + "\n")  # vertex data
        for i in range(len(poreLabel)):
            outputfile.write(" ".join([str(v) for v in vdata[i]]) + "\n")
        outputfile.write("\n#\n")
        outputfile.write("SIMPLEX\n")
        outputfile.write("parameters " + str(len(edata[0]) - 2) + "\n")  # cell data
        for i in range(len(throatLabel)):
            outputfile.write(
                " ".join(
                    [str(v if idx > 1 else int(v)) for idx, v in enumerate(edata[i])]
                )
                + "\n"
            )
        outputfile.write("\n#\n")
        outputfile.write("BOUNDARYDOMAIN\ndefault 1\n#")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Convert .pnm to .dgf")
    parser.add_argument("file", type=str, help="The .pnm file")
    parser.add_argument(
        "-n", "--outputname", type=str, help="The output file name", default=""
    )
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

    # net contains topology and some labels, geo contains pore sizes etc
    net, geo = load(args["file"])
    isDualNetwork = "pore.solid" in net.keys()

    # maybe remove boundary throats
    if args["deleteBoundaryThroats"]:
        print("Removing boundary throats")
        removeBoundaryThroats(net)

    if isDualNetwork:
        net, geo = sanitizeDualNetwork(net, args["splitSinglePores"])

    # begin policy //////////////////////////////////////////
    # TODO: how can we conveniently load / exchange this?
    def throatLength(n, g):
        elements = n["throat.conns"]
        length = copy.deepcopy(g["throat.direct_length"])
        poreRadius = g["pore.extended_diameter"] / 2.0

        for eIdx, element in enumerate(elements):
            vIdx0, vIdx1 = element[0], element[1]
            length[eIdx] = length[eIdx] - poreRadius[vIdx0] - poreRadius[vIdx1]

        return np.clip(length, 1e-9, np.max(length))

    def throatShapeFactor(n, g):
        G = g["throat.cross_sectional_area"] / g["throat.perimeter"] ** 2
        maxG = 1.0 / 16.0  # for circle
        return np.clip(G, 0.0, maxG)

    policy = {
        "pore.radius": lambda n, g: g["pore.extended_diameter"] / 2.0,
        "pore.volume": lambda n, g: g["pore.region_volume"],
        "pore.label": lambda n, g: poreLabel(n),
        "throat.radius": lambda n, g: g["throat.inscribed_diameter"] / 2.0,
        "throat.length": throatLength,
        "throat.area": lambda n, g: g["throat.cross_sectional_area"],
        "throat.throatShapeFactor": throatShapeFactor,
        "throat.label": throatLabel,
    }
    # end policy /////////////////////////////////////////////

    if not args["outputname"]:
        filename = os.path.splitext(args["file"])[0] + ".dgf"
    else:
        filename = os.path.splitext(args["outputname"])[0] + ".dgf"

    print("Writing", filename)
    writeDGF(filename, net, geo, policy)
    # export to vtk
    op.io.VTK.save(network=net, filename=filename)
    print("Writing", filename + ".vtp")
