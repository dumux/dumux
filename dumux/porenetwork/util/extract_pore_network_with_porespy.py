#!/usr/bin/env python3

"""
A simple script for conveniently running the SNOW algorithm of PoreSpy
(watershed image segmentation) in order to extract a (dual) pore network
from a 2D or 3D raster image. Accepts .raw, .png and .nrrd images.
Requires the PoreSpy (https://porespy.org/) and OpenPNM package (https://openpnm.org/).
"""

from argparse import ArgumentParser
import numpy as np
import openpnm as op
import porespy as ps
from porespy.tools import randomize_colors
from loguru import logger


def printReadMessage(fileName, fileType, numVoxels):
    dim = len(numVoxels)
    print(
        "Reading",
        str(dim) + "D",
        fileType,
        "image file",
        fileName,
        "with",
        numVoxels,
        "voxels" if dim == 3 else "pixels",
    )


def openRawFile(fileName, numVoxels):
    assert len(numVoxels) == 2 or len(numVoxels) == 3
    im = np.fromfile(fileName, dtype="uint8", sep="").reshape(numVoxels)
    im = np.array(im, dtype=bool)
    dim = len(im.shape)
    printReadMessage(fileName, "raw", im.shape)
    if dim == 3:
        im = np.swapaxes(im, 0, 2)
    return im


def openNrrdFile(fileName):
    import nrrd

    im, header = nrrd.read(fileName)
    im = np.array(im, dtype=bool)
    printReadMessage(fileName, "nrrd", im.shape)
    return im


def openPngFile(fileName):
    import matplotlib.pyplot as plt

    im = plt.imread(fileName)
    im = np.array(im, dtype=bool)
    printReadMessage(fileName, "png", im.shape)
    return im


def openFile(fileName, numVoxels):
    if fileName.endswith(".raw"):
        return openRawFile(fileName, numVoxels)
    elif fileName.endswith(".nrrd"):
        if len(numVoxels) > 0:
            logger.warning("Argument numVoxels ignored for nrrd files")
        return openNrrdFile(fileName)
    elif fileName.endswith(".png"):
        if len(numVoxels) > 0:
            logger.warning("Argument numVoxels ignored for png files")
        return openPngFile(fileName)
    else:
        raise NotImplementedError(fileName + "not supported")


def extrude2D(im, layers, topAndBottomValues):
    if not len(im.shape) == 2:
        raise IOError("Extrusion only works for 2D images")
    if not len(topAndBottomValues) == 2:
        raise IOError("topAndBottomValues can only have two entries")
    im = np.tile(np.atleast_3d(im), reps=[1, 1, layers])
    im[..., 0:1] = topAndBottomValues[0]
    im[..., -1:] = topAndBottomValues[1]
    print("Extruding image to 3D.")
    return im


def runSnow(im, resolution, sigma, rMax, boundaryWidth, useMarchingCubes):
    return ps.networks.snow2(
        im,
        voxel_size=resolution,
        sigma=sigma,
        r_max=rMax,
        boundary_width=boundaryWidth,
        accuracy=("high" if useMarchingCubes else "standard"),
    )


def runDualSnow(im, resolution, sigma, rMax, boundaryWidth, useMarchingCubes):
    return ps.networks.snow2(
        im + 1,
        voxel_size=resolution,
        sigma=sigma,
        r_max=rMax,
        boundary_width=boundaryWidth,
        accuracy=("high" if useMarchingCubes else "standard"),
        phase_alias={1: "solid", 2: "void"},
    )


def writeSegmentationInfoFile(snowOutput, fileName, resolution):
    dim = len(snowOutput.phases.shape)
    if dim == 2:
        ps.io.to_vtk(
            np.array(randomize_colors(snowOutput.regions), dtype=int)[:, :, np.newaxis],
            fileName,
            voxel_size=resolution,
        )
    else:
        ps.io.to_vtk(
            np.array(randomize_colors(snowOutput.regions), dtype=int)[:, :, :],
            fileName,
            voxel_size=resolution,
        )


def extractNetwork(
    image,
    resolution,
    outputName="default_output",
    sigma=0.4,
    rMax=4,
    boundaryWidth=3,
    dualSnow=False,
    useMarchingCubes=False,
):
    """
    Extracts a porenetwork from an image file and exports
    it as VTK Polydata file (*.vtp) and OpenPNM network file (*.pnm).

    Note: Use the script `openpnm2dgf.py` to convert the OpenPNM network file format
    to the Dune grid file format (*.dgf) that can be read by DuMux.
    """

    print("Porosity is", ps.metrics.porosity(im))

    if not dualSnow:
        snowOutput = runSnow(
            im=im,
            resolution=resolution,
            sigma=sigma,
            rMax=rMax,
            boundaryWidth=boundaryWidth,
            useMarchingCubes=useMarchingCubes,
        )
    else:
        snowOutput = runDualSnow(
            im=im,
            resolution=resolution,
            sigma=sigma,
            rMax=rMax,
            boundaryWidth=boundaryWidth,
            useMarchingCubes=useMarchingCubes,
        )

    # write output for segmentation
    writeSegmentationInfoFile(snowOutput, outputName + "_segmentation", resolution)

    # use PoreSpy to sanitize some parameters (might become obsolete as some point)
    pn, geo = op.io.PoreSpy.import_data(snowOutput.network)

    # trimming pore network to avoid singularity
    print("Number of pores before trimming: ", pn.Np)
    h = pn.check_network_health()
    op.topotools.trim(network=pn, pores=h["trim_pores"])
    print("Number of pores after trimming: ", pn.Np)

    filename = outputName + ("_dual" if dualSnow else "")

    # export to vtk
    op.io.VTK.save(network=pn, filename=filename)
    print("Writing", filename + ".vtp")

    # export network for OpenPNM
    pn.project.save_project(filename)
    print("Writing", filename + ".pnm")


if __name__ == "__main__":
    parser = ArgumentParser(description="extract network from image using PoreSpy")
    parser.add_argument("file", type=str, help="The image file")
    parser.add_argument(
        "-r", "--resolution", type=float, help="The pixel/voxel size", required=True
    )
    parser.add_argument(
        "-e",
        "--extrusionLayers",
        type=int,
        help="Extrude a 2D image to 3D. Specify the number of layers",
        default=0,
    )
    parser.add_argument(
        "-sigma",
        type=float,
        help="The standard deviation of the Gaussian filter",
        default=0.4,
    )
    parser.add_argument("-rMax", type=float, help="Radius for finding peaks", default=4)
    parser.add_argument(
        "-tb",
        "--topAndBottomValues",
        nargs="+",
        type=int,
        help="Define the image value of the first and last vertical row in case of extrusion",
        default=[0, 0],
    )
    parser.add_argument(
        "-v",
        "--numVoxels",
        nargs="+",
        type=int,
        help="The number of pixels/voxels in each direction",
        default=[],
    )
    parser.add_argument(
        "-bw",
        "--boundaryWidth",
        nargs="+",
        type=int,
        help="The number of pixels/voxels in each direction",
        default=[],
    )
    parser.add_argument(
        "-n", "--outputName", type=str, help="The output file name", default=""
    )
    parser.add_argument(
        "-m",
        "--marchingCubesArea",
        action="store_true",
        help="Use marching cubes alg. to calculate throat area",
    )
    parser.add_argument(
        "-d", "--dualSnow", action="store_true", help="Extract dual network"
    )
    parser.add_argument(
        "-i", "--invert", action="store_true", help="Invert the image data"
    )
    args = vars(parser.parse_args())

    im = openFile(args["file"], args["numVoxels"])
    dim = len(im.shape)

    if not args["outputName"]:
        # poreSpy has issues with dots in file names. Replace with underscores here and strip file suffix.
        args["outputName"] = ("_").join((args["file"]).split(".")[0:-1])

    if len(args["boundaryWidth"]) == 0:
        args["boundaryWidth"] = 3

    if args["extrusionLayers"] > 0:
        im = extrude2D(im, args["extrusionLayers"], args["topAndBottomValues"])
        dim = 3

    if args["invert"]:
        print("Inverting data")
        im = ~im

    extractNetwork(
        im,
        args["resolution"],
        outputName=args["outputName"],
        sigma=args["sigma"],
        rMax=args["rMax"],
        boundaryWidth=args["boundaryWidth"],
        dualSnow=args["dualSnow"],
        useMarchingCubes=args["marchingCubesArea"],
    )
