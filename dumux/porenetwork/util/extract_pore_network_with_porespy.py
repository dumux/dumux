#!/usr/bin/env python3
# SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later


"""
A simple script for conveniently running the SNOW algorithm of PoreSpy
(watershed image segmentation) in order to extract a (dual) pore network
from a 2D or 3D raster image. Accepts .raw, .png and .nrrd images.
Requires the PoreSpy (https://porespy.org/) and OpenPNM package (https://openpnm.org/).
"""

from argparse import ArgumentParser
import numpy as np
import openpnm as op  # pylint: disable=import-error
import porespy as ps  # pylint: disable=import-error
from porespy.tools import randomize_colors  # pylint: disable=import-error
from loguru import logger  # pylint: disable=import-error


def printReadMessage(fileName, fileType, numVoxels):
    """DOC ME!"""

    dimension = len(numVoxels)
    print(
        "Reading",
        str(dimension) + "D",
        fileType,
        "image file",
        fileName,
        "with",
        numVoxels,
        "voxels" if dimension == 3 else "pixels",
    )


def openRawFile(fileName, numVoxels):
    """DOC ME!"""

    assert len(numVoxels) == 2 or len(numVoxels) == 3
    image = np.fromfile(fileName, dtype="uint8", sep="").reshape(numVoxels)
    image = np.array(image, dtype=bool)
    dimension = len(image.shape)
    printReadMessage(fileName, "raw", image.shape)
    if dimension == 3:
        image = np.swapaxes(image, 0, 2)
    return image


def openNrrdFile(fileName):
    """DOC ME!"""

    import nrrd  # pylint: disable=import-outside-toplevel,import-error

    image, _ = nrrd.read(fileName)
    image = np.array(image, dtype=bool)
    printReadMessage(fileName, "nrrd", image.shape)
    return image


def openPngFile(fileName):
    """DOC ME!"""

    import matplotlib.pyplot as plt  # pylint: disable=import-outside-toplevel,import-error

    image = plt.imread(fileName)
    image = np.array(image, dtype=bool)
    printReadMessage(fileName, "png", image.shape)
    return image


def openFile(fileName, numVoxels):
    """Opens the image with fileName and number of voxels in each direction
    detecting the file type from extension.

    Returns the image as a numpy array.
    """

    if fileName.endswith(".raw"):
        return openRawFile(fileName, numVoxels)
    if fileName.endswith(".nrrd"):
        if len(numVoxels) > 0:
            logger.warning("Argument numVoxels ignored for nrrd files")
        return openNrrdFile(fileName)
    if fileName.endswith(".png"):
        if len(numVoxels) > 0:
            logger.warning("Argument numVoxels ignored for png files")
        return openPngFile(fileName)

    raise NotImplementedError(fileName + "not supported")


def extrude2D(image, layers, topAndBottomValues):
    """DOC ME!"""

    if not len(image.shape) == 2:
        raise IOError("Extrusion only works for 2D images")
    if not len(topAndBottomValues) == 2:
        raise IOError("topAndBottomValues can only have two entries")
    image = np.tile(
        np.atleast_3d(image), reps=[1, 1, layers]
    )  # reshape in 3D array and repeat layers times
    image[..., 0:1] = topAndBottomValues[0]
    image[..., -1:] = topAndBottomValues[1]
    print("Extruding image to 3D.")
    return image


def runSnow(
    image, resolution, sigma, rMax, boundaryWidth, useMarchingCubes
):  # pylint: disable=too-many-arguments
    """DOC ME!"""

    # pylint: disable=no-member
    return ps.networks.snow2(
        image,
        voxel_size=resolution,
        sigma=sigma,
        r_max=rMax,
        boundary_width=boundaryWidth,
        accuracy=("high" if useMarchingCubes else "standard"),
    )


def runDualSnow(
    image, resolution, sigma, rMax, boundaryWidth, useMarchingCubes
):  # pylint: disable=too-many-arguments
    """DOC ME!"""

    # pylint: disable=no-member
    return ps.networks.snow2(
        image + 1,
        voxel_size=resolution,
        sigma=sigma,
        r_max=rMax,
        boundary_width=boundaryWidth,
        accuracy=("high" if useMarchingCubes else "standard"),
        phase_alias={1: "solid", 2: "void"},
    )


def writeSegmentationInfoFile(snowOutput, fileName, resolution):
    """DOC ME!"""

    dimension = len(snowOutput.phases.shape)
    if dimension == 2:
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
):  # pylint: disable=too-many-arguments
    """
    Extracts a porenetwork from an image file and exports
    it as VTK Polydata file (*.vtp) and OpenPNM network file (*.pnm).

    Note: Use the script `openpnm2dgf.py` to convert the OpenPNM network file format
    to the Dune grid file format (*.dgf) that can be read by DuMux.
    """

    print("Porosity is", ps.metrics.porosity(image))

    if not dualSnow:
        snowOutput = runSnow(
            image=image,
            resolution=resolution,
            sigma=sigma,
            rMax=rMax,
            boundaryWidth=boundaryWidth,
            useMarchingCubes=useMarchingCubes,
        )
    else:
        snowOutput = runDualSnow(
            image=image,
            resolution=resolution,
            sigma=sigma,
            rMax=rMax,
            boundaryWidth=boundaryWidth,
            useMarchingCubes=useMarchingCubes,
        )

    # write output for segmentation
    writeSegmentationInfoFile(snowOutput, outputName + "_segmentation", resolution)

    # use PoreSpy to sanitize some parameters (might become obsolete as some point)
    porenetwork = op.io.network_from_porespy(snowOutput.network)  # porenetwork in OpenPNM format

    # trimming pore network to avoid singularity (remove pores in disconnected clusters,
    # this includes isolated pores).
    # isolated pores could be removed separately by using
    # op.topotools.trim(network=porenetwork, pores=health["isolated_pores"])
    print("Number of pores before trimming: ", porenetwork.Np)
    health = op.utils.check_network_health(porenetwork)
    op.topotools.trim(network=porenetwork, pores=health["disconnected_pores"])
    print("Number of pores after trimming: ", porenetwork.Np)

    filename = outputName + ("_dual" if dualSnow else "")

    # export to vtk
    op.io.project_to_vtk(project=porenetwork.project, filename=filename)
    print("Writing", filename + ".vtp")

    # export network for OpenPNM
    op.utils.Workspace().save_project(project=porenetwork.project, filename=filename)
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
    parser.add_argument("-n", "--outputName", type=str, help="The output file name", default="")
    parser.add_argument(
        "-m",
        "--marchingCubesArea",
        action="store_true",
        help="Use marching cubes alg. to calculate throat area",
    )
    parser.add_argument("-d", "--dualSnow", action="store_true", help="Extract dual network")
    parser.add_argument("-i", "--invert", action="store_true", help="Invert the image data")
    args = vars(parser.parse_args())

    imageData = openFile(args["file"], args["numVoxels"])

    if not args["outputName"]:
        # poreSpy has issues with dots in file names
        # Replace with underscores here and strip file suffix.
        args["outputName"] = ("_").join((args["file"]).split(".")[0:-1])

    if len(args["boundaryWidth"]) == 0:
        args["boundaryWidth"] = 3

    if args["extrusionLayers"] > 0:
        imageData = extrude2D(imageData, args["extrusionLayers"], args["topAndBottomValues"])

    if args["invert"]:
        print("Inverting data")
        imageData = ~imageData

    extractNetwork(
        imageData,
        args["resolution"],
        outputName=args["outputName"],
        sigma=args["sigma"],
        rMax=args["rMax"],
        boundaryWidth=args["boundaryWidth"],
        dualSnow=args["dualSnow"],
        useMarchingCubes=args["marchingCubesArea"],
    )
