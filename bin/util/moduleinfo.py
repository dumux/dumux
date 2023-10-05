# SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

"""
Read information from dune.module files
"""

import os
from util.common import runCommand
from util.common import callFromPath


def extractModuleInfos(moduleFile, keys):
    """Extract information about a Dune module from its dune.module file"""
    results = {}
    with open(moduleFile, "r") as modFile:
        for line in modFile.readlines():
            line = line.strip("\n").split(":")
            if line[0] in keys:
                results[line[0]] = line[1].strip()
            if len(results) == len(keys):
                break

    if len(results) != len(keys):
        errMsg = "Could not extract requested information for all keys.\n"
        errMsg += "Requested keys: " + ", ".join(keys) + "\n"
        errMsg += "Processed keys: " + ", ".join(list(results))
        raise RuntimeError(errMsg)

    return results


def getModuleFile(modulePath):
    """Read the dune.module file"""
    modFile = os.path.join(modulePath, "dune.module")
    if not os.path.exists(modFile):
        raise RuntimeError("Could not find module file")
    return modFile


def getModuleInfo(modulePath, key):
    """Read information about Dune module"""
    return extractModuleInfos(getModuleFile(modulePath), [key])[key]


def parseModuleList(dunecontrolOutput):
    """Determine the module dependencies from dunecontrol terminal output"""
    for line in dunecontrolOutput.split("\n"):
        if "going to build" in line:
            line = line.replace("going to build", "").strip("-").strip("\n").strip().split(" ")
            return line
    return []


def getDependencies(modulePath, verbose=False, includeSelf=False):
    """Get the dependencies of a Dune module"""
    modName = getModuleInfo(modulePath, "Module")
    parentPath = os.path.join(modulePath, "../")
    duneControlPath = os.path.join(parentPath, "dune-common/bin/dunecontrol")
    if not os.path.exists(duneControlPath):
        raise RuntimeError(f"Could not find dunecontrol, expected it to be in {duneControlPath}")

    dcOutput = callFromPath(parentPath)(runCommand)(
        f"./dune-common/bin/dunecontrol --module={modName}"
    )

    if not dcOutput:
        raise RuntimeError("Error: call to dunecontrol failed.")

    dependencyList = parseModuleList(dcOutput)

    if not includeSelf:
        dependencyList.remove(modName)

    if verbose:
        print(" -- Determined the following dependencies: " + ", ".join(dependencyList))
        print(" -- Searching the respective directories...")

    result = []
    parentFiles = [os.path.join(parentPath, d) for d in os.listdir(parentPath)]
    for path in filter(os.path.isdir, parentFiles):
        try:
            depModName = getModuleInfo(path, "Module")
        except RuntimeError:
            if verbose:
                print(
                    f" --- skipping folder '{path}' " "as it could not be identified as dune module"
                )
        else:
            if verbose:
                print(f" --- visited module '{depModName}'")
            if depModName in dependencyList:
                result.append({"name": depModName, "folder": os.path.basename(path), "path": path})

    if len(result) != len(dependencyList):
        raise RuntimeError("Could not find the folders of all dependencies")
    if verbose:
        print(" -- Found all module folders of the dependencies.")
    return result
