import os
from util.common import runCommand
from util.common import callFromPath


def extractModuleInfos(moduleFile, keys):
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
        errMsg += "Processed keys: " + ", ".join([k for k in results])
        raise RuntimeError(errMsg)

    return results


def getModuleFile(modulePath):
    modFile = os.path.join(modulePath, "dune.module")
    if not os.path.exists(modFile):
        raise RuntimeError("Could not find module file")
    return modFile


def getModuleInfo(modulePath, key):
    return extractModuleInfos(getModuleFile(modulePath), [key])[key]


def getDependencies(modulePath, verbose=False, includeSelf=False):
    modName = getModuleInfo(modulePath, "Module")
    parentPath = os.path.join(modulePath, "../")
    duneControlPath = os.path.join(parentPath, "dune-common/bin/dunecontrol")
    if not os.path.exists(duneControlPath):
        raise RuntimeError(
            "Could not find dunecontrol, expected it to be in {}".format(duneControlPath)
        )

    dcOutput = callFromPath(parentPath)(runCommand)(
        "./dune-common/bin/dunecontrol --module={}".format(modName)
    )

    if not dcOutput:
        raise RuntimeError("Error: call to dunecontrol failed.")

    for line in dcOutput.split("\n"):
        if "going to build" in line:
            line = line.replace("going to build", "").strip("-")
            line = line.strip("\n").strip()
            line = line.split(" ")
            deps = line

    if not includeSelf:
        deps.remove(modName)

    if verbose:
        print(" -- Determined the following dependencies: " + ", ".join(deps))
        print(" -- Searching the respective directories...")

    result = []
    parentFiles = [os.path.join(parentPath, d) for d in os.listdir(parentPath)]
    for path in filter(os.path.isdir, parentFiles):
        try:
            depModName = getModuleInfo(path, "Module")
        except Exception:
            if verbose:
                print(
                    f" --- skipping folder '{path}' " "as it could not be identifed as dune module"
                )
        else:
            if verbose:
                print(" --- visited module '{}'".format(depModName))
            if depModName in deps:
                result.append({"name": depModName, "folder": os.path.basename(path), "path": path})

    if len(result) != len(deps):
        raise RuntimeError("Could not find the folders of all dependencies")
    elif verbose:
        print(" -- Found all module folders of the dependencies.")
    return result
