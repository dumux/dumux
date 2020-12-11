import os
from common import runCommand
from common import callFromPath

# extract information (for the given keys) from a module file
def extractModuleInfos(moduleFile, keys):
    results = {}
    with open(moduleFile, 'r') as modFile:
        for line in modFile.readlines():
            line = line.strip('\n').split(':')
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

# get info file of a module
def getModuleFile(modulePath):
    modFile = os.path.join(modulePath, 'dune.module')
    if not os.path.exists(modFile):
        raise RuntimeError("Could not find module file")
    return modFile

# retrieve single information from a module
def getModuleInfo(modulePath, key):
    return extractModuleInfos(getModuleFile(modulePath), [key])[key]

# get the dependencies of a dune module located in the given directory
def getDependencies(modulePath, verbose=False):
    modName = getModuleInfo(modulePath, 'Module')
    parentPath = os.path.join(modulePath, '../')
    duneControlPath = os.path.join(parentPath, 'dune-common/bin/dunecontrol')
    if not os.path.exists(duneControlPath):
        raise RuntimeError('Could not find dunecontrol, expected it to be in {}'.format(duneControlPath))

    run = callFromPath(parentPath)(runCommand)
    for line in run('./dune-common/bin/dunecontrol --module={}'.format(modName)).split('\n'):
        if "going to build" in line:
            line = line.replace('going to build', '').replace('---', '').replace('done', '')
            line = line.strip('\n').strip()
            line = line.split(' ')
            deps = line

    # Now we look for the folders with the modules
    if verbose:
        print("Determined the following dependencies: " + ", ".join(deps))

    result = []
    for dir in [d for d in os.listdir(parentPath) if os.path.isdir(os.path.join(parentPath, d))]:
        try: depModName = getModuleInfo(os.path.join(parentPath, dir), 'Module')
        except: print('--- Note: skipping folder "' + dir + '" as it could not be identifed as dune module')
        else:
            if depModName in deps:
                result.append({'name': depModName, 'folder': dir})

    if len(result) != len(deps):
        raise RuntimeError("Could not find the folders of all dependencies")
    elif verbose:
        print("Found all module folders of the dependencies.")

    return result
