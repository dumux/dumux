#!/usr/bin/env python3

# script to generate an install script for a dune-module,
# accounting for non-published commits and local changes

import os
import sys
import argparse
import subprocess
from getdependencies import *
from getusedversions import *

if sys.version_info[0] < 3:
    sys.exit("\nERROR: Python3 required")

###################
# parse arguments
parser = argparse.ArgumentParser(
    description='This script generates an install script for your dune module,' \
                'taking into account non-published commits and local changes.\n' \
                'This expects that your module is a git repository and that a ' \
                'remote origin exists and has been set already.'
)
parser.add_argument('-p', '--path', required=True, help='the path to your dune module')
parser.add_argument('-f', '--filename', required=False, help='name of the file in which to write the install script')
cmdArgs = vars(parser.parse_args())

####################
# Welcome message
modPath = os.path.join(os.getcwd(), cmdArgs['path'])
modParentPath = os.path.join(modPath, '../')
modName = getModuleName(modPath)
if modName == None:
    sys.exit("\nERROR: Could not determine module name. Make sure the path to your module\n" \
             "         is correct and that a dune.module file is present.")

print("\n-- Creating install script for the module '{}' in the folder '{}'".format(modName, modPath))

###################
# get dependencies
print("\n-- Determining the dependencies")
deps = getDependencies(modPath)
if len(deps) > 0:
    print("-> Found the following dependencies")
    print("\t| {:^50} | {:^50} |".format('module name', 'module folder'))
    print("\t" + 107*'-')
    for dep in deps: print('\t| {:^50} | {:^50} |'.format(dep['name'], dep['folder']))
else:
    print("-> Module has no dependencies")

depNames = [dep['name'] for dep in deps]
depFolders = [dep['folder'] for dep in deps]

#################################
# determine specific commits of all modules
print("\n-- Determining the module versions")
versions = getUsedVersions(depFolders)
if len(versions) != len(depNames):
    sys.exit("ERROR: Could not determine versions of all modules")

print("-> Detected the following versions")
print("\t| {:^50} | {:^50} | {:^50} | {:^30} |".format('module folder', 'branch', 'commit hash', 'commit date'))
print("\t" + 163*'-')
for folder, versionInfo in versions.items():
    print("\t| {:^50} | {:^50} | {:^50} | {:^30} |".format(folder, versionInfo['branch'], versionInfo['revision'], versionInfo['date']))


#################################
# create patches if necessary
print("\n-- Creating patches for unpublished commits and uncommitted changes")
patches = getPatches(depFolders)

# create patch files
if len(patches) > 0:
    patchesPath = os.path.join(modPath, 'patches')
    if not os.path.exists(patchesPath):
        print("-> Creating a folder 'patches' in your module. You should commit this to your\n"\
              "   repository in order for the installation script to work on other machines.")
        os.mkdir(patchesPath)
    else:
        print("-> Adding patches to the 'patches' folder in your module. Make sure to commit\n"\
              "   the newly added patches in order for the install script to work on other machines.")

    # function to get a non-used patch file name (makes sure not to overwrite anything)
    def getPatchFileName(patchesPath, moduleFolder, targetName):
        i = 1
        fileName = os.path.join(patchesPath, targetName + '.patch')
        while os.path.exists(fileName):
            fileName = os.path.join(patchesPath, targetName + '_' + str(i) + '.patch')
            i += 1

        # get relative path from module folder to patch
        return fileName, os.path.relpath(fileName, os.path.join(modParentPath, moduleFolder))

    for modFolder in patches.keys():

        # write the patches to the files
        unpubPatchPath, unpubPatchRelPath = getPatchFileName(patchesPath, modFolder, modFolder + '_unpublished')
        uncommPatchPath, uncommPatchRelPath = getPatchFileName(patchesPath, modFolder, modFolder + '_uncommitted')
        patches[modFolder]['unpublishedpatchrelpath'] = unpubPatchRelPath
        patches[modFolder]['uncommittedpatchrelpath'] = uncommPatchRelPath
        open(unpubPatchPath, 'w').write(patches[modFolder]['unpublished'])
        open(uncommPatchPath, 'w').write(patches[modFolder]['uncommitted'])

else:
    print("-> No Patches required")

##################################
# write installation shell script
instFileName = 'install_' + modName + '.sh' if not cmdArgs['filename'] else cmdArgs['filename']

with open(instFileName, 'w') as installFile:
    installFile.write('#'*100 + '\n')
    installFile.write('# This script installs the module "' + modName + '" together with all dependencies.\n')
    installFile.write('# Everything will be installed into a newly created sub-folder named "DUMUX".\n\n')

    installFile.write('mkdir DUMUX\n')
    installFile.write('cd DUMUX\n\n')

    # function to write installation procedure for a module
    def writeCloneModule(depModName, depModFolder):
        installFile.write('# ' + depModName + '\n')
        installFile.write('# ' + versions[depModFolder]['branch'] + ' # '
                               + versions[depModFolder]['revision'] + ' # '
                               + versions[depModFolder]['date'] + ' # '
                               + versions[depModFolder]['author'] + '\n')
        installFile.write('git clone ' + versions[depModFolder]['remote'] + '\n')
        installFile.write('cd ' + depModFolder + '\n')
        installFile.write('git checkout ' + versions[depModFolder]['branch'] + '\n')
        installFile.write('git reset --hard ' + versions[depModFolder]['revision'] + '\n')

        if depModFolder in patches:
            installFilewrite('git apply ' + patches[depModFolder]['unpublishedpatchrelpath'] + '\n')
            installFilewrite('git apply ' + patches[depModFolder]['uncommittedpatchrelpath'] + '\n')

        installFile.write('cd ..\n\n')

    # write the module clone first in order for the patches to be present
    writeCloneModule(modName, os.path.relpath(modPath, modParentPath))
    for depModName, depModFolder in zip(depNames, depFolders):
        if depModName != modName:
            writeCloneModule(depModName, depModFolder)

    # write configure command
    installFile.write('\n')
    installFile.write('./dune-common/bin/dunecontrol --opts=dumux.cmake.opts all\n')
