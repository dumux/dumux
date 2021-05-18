#!/usr/bin/env python3

# script to generate an install script for a dune-module,
# accounting for non-published commits and local changes

import os
import sys
import argparse
import subprocess
from getmoduleinfo import *
from getusedversions import *

if sys.version_info[0] < 3:
    sys.exit("\nError': Python3 required")


###################
# parse arguments
parser = argparse.ArgumentParser(
    description='This script generates an install script for your dune module,' \
                'taking into account non-published commits and local changes.\n' \
                'This expects that your module is a git repository and that a ' \
                'remote origin exists and has been set already.'
)
parser.add_argument('-p', '--path', required=True, help='The path to your dune module')
parser.add_argument('-f', '--filename', required=False, help='Name of the file in which to write the install script')
parser.add_argument('-i', '--ignoreuntracked', required=False, action='store_true', help='Use this flag to ignore untracked files present in the modules')
parser.add_argument('-t', '--topfoldername', required=False, default='DUMUX',
                    help='Name of the folder that the install script creates upon execution to install the modules in. '\
                         'If you pass an empty string, no folder will be created and installation happens in place.')
parser.add_argument('-o', '--optsfile', required=False,
                    help='Provide custom opts file to be used for the call to dunecontrol. '\
                         'Note that this file is required to be contained and committed within the module or its dependencies.')
parser.add_argument('-s', '--skipfolders', required=False, nargs='*', help='a list of module folders to be skipped')
cmdArgs = vars(parser.parse_args())


####################
# Welcome message
modPath = os.path.abspath( os.path.join(os.getcwd(), cmdArgs['path']) )
modParentPath = os.path.abspath( os.path.join(modPath, '../') )
modFolder = os.path.relpath(modPath, modParentPath)

try: modName = getModuleInfo(modPath, 'Module')
except:
    sys.exit("\Error: Could not determine module name. Make sure the path to\n" \
             "        your module is correct and that a dune.module file is present.")

instFileName = 'install_' + modName + '.sh' if not cmdArgs['filename'] else cmdArgs['filename']
print("\n-- Creating install script '{}' for the module '{}' in the folder '{}'".format(instFileName, modName, modPath))


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
    sys.exit("Error: Could not determine dependencies. At least the module itself should appear.")

depNames = [dep['name'] for dep in deps]
depFolders = [dep['folder'] for dep in deps]
depFolderPaths = [os.path.abspath(os.path.join(modParentPath, d)) for d in depFolders]
if cmdArgs['skipfolders']:
    cmdArgs['skipfolders'] = [f.strip('/') for f in cmdArgs['skipfolders']]
    depFolders = [d for d in depFolders if d not in cmdArgs['skipfolders']]
    depNames = [getModuleInfo(d, 'Module') for d in depFolders]
    depFolderPaths = [d for d in depFolderPaths if os.path.basename(d.strip('/')) not in cmdArgs['skipfolders']]

#################################
# determine specific commits of all modules
print("\n-- Determining the module versions")
try: versions = getUsedVersions(depFolderPaths, cmdArgs['ignoreuntracked'])
except Exception as err:
    print('\nCaught exception: ' + str(err))
    if 'untracked' in str(err):
        print('If you are certain that the untracked files are not needed for the ' \
              'installation of the modules, run this script with the -i flag.')
    sys.exit(1)

if len(versions) != len(depNames):
    sys.exit("Error': Could not determine versions of all modules")

print("-> Detected the following versions")
printVersionTable(versions)

#################################
# create patches if necessary
print("\n-- Creating patches for unpublished commits and uncommitted changes")
patches = getPatches(depFolderPaths, cmdArgs['ignoreuntracked'])

# create patch files
if len(patches) > 0:
    patchesPath = os.path.join(modPath, 'patches')
    if not os.path.exists(patchesPath):
        os.mkdir(patchesPath)
        print("-> Created a folder 'patches' in your module. You should commit this to your\n"\
              "   repository in order for the installation script to work on other machines.")
    else:
        print("-> Adding patches to the 'patches' folder in your module. Make sure to commit\n"\
              "   the newly added patches in order for the install script to work on other machines.")

    # get a non-used patch file name (makes sure not to overwrite anything)
    def getPatchFileName(depModPath, targetName):
        i = 1
        fileName = os.path.join(patchesPath, targetName + '.patch')
        while os.path.exists(fileName):
            fileName = os.path.join(patchesPath, targetName + '_' + str(i) + '.patch')
            i += 1

        return fileName, os.path.relpath(fileName, depModPath)

    # function to write a new patch file (returns the path to the new file)
    def writeDepModPatchFile(depModPath, depModName, type):
        patchPath, patchRelPath = getPatchFileName(depModPath, depModName + '_' + type)
        patches[depModPath][type + '_relpath'] = patchRelPath
        open(patchPath, 'w').write(patches[depModPath][type])
        return patchPath

    print("-> Created patch files:")
    for depModPath in patches.keys():
        depModName = getModuleInfo(depModPath, "Module")
        if 'unpublished' in patches[depModPath]:
            print(' '*3 + writeDepModPatchFile(depModPath, depModName, 'unpublished'))
        if 'uncommitted' in patches[depModPath]:
            print(' '*3 + writeDepModPatchFile(depModPath, depModName, 'uncommitted'))

else:
    print("-> No Patches required")


##################################
# write installation shell script

# in the install script, we work with relative paths from the parent to module folder
def getModRelPath(path): return os.path.relpath(path, modParentPath)
versions = {getModRelPath(modPath): val for modPath, val in versions.items()}
patches = {getModRelPath(modPath): val for modPath, val in patches.items()}

topFolderName = cmdArgs['topfoldername']
optsRelPath = 'dumux/cmake.opts'
if cmdArgs['optsfile']:
    optsPath = os.path.abspath( os.path.join(os.getcwd(), cmdArgs['optsfile']) )
    optsRelPath = os.path.relpath(optsPath, modParentPath)

# TODO: add support for different install script languages (e.g. Python).
with open(instFileName, 'w') as installFile:
    installFile.write('#!/bin/bash\n\n')
    installFile.write('#'*80 + '\n')
    installFile.write('# This script installs the module "' + modName + '" together with all dependencies.\n')
    installFile.write('\n')

    exitFunc = 'exitWith'
    installFile.write('# defines a function to exit with error message\n')
    installFile.write(exitFunc + ' ()\n')
    installFile.write('{\n')
    installFile.write('    echo ' + r'"\n$1"' +'\n')
    installFile.write('    exit 1\n')
    installFile.write('}\n\n')

    # function to write a command with error check into the script file
    def writeCommand(command, errorMessage, indentationLevel = 0):
        installFile.write(' '*indentationLevel + 'if ! ' + command + ';')
        installFile.write(' then ' + exitFunc + ' "' + errorMessage + '"; fi\n')

    # write section on creating the top folder
    if topFolderName:
        installFile.write('# Everything will be installed into a newly created sub-folder named "' + topFolderName + '".\n')
        installFile.write('echo "Creating the folder ' + topFolderName + ' to install the modules in"\n')
        writeCommand('mkdir -p ' + topFolderName , '--Error: could not create top folder ' + topFolderName)
        writeCommand('cd ' + topFolderName, '--Error: could not enter top folder ' + topFolderName)
        installFile.write('\n')
    else:
        installFile.write('# Everything will be installed inside the folder from which the script is executed.\n\n')

    # function to write installation procedure for a module
    def writeCloneModule(depModName, depModFolder):
        installFile.write('# ' + depModName + '\n')
        installFile.write('# ' + versions[depModFolder]['branch'] + ' # '
                               + versions[depModFolder]['revision'] + ' # '
                               + versions[depModFolder]['date'] + ' # '
                               + versions[depModFolder]['author'] + '\n')

        writeCommand('git clone ' + versions[depModFolder]['remote'],
                     '-- Error: failed to clone ' + depModName + '.')
        writeCommand('cd ' + depModFolder,
                     '-- Error: could not enter folder ' + depModFolder + '.')
        writeCommand('git checkout ' + versions[depModFolder]['branch'],
                     '-- Error: failed to check out branch ' + versions[depModFolder]['branch'] + ' in module ' + depModName + '.')
        writeCommand('git reset --hard ' + versions[depModFolder]['revision'],
                     '-- Error: failed to check out commit ' + versions[depModFolder]['revision'] + ' in module ' + depModName + '.')

        # write section on application of patches
        def writeApplyPatch(patchRelPath):
            installFile.write('if [ -f ' + patchRelPath + ' ]; then\n')
            writeCommand('git apply ' + patchRelPath,'--Error: failed to apply patch ' + patchRelPath + ' in module ' + depModName + '.', 4)
            installFile.write('else\n')
            installFile.write(' '*4 + exitFunc + ' "--Error: patch ' + patchRelPath + ' was not found."\n')
            installFile.write('fi\n')

        if depModFolder in patches and 'unpublished_relpath' in patches[depModFolder]:
            writeApplyPatch(patches[depModFolder]['unpublished_relpath'])
        if depModFolder in patches and 'uncommitted_relpath' in patches[depModFolder]:
            writeApplyPatch(patches[depModFolder]['uncommitted_relpath'])

        installFile.write('echo "-- Successfully set up the module ' + depModName + r'\n"' + '\n')
        installFile.write('cd ..\n\n')

    # write the module clone first in order for the patches to be present
    if modName in depNames:
        writeCloneModule(modName, modFolder)
    for depModName, depModFolder in zip(depNames, depFolders):
        if depModName != modName:
            writeCloneModule(depModName, depModFolder)

    # write configure command
    installFile.write('echo "-- All modules haven been cloned successfully. Configuring project..."\n')
    writeCommand('./dune-common/bin/dunecontrol --opts=' + optsRelPath + ' all', '--Error: could not configure project')

    # build the tests of the module
    installFile.write('\n')
    installFile.write('echo "-- Configuring successful. Compiling applications..."\n')
    writeCommand('cd ' + modFolder + "/build-cmake", '--Error: could not enter build directory at ' + modFolder + '/build-cmake')
    writeCommand('make build_tests', '--Error: applications could not be compiled. Please try to compile them manually.')

print("\n-- Successfully created install script file " + instFileName)

if len(patches) > 0:
    print("-> It is recommended that you now commit and publish the 'patches' folder and this install script in your module such that others can use it.")
    print("   IMPORTANT: After you committed the patches, you have to adjust the line of the install script in which your module is checked out to a specific commit.")
    print("              That is, in the line 'git reset --hard COMMIT_SHA' for your module, replace COMMIT_SHA by the commit in which you added the patches.")
    print("              If patches had to be created for your own module, please think about comitting and pushing your local changes and rerunning this script again.")

if modName in depNames: # print gudience to installation if the module is not skipped
    print("\n-- You might want to put installation instructions into the README.md file of your module, for instance:\n")
    print("   ## Installation\n")
    print("   The easiest way of installation is to use the script `" + instFileName + "` provided in this repository.")
    print("   Using `wget`, you can simply install all dependent modules by typing:\n")
    print("   ```sh")
    print("   wget " + versions[modFolder]['remote'] + "/" + instFileName)
    print("   chmod u+x " + instFileName)
    print("   ./" + instFileName)
    print("   ```\n")

    if topFolderName: print("   This will create a sub-folder `" + topFolderName + "`, clone all modules into it, configure the entire project and build the applications contained in this module.")
    else: print("   This will clone all modules into the folder from which the script is called, configure the entire project and build the applications contained in this module.")
