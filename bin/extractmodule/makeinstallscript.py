#!/usr/bin/env python3

# script to generate an install script for a dune-module,
# accounting for non-published commits and local changes

import os
import sys
import argparse
import subprocess

from util import getPersistentVersions, versionTable, getPatches
from installscriptwriter import InstallScriptWriterBash
from installscriptwriter import InstallScriptWriterPython

try:
    path = os.path.split(os.path.abspath(__file__))[0]
    sys.path.append(os.path.join(path, '../bin/util'))

    from getmoduleinfo import getModuleInfo, getDependencies
except Exception:
    sys.exit('Could not import getModuleInfo')

if sys.version_info[0] < 3:
    sys.exit("\nError': Python3 required")


def supportedLanguages():
    return ['python', 'bash']


def getScriptExtension(language):
    assert language in supportedLanguages()
    ext = {
        'python': '.py',
        'bash': '.sh'
    }
    return ext[language]


def makeScriptWriter(language):
    if language == 'bash':
        return InstallScriptWriterBash()
    elif language == 'python':
        return InstallScriptWriterPython()
    raise ValueError(f'Could not create writer for language {language}')


def getDefaultScriptName(modName, language):
    return 'install_{}{}'.format(
        modName,
        getScriptExtension(language)
    )


def printProgressInfo(infoLines, indLevel=0):
    firstPrefix = '\n' + '--'*(indLevel+1)
    emptyPrefix = firstPrefix.replace('-', ' ').strip('\n')
    print(f"{firstPrefix} {infoLines[0]}")
    for line in infoLines[1:]:
        print(f"{emptyPrefix} {line}")


def filterDependencies(dependencies, skipFolders=[]):
    if not skipFolders:
        return dependencies
    else:
        return [
            dep for dep in dependencies if dep['folder'] not in skipFolders
        ]


def addDependencyVersions(dependencies, ignoreUntracked=False):
    def getKey(dependency):
        return os.path.abspath(dependency['folder'])

    versions = getPersistentVersions(
        [getKey(d) for d in dependencies], ignoreUntracked
    )
    if len(versions) != len(dependencies):
        raise Exception("Not all versions of all modules could be found.")

    mergedResult = []
    for depInfo in dependencies:
        versionInfo = versions[getKey(depInfo)]
        mergedResult.append({**depInfo, **versionInfo})
    return mergedResult


def addDependencyPatches(dependenciesWithVersions):
    def getKey(dependency):
        return os.path.abspath(dependency['folder'])

    patches = getPatches({
        getKey(d): d for d in dependenciesWithVersions
    })

    mergedResult = []
    for depInfo in dependenciesWithVersions:
        patch = patches[getKey(depInfo)]
        mergedResult.append({**depInfo, **patch})
    return mergedResult


def makeInstallScript(modPath,
                      dependencies,
                      scriptName,
                      writer,
                      topFolderName='DUMUX',
                      optsFile=None):

    modPath = os.path.abspath(modPath)
    modParentPath = os.path.abspath(os.path.join(modPath, '../'))
    modName = getModuleInfo(modPath, 'Module')

    # make sure we have relative paths for installation instructions
    def makeRelPath(p):
        return os.path.relpath(p, modParentPath)

    def modify(dep):
        result = dep
        result['folder'] = makeRelPath(result['folder'])
        return result
    dependencies = [modify(d) for d in dependencies]

    if not optsFile:
        optsFile = os.path.join(modParentPath, 'dumux/cmake.opts')
    optsFile = os.path.abspath(optsFile)
    optsRelPath = makeRelPath(optsFile)

    with open(scriptName, 'w') as script:
        writer.writeSheBang(script)

        script.write('\n')
        writer.writeImports(script)

        script.write('\n')
        writer.writeComment(
            "\n"
            f"This installs the module {modName} and its dependencies.\n"
            "The exact revivions used are listed in the table below.\n"
            "However, note that this script may also apply further patches.\n"
            "\n",
            script
        )

        script.write('\n')
        writer.writeComment(
            versionTable({d['folder']: d for d in dependencies}),
            script
        )

        script.write('\n')
        writer.writeInstallation(dependencies, script)

        script.write('\n')
        writer.writeConfiguration(optsRelPath, script)


def printFoundDependencies(deps):
    if len(deps) > 0:
        infoText = [
            "Found the following dependencies",
            "\t| {:^50} | {:^50} |".format('module name', 'module folder'),
            "\t" + 107*'-']
        for dep in deps:
            infoText.append(
                '\t| {:^50} | {:^50} |'.format(dep['name'], dep['folder'])
            )
        printProgressInfo(infoText)


def printFoundVersionInfo(dependenciesWithVersions):
    table = versionTable({
        d['folder']: d for d in dependenciesWithVersions
    })
    printProgressInfo(
        ["The following (remotely available) versions are used as a basis",
         "on top of which the required patches will be automatically created:",
         "\n{}".format(table)]
    )


def printInstallationInstruction(scriptName, remote, topFolderName=None):

    if topFolderName:
        description = f"""
This will create a folder `{topFolderName}`, clone all modules into it,
configure the entire project and build the contained applications"
"""
    else:
        description = """
This will clone all modules into the folder from which the script is called,
configure the entire project and build the contained applications"
"""

    if not remote:
        remote = "$REMOTE_URL"

    printProgressInfo(['Info:', f"""

You might want to put installation instructions into the README.md file of your
module, for instance:

     ## Installation
     The easiest way of installation is to use the script `{scriptName}`
     provided in this repository. Using `wget`, you can simply install all
     dependent modules by typing:

     ```sh
     wget {remote}/{scriptName}
     ./{scriptName}
     ```

     {description}
"""])


if __name__ == '__main__':

    ###################
    # parse arguments
    parser = argparse.ArgumentParser(
        description='This script generates an install script for your module, '
                    'taking into account non-published commits & changes.\n'
                    'This expects that all modules are git repositories and '
                    'have a remote origin URL defined.'
    )
    parser.add_argument('-p', '--path',
                        required=True,
                        help='The path to your dune module')
    parser.add_argument('-f', '--filename',
                        required=False,
                        help='File in which to write the install script')
    parser.add_argument('-i', '--ignoreuntracked',
                        required=False, action='store_true',
                        help='Use this to ignore untracked files present')
    parser.add_argument('-t', '--topfoldername',
                        required=False, default='DUMUX',
                        help='Folder that the install script creates upon'
                             'execution to install the modules in. If you '
                             'pass an empty string, no folder will be created '
                             'and installation happens in place.')
    parser.add_argument('-o', '--optsfile',
                        required=False,
                        help='Provide custom opts file to be used for project '
                             'configuration. Note that this file is required '
                             'to be contained and committed in the module or '
                             'one of its dependencies.')
    parser.add_argument('-s', '--skipfolders',
                        required=False, nargs='+',
                        help='a list of module folders to be skipped')
    parser.add_argument('-l', '--language',
                        required=False, default='python',
                        choices=supportedLanguages(),
                        help='Language in which to write the install script')

    cmdArgs = vars(parser.parse_args())

    modPath = cmdArgs['path']
    skipFolders = cmdArgs['skipfolders']
    if skipFolders:
        skipFolders = list(set(skipFolders))

    printProgressInfo(["Determining the module dependencies"])
    deps = getDependencies(modPath)
    deps = filterDependencies(deps, skipFolders)
    printFoundDependencies(deps)
    if not deps:
        sys.exit("No dependencies found. Exiting.")

    printProgressInfo(["Determining the module versions"])
    deps = addDependencyVersions(deps, cmdArgs.get('ignoreuntracked', False))
    printFoundVersionInfo(deps)

    printProgressInfo(["Making patches for unpublished & uncommited changes"])
    deps = addDependencyPatches(deps)

    # actual script generation
    modPath = os.path.abspath(modPath)
    modName = getModuleInfo(modPath, 'Module')
    printProgressInfo(
        ["Creating install script for module '{}' in folder '{}'"
         .format(modName, modPath)]
    )

    language = cmdArgs['language']
    scriptName = cmdArgs.get('filename', None)
    if not scriptName:
        scriptName = getDefaultScriptName(modName, language)

    makeInstallScript(
        modPath=modPath,
        dependencies=deps,
        scriptName=scriptName,
        writer=makeScriptWriter(language),
        topFolderName=cmdArgs.get('topfoldername', None),
        optsFile=cmdArgs.get('optsFile', None)
    )

    subprocess.call(['chmod', 'u+x', scriptName])

    for d in deps:
        if d['name'] == modName:
            remote = d['remote']

    printProgressInfo([f"Successfully created install script '{scriptName}'"])
    printInstallationInstruction(
        scriptName,
        remote,
        cmdArgs.get('topfoldername', None)
    )
