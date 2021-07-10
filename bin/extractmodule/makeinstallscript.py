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

    from common import getCurrentTimeStamp, addPrefixToLines, indent
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
        return dependency['path']

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
        return dependency['path']

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
    modName = getModuleInfo(modPath, 'Module')

    if not optsFile:
        optsFile = 'dumux/cmake.opts'
    if os.path.isabs(optsFile):
        raise ValueError("Opts file must be given as relative path")
    if not any(optsFile.startswith(d['folder']) for d in dependencies):
        print("Warning: opts file is not contained in any of the dependencies")

    patches = []
    scriptFolder = os.path.dirname(scriptName)

    def makePatchFile(patchBaseName, content):
        now = getCurrentTimeStamp()
        now = now.replace(' ', '_')
        patchName = f'{patchBaseName}_{now}'
        patchPath = os.path.join(scriptFolder, patchName)

        if os.path.exists(patchPath):
            raise IOError(f"Patch {patchPath} because it already exists")
        with open(patchPath, 'w') as patch:
            patch.write(content)

        patches.append(patchPath)
        return patchName

    with open(scriptName, 'w') as script:

        writer.setOutputStream(script)
        writer.writeSheBang()

        script.write('\n')
        writer.writeComment(
            "\n"
            f"This installs the module {modName} and its dependencies.\n"
            "The exact revisions used are listed in the table below.\n"
            "However, note that this script may also apply further patches.\n"
            "If so, all patches are required to be the current folder, or, \n."
            "in the one that you specified as argument to this script"
            "\n"
        )

        script.write('\n')
        writer.writeComment(
            versionTable({d['folder']: d for d in dependencies})
        )

        script.write('\n')
        writer.writePreamble(topFolderName)

        script.write('\n')
        writer.writeRuntimeArgsParse()

        for dep in dependencies:
            script.write('\n')
            writer.writeMessageOutput('Installing {}'.format(dep['name']))
            writer.writeInstallation(dep)

            def writePatch(dep, patchKey):
                baseName = '{}_{}'.format(dep['name'], patchKey)
                patchName = makePatchFile(baseName, dep[patchKey])

                script.write('\n')
                writer.writeMessageOutput(f'Applying patch {patchName}')
                writer.writePatchApplication(dep['folder'], patchName)

            if dep['unpublished'] is not None:
                writePatch(dep, 'unpublished')
            if dep['uncommitted'] is not None:
                writePatch(dep, 'uncommitted')

        script.write('\n')
        writer.writeMessageOutput('Configuring project')
        writer.writeConfiguration(optsFile)

    return {'patches': patches}


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


def printFinalMessage(scriptName,
                      patches=[],
                      topFolderName=None):

    if topFolderName:
        description = f"""
Running this script will create a folder `{topFolderName}`, clone all modules
into it, configure the entire project and build the contained applications"
"""
    else:
        description = """
Running this script will clone all modules into the folder from which it is
called, configure the entire project and build the contained applications"
"""

    if not patches:
        printProgressInfo(['Info:', description])
    else:
        patchesList = '\n'.join(patches)
        patchesList = addPrefixToLines('- ', patchesList)
        patchesList = indent(patchesList)
        printProgressInfo(['Info:', f"""
{description}

Note that this script requires several patches:
{patchesList}

Make sure that these are in the same folder as the script or pass the folder
containing them as runtime argument.
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
    deps = getDependencies(modPath, verbose=True)
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

    config = makeInstallScript(
        modPath=modPath,
        dependencies=deps,
        scriptName=scriptName,
        writer=makeScriptWriter(language),
        topFolderName=cmdArgs.get('topfoldername', None),
        optsFile=cmdArgs.get('optsFile', None)
    )

    subprocess.call(['chmod', 'u+x', scriptName])
    printProgressInfo([f"Successfully created install script '{scriptName}'"])
    printFinalMessage(
        scriptName,
        config['patches'],
        cmdArgs.get('topfoldername', None)
    )
