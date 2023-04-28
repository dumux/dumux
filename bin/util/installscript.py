# SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

""""
Helper functions to generate an install script for a dune-module,
accounting for non-published commits and local changes
"""

import os
import sys
import textwrap

from util.common import getPersistentVersions, hasUntrackedFiles, versionTable, getPatches
from util.moduleinfo import getModuleInfo
from util.installscript_writer import InstallScriptWriterBash
from util.installscript_writer import InstallScriptWriterPython

if sys.version_info[0] < 3:
    sys.exit("\nError': Python3 required")


def supportedLanguages():
    """Supported languages for the install script output"""
    return ["python", "bash"]


def getScriptExtension(language):
    """Default script extension for the given language"""
    assert language in supportedLanguages()
    ext = {"python": ".py", "bash": ".sh"}
    return ext[language]


def getScriptLanguageFromExtension(ext):
    """Default script language for the given extension"""
    language = {".py": "python", ".sh": "bash"}
    return language[ext]


def makeScriptWriter(language):
    """Create a new install script writer instance"""
    if language == "bash":
        return InstallScriptWriterBash()
    if language == "python":
        return InstallScriptWriterPython()
    raise ValueError(f"Could not create writer for language {language}")


def getDefaultScriptName(modName, language):
    """The default script name"""
    return f"install_{modName}{getScriptExtension(language)}"


def printProgressInfo(infoLines, indLevel=0):
    """Inform user about progress"""
    firstPrefix = "\n" + "--" * (indLevel + 1)
    emptyPrefix = firstPrefix.replace("-", " ").strip("\n")
    print(f"{firstPrefix} {infoLines[0]}")
    for line in infoLines[1:]:
        print(f"{emptyPrefix} {line}")


def filterDependencies(dependencies, skipFolders=None):
    """Filter dependencies to skip given folders"""
    if skipFolders is None:
        return dependencies

    def skipFolder(folderName):
        return any(folderName == os.path.basename(path) for path in skipFolders)

    return [dep for dep in dependencies if not skipFolder(dep["folder"])]


def addDependencyVersions(dependencies):
    """Add version info to all dependencies"""

    def getKey(dependency):
        return dependency["path"]

    versions = getPersistentVersions([getKey(d) for d in dependencies])
    if len(versions) != len(dependencies):
        raise Exception("Not all versions of all modules could be found.")

    mergedResult = []
    for depInfo in dependencies:
        versionInfo = versions[getKey(depInfo)]
        mergedResult.append({**depInfo, **versionInfo})
    return mergedResult


def modulesWithUntrackedFiles(dependencies):
    """Find modules with untracked files"""
    return [dep for dep in dependencies if hasUntrackedFiles(dep["path"])]


def addDependencyPatches(dependenciesWithVersions):
    """Add patch info to all dependencies"""

    def getKey(dependency):
        return dependency["path"]

    patches = getPatches({getKey(d): d for d in dependenciesWithVersions})

    mergedResult = []
    for depInfo in dependenciesWithVersions:
        patch = patches[getKey(depInfo)]
        mergedResult.append({**depInfo, **patch})
    return mergedResult


def makeInstallScript(modPath, dependencies, scriptName, topFolderName="DUMUX", optsFile=None):
    """Main driver: create installation script for a dune module"""
    _, extension = os.path.splitext(scriptName)
    writer = makeScriptWriter(getScriptLanguageFromExtension(extension))
    modPath = os.path.abspath(modPath)
    modName = getModuleInfo(modPath, "Module")

    modOptsFile = f"{modPath}/cmake.opts"
    if not optsFile:
        if os.path.isfile(modOptsFile):
            optsFile = f"{os.path.relpath(modPath)}/cmake.opts"
        else:
            optsFile = "dumux/cmake.opts"
    if os.path.isabs(optsFile):
        raise ValueError("Opts file must be given as relative path")
    if not any(optsFile.startswith(d["folder"]) for d in dependencies):
        print("Warning: opts file is not contained in any of the dependencies")

    with open(scriptName, "w") as script:

        writer.setOutputStream(script)
        writer.writeSheBang()

        script.write("\n")
        writer.writeComment(
            textwrap.dedent(
                f"""\

            This installs the module {modName} and its dependencies.
            The exact revisions used are listed in the table below.
            However, note that this script may also apply further patches.
            If so, all patches are required to be the current folder, or,
            in the one that you specified as argument to this script.

        """
            )
        )

        script.write("\n")
        writer.writeComment(versionTable(dependencies))

        script.write("\n")
        writer.writePreamble(topFolderName)

        for dep in dependencies:
            script.write("\n")
            writer.writeMessageOutput(f"Installing {dep['name']}")
            writer.writeInstallation(dep)

        def writePatch(patch, description, moduleName, folder):
            script.write("\n")
            writer.writeMessageOutput(f"Applying patch for {description} in {moduleName}")
            writer.writePatchApplication(folder, patch)

        for dep in dependencies:
            if dep["untracked"] is not None:
                description = "untracked files"
                writePatch(dep["untracked"], description, dep["name"], dep["folder"])
            if dep["unpublished"] is not None:
                description = "unpublished commits"
                writePatch(dep["unpublished"], description, dep["name"], dep["folder"])
            if dep["uncommitted"] is not None:
                description = "uncommitted changes"
                writePatch(dep["uncommitted"], description, dep["name"], dep["folder"])

        script.write("\n")
        writer.writeMessageOutput("Configuring project")
        writer.writeConfiguration(optsFile)


def printFoundDependencies(deps):
    """Output found dependencies"""
    if len(deps) > 0:
        infoText = ["Found the following dependencies"]
        infoText.extend(versionTable(deps, {"name": "module name", "path": "folder"}).split("\n"))
        printProgressInfo(infoText)


def printFoundVersionInfo(dependenciesWithVersions):
    """Output found versions"""
    table = versionTable(dependenciesWithVersions)
    printProgressInfo(
        [
            "The following (remotely available) versions are used as a basis",
            "on top of which the required patches will be automatically created:",
            f"\n{table}",
        ]
    )


def printFinalMessage(topFolderName=None):
    """Final message after the install script has been created"""
    if topFolderName:
        description = textwrap.dedent(
            f"""\
            Running this script will create a folder `{topFolderName}`, clone all modules
            into it, configure the entire project and build the contained applications
        """
        )
    else:
        description = textwrap.dedent(
            """\
            Running this script will clone all modules into the folder from which it is
            called, configure the entire project and build the contained applications
        """
        )

    printProgressInfo(["Info:", description])
