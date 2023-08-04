#!/usr/bin/env python3
# SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later


""""
Script to generate an install script for a dune-module,
accounting for non-published commits and local changes
"""

import sys
import os
import argparse
import subprocess

from util.moduleinfo import getDependencies, getModuleInfo
from util.installscript import (
    modulesWithUntrackedFiles,
    addDependencyPatches,
    addDependencyVersions,
    getDefaultScriptName,
    filterDependencies,
    makeInstallScript,
    printProgressInfo,
    printFoundDependencies,
    printFoundVersionInfo,
    printFinalMessage,
    supportedLanguages,
)


def runMakeInstallScript():
    """ "Generate an install script for a dune-module"""

    parser = argparse.ArgumentParser(
        description="This script generates an install script for your module, "
        "taking into account non-published commits & changes.\n"
        "This expects that all modules are git repositories and "
        "have a remote origin URL defined."
    )
    parser.add_argument("-p", "--path", required=True, help="The path to your dune module")
    parser.add_argument(
        "-f",
        "--filename",
        required=False,
        help="File in which to write the install script",
    )
    parser.add_argument(
        "-a",
        "--allow-untracked",
        required=False,
        action="store_true",
        help="Use this to include untracked files into patches",
    )
    parser.add_argument(
        "-t",
        "--topfoldername",
        required=False,
        default="DUMUX",
        help="Folder that the install script creates upon"
        "execution to install the modules in. If you "
        "pass an empty string, no folder will be created "
        "and installation happens in place.",
    )
    parser.add_argument(
        "-o",
        "--optsfile",
        required=False,
        help="Provide custom opts file to be used for project "
        "configuration. Note that this file is required "
        "to be contained and committed in the module or "
        "one of its dependencies.",
    )
    parser.add_argument(
        "-s",
        "--skipfolders",
        required=False,
        nargs="+",
        help="a list of module folders to be skipped",
    )
    parser.add_argument(
        "-l",
        "--language",
        required=False,
        default="python",
        choices=supportedLanguages(),
        help="Language in which to write the install script",
    )

    cmdArgs = vars(parser.parse_args())

    modPath = cmdArgs["path"]
    skipFolders = cmdArgs["skipfolders"]
    if skipFolders:
        skipFolders = list(set(skipFolders))

    printProgressInfo(["Determining the module dependencies"])
    deps = getDependencies(modPath, verbose=True, includeSelf=True)
    deps = filterDependencies(deps, skipFolders)
    printFoundDependencies(deps)
    if not deps:
        sys.exit("No dependencies found. Exiting.")

    printProgressInfo(["Determining the module versions"])
    deps = addDependencyVersions(deps)
    printFoundVersionInfo(deps)

    untracked = modulesWithUntrackedFiles(deps)
    if not cmdArgs.get("allow_untracked", False) and untracked:
        untrackedPaths = "\n".join(f" - {dep['path']}" for dep in untracked)
        raise Exception(
            "Found untracked files in the following modules:\n"
            f"{untrackedPaths}\n"
            "Please commit, stash, or remove them. Alternatively, you can "
            "use --allow-untracked to include them in the patches."
        )

    printProgressInfo(["Making patches for unpublished, uncommitted & untracked changes"])
    deps = addDependencyPatches(deps)

    # actual script generation
    modPath = os.path.abspath(modPath)
    modName = getModuleInfo(modPath, "Module")
    printProgressInfo([f"Creating install script for module '{modName}' in folder '{modPath}'"])

    scriptName = cmdArgs["filename"]
    if not scriptName:
        scriptName = getDefaultScriptName(modName, cmdArgs["language"])

    makeInstallScript(
        modPath=modPath,
        dependencies=deps,
        scriptName=scriptName,
        topFolderName=cmdArgs.get("topfoldername", None),
        optsFile=cmdArgs.get("optsFile", None),
    )

    subprocess.call(["chmod", "u+x", scriptName])
    printProgressInfo([f"Successfully created install script '{scriptName}'"])
    printFinalMessage(cmdArgs.get("topfoldername", None))


if __name__ == "__main__":
    runMakeInstallScript()
