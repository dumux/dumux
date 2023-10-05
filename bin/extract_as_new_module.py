#!/usr/bin/env python3
# SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

# pylint: disable=redefined-outer-name

"""
This script extracts some specified applications into a separate Dune module.
For example make a dumux-pub repository accompanying a scientific paper.
"""

import os
import os.path
import subprocess
import argparse
import shutil
import textwrap
import itertools
import glob
import multiprocessing as mp
from shutil import copytree
from pathlib import Path
from functools import partial

from util.moduleinfo import getDependencies
from util.common import (
    callFromPath,
    runCommand,
    userQuery,
    queryYesNo,
    includedCppProjectHeaders,
    findMatchingFiles,
    versionTable,
)
from util.installscript import (
    makeInstallScript,
    supportedLanguages,
    getScriptExtension,
    filterDependencies,
    addDependencyVersions,
    addDependencyPatches,
)

MAIN_BRANCH_NAME = "main"


def readmeFileName():
    """The default readme filename"""
    return "README.md"


def makeStringList(items, indentation="   "):
    """Make a markdown list given the items"""
    return "\n".join([(indentation + "- " + str(it)) for it in items])


def replaceFileContent(file, newContent):
    """Replace the files content with the new content"""
    with open(file, "w") as newFile:
        newFile.write(newContent)


def appendFileContent(file, content):
    """Append some content to the given file"""
    with open(file, "a") as newFile:
        newFile.write(content)


def isInSubTree(file, base):
    """Check if a file is in the folder tree starting from base"""
    return Path(base).resolve() in Path(file).resolve().parents


def removeRedundantFolders(folders):
    """Remove folders that are duplicates or that are contained in a parent folder"""
    uniqueFolders = list(set(folders))
    return [sf for sf in uniqueFolders if not any(isInSubTree(sf, base) for base in uniqueFolders)]


def checkModuleFolder(moduleDirectory):
    """Verify a module folder (module to be extracted from)"""
    if not os.path.isdir(moduleDirectory):
        raise Exception(
            f"Module folder {moduleDirectory} not found. "
            f"Make sure to run this script from one level above {moduleDirectory}."
        )


def checkSubFolders(moduleDirectory, subDirectories):
    """Verify that subfolders exist in the module"""
    for subDirectory in subDirectories:
        path = Path(moduleDirectory) / Path(subDirectory)
        errMsg = f"Cannot handle the given folder {str(path)}"
        if not path.exists():
            raise Exception(errMsg + " because it does not exist.")
        if not path.is_dir():
            raise Exception(errMsg + " because it is not a directory.")


def extractSourceFiles(moduleDirectory, subFolder):
    """Find all applications (*.cc files) in the given subfolder"""
    sources = []
    for folder in subFolder:
        folderPath = os.path.abspath(os.path.join(moduleDirectory, folder))
        curSources = findMatchingFiles(folderPath, "*.cc")
        sources += [os.path.normpath(os.path.join(folderPath, s)) for s in curSources]

    if not sources:
        raise Exception(
            "No sources found in the provided subFolder.",
            f" Run '{os.path.abspath(__file__)} --help' for details.",
        )

    return sources


def runDuneProject():
    """Run the duneproject script"""
    duneProjectCommand = shutil.which("duneproject", path="dune-common/bin")
    if duneProjectCommand is None:
        raise IOError("Could not find 'duneproject' in dune-common/bin.")
    try:
        subprocess.run([duneProjectCommand], check=True)
    except Exception as exc:
        raise Exception("Error during creation of the new module.") from exc


def detectNewModule():
    """
    Detect a new module created by dunecontrol (using the timestamp)
    This has to be used right after runDuneProject.
    """

    print("\nDetecting the newly created module")
    directories = [d for d in os.listdir() if os.path.isdir(d)]
    newModule = max(directories, key=os.path.getmtime)

    isCorrectModule = queryYesNo(
        f"Found '{newModule}' to be the new module. Is this correct?", default=None
    )
    if not isCorrectModule:
        newModule = userQuery("Please provide the name of the new module:")

    duneModuleFile = f"{newModule}/dune.module"
    if not os.path.exists(duneModuleFile):
        print(f"Could not find module file {duneModuleFile}")
    return newModule


def copySubFolders(subFolder, oldPath, newPath):
    """Copy folders from old path to new path"""
    for sub in subFolder:
        copytree(os.path.join(oldPath, sub), os.path.join(newPath, sub))


def addFoldersToCMakeLists(modulePath, subFolder):
    """Make sure CMakeLists.txt files exist and contain the necessary add_subdirectory calls"""
    for subFolderPATH in subFolder:
        nestedLevel = subFolderPATH.count(os.sep) + 1
        cmakeListsFile = os.path.join(modulePath, "CMakeLists.txt")
        for i in range(nestedLevel):
            if os.path.exists(cmakeListsFile):
                with open(cmakeListsFile, "r") as cml:
                    lines = list(reversed(cml.readlines()))

                idx = -1
                for line in lines:
                    if line.startswith("add_subdirectory"):
                        idx = lines.index(line)
                        break

                newLines = lines[0:idx]
                newLines += [f"add_subdirectory({subFolderPATH.split(os.path.sep)[i]})\n"]
                newLines += lines[idx:]
                newContent = "".join(line for line in reversed(newLines))

                replaceFileContent(cmakeListsFile, newContent)
            else:
                with open(cmakeListsFile, "w") as cml:
                    cml.write(f"add_subdirectory({subFolderPATH.split(os.path.sep)[i]})\n")
            cmakeListsFile = os.path.join(
                os.path.dirname(cmakeListsFile),
                subFolderPATH.split(os.path.sep)[i],
                "CMakeLists.txt",
            )


def findHeaders(modulePath, sourceFiles):
    """Find header included (recursively) in the given source files"""
    with mp.Pool() as pool:
        headers = itertools.chain.from_iterable(
            pool.map(partial(includedCppProjectHeaders, projectBase=modulePath), sourceFiles)
        )
    return list(set(headers))


def copyFiles(filePaths, oldModulePath, newModulePath):
    """Copy files from the old to the new module"""
    newFiles = []
    for filePath in filePaths:
        newDirectory = os.path.join(
            newModulePath, os.path.relpath(os.path.dirname(filePath), oldModulePath)
        )
        newFile = os.path.join(newModulePath, os.path.relpath(filePath, oldModulePath))

        newFiles.append(newFile)
        os.makedirs(newDirectory, exist_ok=True)
        shutil.copy(filePath, newFile)
    return newFiles


def foldersWithoutSourceFiles(modulePath, checkSubFolder, sources):
    """Find those folders that do not contain a *.cc source file (application)"""
    sourceDirectories = [os.path.dirname(s) for s in sources]
    sourceDirectories = list(set(sourceDirectories))

    def hasChildSourceDirectory(directory):
        return any(isInSubTree(s, directory) for s in sourceDirectories)

    def isNotASourceDirectory(directory):
        return directory not in sourceDirectories and not hasChildSourceDirectory(directory)

    noSourceDirectories = []
    for sub in checkSubFolder:
        for root, dirs, _ in os.walk(os.path.join(modulePath, sub)):
            for directory in dirs:
                directory = os.path.join(root, directory)
                if isNotASourceDirectory(directory):
                    noSourceDirectories.append(directory)
    noSourceDirectories = removeRedundantFolders(noSourceDirectories)

    def removeEmptyParents():
        folderMap = {}
        for folder in noSourceDirectories:
            parent = os.path.dirname(folder)
            if parent not in folderMap:
                folderMap[parent] = []
            folderMap[parent].append(folder)

        for parent, folders in folderMap.items():
            found = set(folders)
            for root, dirs, _ in os.walk(parent):
                dirs = [os.path.join(root, d) for d in dirs]
                if set(dirs) == found and isNotASourceDirectory(parent):
                    for entry in found:
                        noSourceDirectories.remove(entry)
                    noSourceDirectories.append(parent)
        return noSourceDirectories

    return removeEmptyParents()


def removeFolder(modulePath, subFolder):
    """Remove a folder from the newly created module (adjust CMakeLists.txt)"""

    def removeAddSubdirectoryCommand(cmakeLists, folder):
        with open(cmakeLists, "r") as cml:
            content = cml.read()

        key = f"add_subdirectory({folder})\n"
        if key in content:
            replaceFileContent(cmakeLists, content.replace(key, ""))
            return True
        return False

    subFolderPath = os.path.abspath(os.path.join(modulePath, subFolder))
    subFolderName = os.path.basename(subFolderPath.rstrip(os.sep))

    mainCML = os.path.join(modulePath, "CMakeLists.txt")
    parentCML = os.path.join(subFolderPath, "../CMakeLists.txt")
    if not removeAddSubdirectoryCommand(mainCML, subFolder):
        if os.path.exists(parentCML):
            removeAddSubdirectoryCommand(parentCML, subFolderName)
    shutil.rmtree(subFolderPath)


def guideFolderDeletion(modulePath, candidates):
    """Interactive process to delete folders that may not be needed"""
    candidateList = makeStringList(candidates)
    print(
        "\n"
        "Could not automatically determine if the following directories\n"
        "contain data that is essential for the extracted applications:\n"
        "\n"
        f"{candidateList}\n"
    )

    deleted = []
    if queryYesNo(
        "Do you want to remove some of them (by choosing 'no' they are all preserved)?",
        default="no",
    ):
        for folder in foldersWithoutSources:
            if queryYesNo(f"Do you want to delete the folder {folder}?", default="yes"):
                relativePath = os.path.relpath(folder, modulePath)
                removeFolder(newModulePath, relativePath)
                deleted.append(relativePath)
    return deleted


def pubProjectURL(projectName):
    """Default URL for dumux-pub modules"""
    baseURL = "https://git.iws.uni-stuttgart.de/dumux-pub/"
    return baseURL + f"{projectName.lower()}.git"


def queryEmptyRemoteURL():
    """Interactive process to determine the remote URL of the new repository"""
    while True:
        if queryYesNo("Is your repository hosted in dumux-pub?"):
            project = input(
                "Please provide the name of your project "
                "(usually AuthorLastNameYearX, e.g. Luigi2020b)\n"
            )
            remote = pubProjectURL(project)
        else:
            remote = input("Please provide the URL of your repository:\n")

        try:
            print("Checking the repo (you may have to introduce credentials):")
            remoteContent = runCommand(f"git ls-remote {remote}", suppressTraceBack=True)
        except subprocess.CalledProcessError:
            print(f" - Could not find your repo at {remote}. ")
            print(" - Please revisit the provided information.")
            continue

        if remoteContent == "":
            return remote

        print(
            "- Remote repository is not empty. Cannot push to non-empty repository.\n"
            "- Please provide a URL of an empty/bare Git repository."
        )


def runGitCommand(path, cmd):
    """Helper function to call and log git commands"""
    print(f"Running {cmd} (you may have to provide your credentials)")
    callFromPath(path)(runCommand)(cmd)


def pushMainBranch(modulePath, remoteURL):
    """Push to the main branch of the new repository"""
    runGitCommand(modulePath, f"git push -u {remoteURL} {MAIN_BRANCH_NAME}")


def guideRepositoryInitialization(modulePath):
    """Initialize new module as git repo (and push if remote is provided)"""
    hasRepo = queryYesNo(
        "Do you have an empty remote repository to push the code to?"
        " (If not, we recommend that you create one and answer with 'yes')"
    )
    remoteURL = None if not hasRepo else queryEmptyRemoteURL()

    runGitCommand(modulePath, "git init")
    runGitCommand(modulePath, f"git switch -C {MAIN_BRANCH_NAME}")
    runGitCommand(modulePath, "git add .")
    runGitCommand(modulePath, 'git commit -m "Initial commit"')

    if hasRepo:
        runGitCommand(modulePath, f"git remote add origin {remoteURL}")
        pushMainBranch(modulePath, remoteURL)

    return remoteURL


def guideAddLicenseFile(modulePath, readme):
    """Optionally add a GPLv3 license to the new module"""
    addLicense = queryYesNo(
        "Do you want to add a GPLv3 license file as LICENSES/GPL-3.0-or-later.txt?"
    )
    if addLicense:
        moduleLicenseLocation = os.path.join(newModulePath, "LICENSES/GPL-3.0-or-later.txt")
        # download the GPL license file from the GNU website to the LICENSES directory
        os.makedirs(os.path.dirname(moduleLicenseLocation), exist_ok=True)
        licenseURL = "https://www.gnu.org/licenses/gpl-3.0.txt"
        downloadCommand = ["wget", "-O", moduleLicenseLocation, licenseURL]
        try:
            subprocess.run(downloadCommand, check=True)
        except subprocess.CalledProcessError:
            print(
                f"WARNING: Could not download the GPLv3 license file from {licenseURL}."
                " Skipping license file addition."
            )
            # remove the LICENSES directory if it only contains the empty license file
            fileCount = len(os.listdir(os.path.dirname(moduleLicenseLocation)))
            if os.path.exists(moduleLicenseLocation) and fileCount == 1:
                shutil.rmtree(os.path.dirname(moduleLicenseLocation))
            # remove only the license file if it is not the only file in the LICENSES directory
            elif os.path.exists(moduleLicenseLocation):
                os.remove(moduleLicenseLocation)
        else:
            print(
                "License file was downloaded and added to the new module."
                " Please check if the license file is correct."
            )
            # commit license to git
            runGitCommand(modulePath, "git add LICENSES/GPL-3.0-or-later.txt")
            runGitCommand(modulePath, 'git commit -m "Add GPLv3 as GPL-3.0-or-later.txt file"')
            # include license information to README.md
            appendFileContent(readme, infoReadmeLicense())
            runGitCommand(modulePath, "git add README.md")
            runGitCommand(modulePath, 'git commit -m "Add license information to README.md"')
    print("\n")


def dependenciesAndPatches(modulePath, skip=None):
    """Determine the module's dependencies"""
    try:
        print(
            f"Determining dependencies of the module: {modulePath}..."
            " this may take several minutes"
        )
        deps = getDependencies(modulePath)
        deps = filterDependencies(deps, skip or [])
        deps = addDependencyVersions(deps)
        deps = addDependencyPatches(deps)
    except Exception as exc:
        raise Exception("Error getting the dependencies.") from exc
    return deps


def guideVersionsReadme(modulePath, dependencies, readme, remoteURL=None):
    """Write detailed version information to the new readme file"""
    writeVersionInfo = queryYesNo(
        "Write detailed version information"
        f" (folder/branch/commits/dates) into {readmeFileName()}?"
    )
    if writeVersionInfo:
        table = versionTable(dependencies)
        appendFileContent(readme, f"\n## Version Information\n\n{table}\n")
        runGitCommand(modulePath, f'git commit {readme} -m "[readme] Update version information"')

        if remoteURL:
            pushMainBranch(modulePath, remoteURL)


def guideInstallScriptGeneration(modulePath, dependencies, scriptNameBody):
    """Interactively add an install script (or not)"""
    language = userQuery(
        "In which language would you like to generate the install script?"
        " (choose 'none' if you don't want an installation script)",
        supportedLanguages() + ["none"],
    )

    if language == "none":
        return None

    installScriptName = scriptNameBody + getScriptExtension(language)

    try:
        makeInstallScript(
            modPath=modulePath,
            dependencies=dependencies,
            scriptName=installScriptName,
            topFolderName="",
        )
    except Exception as exc:  # pylint: disable=broad-except
        print(f"Error during install script generation: {exc}")

    return installScriptName


def processInstallScript(script, modulePath, readme, remoteURL=None):
    """Add install script and installation instructions to the new module"""
    newScript = os.path.join(modulePath, script)
    shutil.move(script, newScript)
    subprocess.call(["chmod", "u+x", newScript])
    runGitCommand(modulePath, f"git add {script}")
    runGitCommand(modulePath, 'git commit -m "Add installation script"')

    appendFileContent(readme, infoReadmeInstallation(remoteURL, script, newModuleName))
    runGitCommand(
        modulePath,
        f'git commit {readme} -m "[readme] Update installation instructions"',
    )

    if remoteURL:
        pushMainBranch(modulePath, remoteURL)
    else:
        print(
            "\n"
            "Please remember to manually fix the installation instructions\n"
            "once you know the remote URL where your repository will be hosted"
        )


def noRemoteURLInfo(newModule):
    """Print information if not remote is present"""
    return textwrap.dedent(
        f"""\
        No remote URL given for new module {newModule}.
        Please remember to replace the placeholder `$remoteurl`
        in README.md of the new module {newModule} once a remote is available.
    """
    )


def infoInitial(moduleDirectory, subFolder, sourceFiles):
    """Some general information for users of the script"""

    sourcesList = makeStringList(sourceFiles)
    subFolderList = makeStringList(subFolder)

    return textwrap.dedent(
        """\
        This script will extract the following subfolder(s) of
        the module '{0}':
        {1}

        and all headers contained in '{0}'
        that are required to build the executables from the sources:
        {2}

        The extracted files are copied into a new DUNE module retaining the directory
        structure. The files required for creating a working DUNE module (such as
        CMakeLists.txt) will be created and/or updated.

        In the next step, the 'duneproject' script will be executed to guide the
        creation of your new DUNE module. Please answer all upcoming queries to the
        best of your knowledge.

        Important: the new module should NOT depend on the module '{0}'
    """
    ).format(moduleDirectory, subFolderList, sourcesList)


def infoReadmeMain(moduleDirectory, subFolder, sourceFiles):
    """Main part of the README.md document"""

    def relativePath(path):
        return os.path.relpath(path, moduleDirectory)

    subFolderString = "".join([f"*   `{d}`\n" for d in subFolder])
    sourceString = "".join([f"*   `{relativePath(s)}`\n" for s in sourceFiles])

    return textwrap.dedent(
        """\
        This file has been created automatically. Please adapt it to your needs.

        ## Content

        The content of this DUNE module was extracted from the module `{0}`.
        In particular, the following subFolder of `{0}` have been extracted:
        {1}

        Additionally, all headers in `{0}` that are required to build the
        executables from the sources
        {2}

        have been extracted. You can configure the module just like any other DUNE
        module by using `dunecontrol`. For building and running the executables,
        please go to the build folders corresponding to the sources listed above.

    """
    ).format(moduleDirectory, subFolderString, sourceString)


def infoReadmeLicense():
    """License part of the README.md document"""
    return textwrap.dedent(
        """\

        ## License

        This project is licensed under the terms and conditions of the GNU General Public
        License (GPL) version 3 or - at your option - any later version.
        The GPL can be found under [GPL-3.0-or-later.txt](LICENSES/GPL-3.0-or-later.txt)
        provided in the `LICENSES` directory located at the topmost of the source code tree.

    """
    )


def infoReadmeInstallation(remoteURL, installScriptName, newModuleName):
    """Installation part of the README.md document"""

    installScriptPath = os.path.join(newModuleName, installScriptName)
    remoteHints = ""
    if not remoteURL:
        remoteURL = "$remoteurl"
        remoteHints = "\nImportant: $remoteurl has to be corrected!\n"
    return textwrap.dedent(
        """\

        ## Installation

        The installation procedure is done as follows:
        Create a root folder, e.g. `DUMUX`, enter the previously created folder,
        clone this repository and use the install script `{0}`
        provided in this repository to install all dependent modules.
        {1}
        ```sh
        mkdir DUMUX
        cd DUMUX
        git clone {2} {3}
        ./{4}
        ```

        This will clone all modules into the directory `DUMUX`,
        configure your module with `dunecontrol` and build tests.

    """
    ).format(installScriptName, remoteHints, remoteURL, newModuleName, installScriptPath)


def infoFinal(newModuleName):
    """Print success message with information"""
    return textwrap.dedent(
        f"""\
        ========================================================================

        The module was extracted successfully!

        The extracted module is contained in the subfolder '{newModuleName}'.
        You can configure it with
        ./dune-common/bin/dunecontrol --opts={newModuleName}/cmake.opts --only={newModuleName} all
    """
    )


###############
# Main script #
if __name__ == "__main__":

    # set script parameters
    EPILOG = """
    -----------------------------------------------------------
    The script has to be called one level above moduleDirectory.
    At least one of the subFolder (FOLDER_1 [FOLDER_2 ...]) has
    to contain a source file *.cc of an executable for which
    you would like to timber a table in dumux-pub.)
    -----------------------------------------------------------

    Example usage:
    python3 dumux/bin/extract_as_new_module.py dumux-fracture appl test

    (extracts the subFolder appl and test from the module dumux-fracture)

    """

    parser = argparse.ArgumentParser(
        prog="extract_as_new_module.py",
        usage="./dumux/bin/extractmodule/extract_as_new_module.py"
        " module_dir SUBFOLDER_1 [SUBFOLDER_2 ...]",
        description="This script extracts subFolder of a given DUNE module"
        " into a new DUNE module.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=EPILOG,
    )
    parser.add_argument("module_dir", help="Module from which the subfolder is extracted")
    parser.add_argument("subfolder", nargs="+", help='subfolder(s) of "module_dir" to be extracted')

    # prepare input
    args = vars(parser.parse_args())
    moduleDirectory = args["module_dir"].strip(os.sep)
    modulePath = os.path.abspath(moduleDirectory)
    baseFolder = os.path.abspath(os.path.join(moduleDirectory, "../"))
    subFolder = removeRedundantFolders(args["subfolder"])

    # make sure module directory is valid
    checkModuleFolder(moduleDirectory)
    # make sure subdirectories exist
    checkSubFolders(moduleDirectory, subFolder)

    # find main source files (i.e. applications) to extract
    sourceFiles = extractSourceFiles(moduleDirectory, subFolder)

    # guide user through new module creation
    print(infoInitial(moduleDirectory, subFolder, sourceFiles))
    input("Please read the above carefully and press [Enter] to proceed or abort with [Ctrl-C]...")

    # duneproject creates a new Dune module
    runDuneProject()

    # find the new module (newest folder) created by duneproject
    # we could also parse the output of duneproject but this solution
    # seems easier to code
    newModuleName = detectNewModule()
    newModulePath = os.path.join(baseFolder, newModuleName)

    # prepare all data in new module
    copySubFolders(subFolder, modulePath, newModulePath)
    addFoldersToCMakeLists(newModulePath, subFolder)

    # find all headers necessary to build the applications
    headers = findHeaders(modulePath, sourceFiles)
    newHeaders = copyFiles(headers, modulePath, newModulePath)
    newSourceFiles = copyFiles(sourceFiles, modulePath, newModulePath)

    # copy the gitignore file from dumux
    copyFiles(
        [os.path.join(baseFolder, "dumux/.gitignore")],
        os.path.join(baseFolder, "dumux"),
        newModulePath,
    )

    # copy the cmake.opts file from dumux
    copyFiles(
        [os.path.join(baseFolder, "dumux/cmake.opts")],
        os.path.join(baseFolder, "dumux"),
        newModulePath,
    )

    # get the CMake macro file name of the new module
    cmakeMacroPattern = os.path.join(newModulePath, "cmake/modules/*Macros.cmake")
    sourceCMakeMacroName = glob.glob(cmakeMacroPattern)[0]
    # remove the generated CMake directory and copy the one from the source module
    shutil.rmtree(os.path.join(newModulePath, "cmake"))
    copySubFolders(["cmake"], modulePath, newModulePath)
    # get the name of the extracted CMake macro file name and rename it to match the new module name
    newModuleCMakeMacroFile = glob.glob(cmakeMacroPattern)[0]
    os.rename(newModuleCMakeMacroFile, sourceCMakeMacroName)

    # guide user through removal of possibly unnecessary folders
    foldersWithoutSources = foldersWithoutSourceFiles(
        newModulePath, subFolder, newHeaders + newSourceFiles
    )

    actualSubFolder = subFolder
    if foldersWithoutSources:
        deletedFolders = guideFolderDeletion(newModulePath, foldersWithoutSources)
        actualSubFolder = [s for s in subFolder if s not in deletedFolders]

    # remove stuff that is created when running duneproject
    if os.path.join(newModulePath, "dune") not in foldersWithoutSources:
        removeFolder(newModulePath, "dune")
    if os.path.join(newModulePath, "src") not in foldersWithoutSources:
        removeFolder(newModulePath, "src")

    # delete auto-generated README file
    duneReadme = os.path.join(newModulePath, "README")
    if os.path.exists(duneReadme):
        os.remove(duneReadme)

    # prepare new README.md file
    newReadme = os.path.join(newModulePath, readmeFileName())
    replaceFileContent(newReadme, infoReadmeMain(moduleDirectory, actualSubFolder, sourceFiles))

    # try to initialize repo (to use its URL in later steps)
    remoteURL = guideRepositoryInitialization(newModulePath)
    if not remoteURL:
        print(noRemoteURLInfo(newModuleName))

    # optionally add a license file
    guideAddLicenseFile(newModulePath, newReadme)

    # make install script & finalize readme
    deps = dependenciesAndPatches(newModulePath, [modulePath])
    if deps:
        guideVersionsReadme(newModulePath, deps, newReadme, remoteURL)

        installScript = "install_" + newModuleName
        installScript = guideInstallScriptGeneration(newModulePath, deps, installScript)
        if installScript is not None:
            processInstallScript(installScript, newModulePath, newReadme, remoteURL)
    else:
        print("No dependencies found. Skipping install script generation")

    print(infoFinal(newModuleName))
