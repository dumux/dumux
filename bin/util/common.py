# SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

"""
Helper functions used in several utility scripts (e.g. extract_module.py, ...)
"""

import os
import re
import sys
import fnmatch
import functools
import subprocess
import traceback
import shlex


TERM_FORMATTING = {
    "warning": "\033[1;31m",
    "highlight": "\033[1;34m",
    "reset": "\033[0m",
    "none": "",
}


def styledBotPrint(string, style="none", **kwargs):
    """Styled print to terminal when asking user for input"""
    sys.stdout.write("\n🤖 ")
    sys.stdout.write(TERM_FORMATTING[style])
    print(string, **kwargs)
    sys.stdout.write(TERM_FORMATTING["reset"])


def addPrefix(prefix, text, separator=" "):
    """Add prefix to a string"""
    return prefix + separator + text


def addPrefixToLines(prefix, text, separator=" "):
    """Add prefix every line of a multiline string (separated by endline character)"""
    return "\n".join(addPrefix(prefix, line, separator) for line in text.split("\n"))


def escapeCharacter(text, character, escCharacter="\\"):
    """Escape a given character with backslashes"""
    return text.replace(character, f"{escCharacter}{character}")


def escapeCharacters(text, characters, escCharacter="\\"):
    """Escape given characters with backslashes"""
    for char in characters:
        text = escapeCharacter(text, char, escCharacter)
    return text


def indent(text, indentation="  "):
    """Add space indentation to text"""
    text = text.split("\n")
    text = [indentation + line for line in text]
    return "\n".join(text)


def makeTable(dictList, config=None, padding=2):
    """Format as a table"""
    if config is None:
        config = {key: key for d in dictList for key in d}

    def getColWidth(row):
        return max(len(str(r)) for r in row) + padding * 2

    def getCol(key):
        return [config[key]] + [d.get(key, "") for d in dictList]

    widths = {key: getColWidth(getCol(key)) for key in config}

    def makeRow(rowValues):
        row = "|"
        for key in config:
            row += f"{rowValues.get(key, '').center(widths[key])}|"
        return row

    table = [makeRow({key: config[key] for key in config})]
    table.append("|" + "|".join("-" * widths[key] for key in config) + "|")
    table.extend(makeRow(row) for row in dictList)
    return "\n".join(table)


def getCommandErrorHints(command):
    """Maybe give a hint matching command"""
    if "git " in command:
        return (
            "It seems that a git command failed. Please check:\n"
            "    -- is the module registered as git repository?\n"
            "    -- is upstream defined for the branch?"
        )
    return None


def runCommand(command, check=True, suppressTraceBack=False, errorMessage=""):
    """execute a command and retrieve the output"""

    try:
        return subprocess.run(
            shlex.split(command), check=check, text=True, capture_output=True
        ).stdout
    except subprocess.CalledProcessError:
        eType, eValue, eTraceback = sys.exc_info()
        if suppressTraceBack:
            traceback.print_exception(eType, eType(errorMessage), None)
        elif errorMessage:
            traceback.print_exception(eType, eType(errorMessage), eTraceback)
        else:
            print("An error occurred during subprocess run:")
            print(f"-- command: {command}")
            print(f"-- folder: {os.getcwd()}")
            traceback.print_exception(eType, eValue, eTraceback)
            hints = getCommandErrorHints(command)
            if hints is not None:
                print(hints)
        return ""


def callFromPath(path):
    """decorator to call function from within the given path"""

    def decorateCallFromPath(callFunc):
        @functools.wraps(callFunc)
        def wrapCallFromPath(*args, **kwargs):
            curPath = os.getcwd()
            os.chdir(path)
            result = callFunc(*args, **kwargs)
            os.chdir(curPath)
            return result

        return wrapCallFromPath

    return decorateCallFromPath


def userQuery(query, choices=None):
    """query something from the user"""

    choicesString = ", ".join(str(c) for c in choices) if choices else ""
    querySuffix = f" (choices: {choicesString})\n" if choices else " "

    while True:
        styledBotPrint(f"{query.strip()}{querySuffix}", style="highlight")
        inp = input()

        if choices and inp not in choices:
            styledBotPrint(
                f"Invalid answer: '{inp}'. Choose from {choicesString}.", style="warning"
            )
        else:
            return inp


def queryYesNo(question, default="yes"):
    """query a yes/no answer from the user"""

    affirmative = ["yes", "y", "ye"]
    negative = ["no", "n"]

    def getChoices():
        return ", ".join(c for c in affirmative + negative)

    def isAffirmative(choice):
        return choice in affirmative

    def isNegative(choice):
        return choice in negative

    def isValid(choice):
        return isAffirmative(choice) or isNegative(choice)

    if default is not None and not isValid(default):
        raise ValueError(f"\nInvalid default answer: '{default}', choices: '{getChoices()}'\n")

    if default is None:
        prompt = " [y/n] "
    else:
        prompt = " [Y/n] " if isAffirmative(default) else " [y/N] "

    while True:
        styledBotPrint(f"{question.strip()}{prompt}", style="highlight", end="")
        choice = input().lower()

        if default is not None and choice == "":
            return isAffirmative(default)

        if not isValid(choice):
            styledBotPrint(
                f"Invalid answer: '{choice}'. Choose from '{getChoices()}'", style="warning"
            )
        else:
            return isAffirmative(choice)


def cppHeaderFilter():
    """
    Filter out source files that are not headers
    (consider everything except for config.h)
    """
    return lambda fileName: fileName != "config.h"


def includedCppProjectHeaders(file, projectBase, headers=None, headerFilter=cppHeaderFilter()):
    """get all project headers included by a cpp file"""
    headers = headers or []

    filePath = os.path.join(projectBase, file)
    if not os.path.exists(filePath):
        raise IOError(f"Cpp file {filePath} does not exist")

    with open(filePath, "r") as sourceFile:
        content = sourceFile.read()
        headerInBracket = re.findall(r"#include\s+<(.+?)>", content)
        headerInQuotation = re.findall(r'#include\s+"(.+?)"', content)

        def process(pathInProject):
            headerPath = os.path.join(projectBase, pathInProject)
            if os.path.exists(headerPath):
                if headerFilter(pathInProject):
                    if headerPath not in headers:
                        headers.append(headerPath)
                        includedCppProjectHeaders(headerPath, projectBase, headers, headerFilter)

        for header in headerInBracket:
            process(header)
        for header in headerInQuotation:
            absHeaderPath = os.path.join(os.path.dirname(file), header)
            projectPath = os.path.relpath(absHeaderPath, projectBase)
            process(projectPath)
    return headers


def findMatchingFiles(path, pattern):
    """find all files below the given folder that match the given pattern"""

    result = []
    for root, _, files in os.walk(path):
        relativeRootPath = os.path.relpath(root, path)
        for file in files:
            if fnmatch.fnmatch(file, pattern):
                result.append(os.path.join(relativeRootPath, file))
    return result


def isGitRepository(pathToRepo="."):
    """Check if git repository exists at given path"""
    try:
        run = callFromPath(pathToRepo)(runCommand)
        run("git status")
        return True
    except subprocess.CalledProcessError:
        return False


def getRemote(pathToRepo="."):
    """Check if git remote exists at given path"""
    run = callFromPath(pathToRepo)(runCommand)
    return run("git ls-remote --get-url").strip("\n")


def fetchRepo(remote, pathToRepo="."):
    """Fetch repo"""
    run = callFromPath(pathToRepo)(runCommand)
    run(f"git fetch {remote}")


def hasUntrackedFiles(pathToRepo="."):
    """Check for untracked (by git) file for given repo"""
    run = callFromPath(pathToRepo)(runCommand)
    return run("git ls-files --others --exclude-standard") != ""


def isPersistentBranch(branchName):
    """Check if a branch is what we consider persistent (default protected branches)"""
    if branchName in ("origin/master", "origin/main"):
        return True
    if branchName.startswith("origin/releases/"):
        return True
    return False


# get the most recent commit that also exists on remote master/release branch
# may be used to find a commit we can use as basis for a pub module
def mostRecentCommonCommitWithRemote(modFolderPath, branchFilter=isPersistentBranch):
    """Find most recent commit shared with remote"""
    run = callFromPath(modFolderPath)(runCommand)

    def findBranches(sha):
        candidates = run(f"git branch -r --contains {sha}").split("\n")
        candidates = [branch.strip().split(" ->")[0] for branch in candidates]
        return list(filter(branchFilter, candidates))

    revList = run("git rev-list HEAD").split("\n")
    for rev in revList:
        branches = findBranches(rev)
        if branches:
            return branches[0], rev

    raise RuntimeError(
        "Could not find suitable ancestor commit on a branch that matches the given filter"
    )


# function to extract persistent, remotely available git versions for all
def getPersistentVersions(modFolderPaths):
    """Get versions of last commit on a persistent branch"""
    result = {}
    for modFolderPath in modFolderPaths:

        if not isGitRepository(modFolderPath):
            raise Exception("Folder is not a git repository")

        result[modFolderPath] = {}
        result[modFolderPath]["remote"] = getRemote(modFolderPath)

        # update remote to make sure we find all upstream commits
        fetchRepo(result[modFolderPath]["remote"], modFolderPath)

        branch, rev = mostRecentCommonCommitWithRemote(modFolderPath)
        run = callFromPath(modFolderPath)(runCommand)

        result[modFolderPath]["revision"] = rev
        result[modFolderPath]["date"] = run(f"git log -n 1 --format=%ai {rev}").strip("\n")
        result[modFolderPath]["author"] = run(f"git log -n 1 --format=%an {rev}").strip("\n")

        # this may return HEAD if we are on some detached HEAD tree
        result[modFolderPath]["branch"] = branch

    return result


def getPatches(persistentVersions):
    """Generate patches"""
    result = {}
    for path, gitInfo in persistentVersions.items():
        run = callFromPath(path)(runCommand)

        uncommittedPatch = run("git diff")
        unpublishedPatch = run(f"git format-patch --stdout {gitInfo['revision']}")
        untrackedPatch = ""
        untrackedFiles = run("git ls-files --others --exclude-standard")
        binaryExtension = (
            ".png",
            ".gif",
            ".jpg",
            ".tiff",
            ".bmp",
            ".DS_Store",
            ".eot",
            ".otf",
            ".ttf",
            ".woff",
            ".rgb",
            ".pdf",
        )
        if untrackedFiles:
            for file in untrackedFiles.splitlines():
                if not str(file).endswith(binaryExtension):
                    untrackedPatch += run(f"git --no-pager diff /dev/null {file}", check=False)

        result[path] = {}
        result[path]["untracked"] = untrackedPatch if untrackedPatch else None
        result[path]["unpublished"] = unpublishedPatch if unpublishedPatch else None
        result[path]["uncommitted"] = uncommittedPatch if uncommittedPatch else None
    return result


DEFAULT_VERSION_TABLE_CONFIG = {
    "name": "module name",
    "branch": "branch name",
    "revision": "commit sha",
    "date": "commit date",
}


def versionTable(
    versions,
    config=None,
    padding=2,
):
    """Make a table containing module versions"""
    if config is None:
        config = DEFAULT_VERSION_TABLE_CONFIG
    return makeTable(versions, config, padding)


def printVersionTable(versions):
    """Print a table containing module versions"""
    print(versionTable(versions=versions))
