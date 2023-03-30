# SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

"""Language-specific backends for install script generation"""

import textwrap
from abc import ABC, abstractmethod
from util.common import addPrefixToLines, escapeCharacters


def getRawString(text):
    """Get a raw string that can be written to file"""

    def makeRaw(text):
        return repr(text)

    def removeEnclosingQuotes(text):
        return text[1:-1]

    return removeEnclosingQuotes(makeRaw(text))


class InstallScriptWriterInterface(ABC):
    """Abstract writer interface to be implemented for dune module install script writers"""

    def __init__(self):
        """Initialize"""
        super().__init__()
        self.ostream = None

    def setOutputStream(self, stream):
        """Where do we write to?"""
        self.ostream = stream

    @abstractmethod
    def writeSheBang(self):
        """
        Write the she bang (first line of the script that
        specifies program executing the script)
        """

    @abstractmethod
    def writeComment(self, comment):
        """Write a code comment"""

    @abstractmethod
    def writeMessageOutput(self, message):
        """Write a message"""

    @abstractmethod
    def writePreamble(self, topFolderName=None):
        """Write the preamble of the script"""

    @abstractmethod
    def writeInstallation(self, dependency):
        """Write the installation process of module dependencies"""

    @abstractmethod
    def writePatchApplication(self, folder, patchContent):
        """Write the part that applies patches"""

    @abstractmethod
    def writeConfiguration(self, opts):
        """Write the configuration part"""


class InstallScriptWriterBash(InstallScriptWriterInterface):
    """Write a bash install script"""

    def writeSheBang(self):
        """Shebang for bash"""
        self.ostream.write("#!/bin/bash\n")

    def writeComment(self, comment):
        """Write a code comment"""
        comment = addPrefixToLines("#", comment)
        self.ostream.write(comment)

    def writeMessageOutput(self, message):
        """Write a message"""
        self.ostream.write(f'echo "{message}"\n')

    def writePreamble(self, topFolderName=None):
        """Write preable of the script (utility functions)"""
        self.ostream.write(
            textwrap.dedent(
                """\

            exitWithError()
            {
                MSG=$1
                echo "$MSG"
                exit 1
            }

            installModule()
            {
                FOLDER=$1
                URL=$2
                BRANCH=$3
                REVISION=$4

                if [ ! -d "$FOLDER" ]; then
                    if ! git clone $URL; then exitWithError "clone failed"; fi
                    pushd $FOLDER
                        if ! git checkout $BRANCH; then exitWithError "checkout failed"; fi
                        if ! git reset --hard $REVISION; then exitWithError "reset failed"; fi
                    popd
                else
                    echo "Skip cloning $URL since target folder "$FOLDER" already exists."
                fi
            }

            applyPatch()
            {
                FOLDER=$1
                PATCH=$2

                pushd $FOLDER
                    echo "$PATCH" > tmp.patch
                    if ! git apply tmp.patch; then exitWithError "patch failed"; fi
                    rm tmp.patch
                popd
            }

        """
            )
        )
        top = topFolderName if topFolderName else "."
        self.ostream.write(f'TOP="{top}"\n')
        self.ostream.write("mkdir -p $TOP\n")
        self.ostream.write("cd $TOP\n")

    def writeInstallation(self, dependency):
        """Write installation part of the script"""
        self.ostream.write(
            "installModule "
            f"{dependency['folder']} "
            f"{dependency['remote']} "
            f"{dependency['branch']} "
            f"{dependency['revision']}"
        )

    def writePatchApplication(self, folder, patchContent):
        """Write patch application part of the script"""

        def removeEscapedSingleQuotes(line):
            return line.replace(r"\'", "'")

        self.ostream.write('PATCH="\n')
        for line in patchContent.rstrip("\n").split("\n"):
            line = getRawString(line)
            line = removeEscapedSingleQuotes(line)
            line = escapeCharacters(line, ['"', "$", "`"])
            self.ostream.write(line)
            self.ostream.write("\n")
        self.ostream.write('"\n')
        self.ostream.write(f'applyPatch {folder} "$PATCH"')

    def writeConfiguration(self, opts):
        """Write configure part of the script"""
        self.ostream.write(
            f"if ! ./dune-common/bin/dunecontrol --opts={opts} all; then\n"
            '    echo "Configuration of the project failed"\n'
            "    exit 1\n"
            "fi\n"
        )


class InstallScriptWriterPython(InstallScriptWriterInterface):
    """Write a Python install script"""

    def writeSheBang(self):
        """Shebang for python3"""
        self.ostream.write("#!/usr/bin/env python3\n")

    def writeComment(self, comment):
        """Write a code comment"""
        comment = addPrefixToLines("#", comment)
        self.ostream.write(comment)

    def writeMessageOutput(self, message):
        """Write a message"""
        self.ostream.write(f'print("{message}")\n')

    def writePreamble(self, topFolderName=None):
        """Write the preamble of the script"""
        top = topFolderName if topFolderName else "."
        self.ostream.write(
            textwrap.dedent(
                f"""\

            import os
            import sys
            import subprocess

            top = "{top}"
            os.makedirs(top, exist_ok=True)


            def runFromSubFolder(cmd, subFolder):
                folder = os.path.join(top, subFolder)
                try:
                    subprocess.run(cmd, cwd=folder, check=True)
                except Exception as e:
                    cmdString = ' '.join(cmd)
                    sys.exit(
                        "Error when calling:\\n{{}}\\n-> folder: {{}}\\n-> error: {{}}"
                        .format(cmdString, folder, str(e))
                    )


            def installModule(subFolder, url, branch, revision):
                targetFolder = url.split("/")[-1]
                if targetFolder.endswith(".git"):
                    targetFolder = targetFolder[:-4]
                if not os.path.exists(targetFolder):
                    runFromSubFolder(['git', 'clone', url, targetFolder], '.')
                    runFromSubFolder(['git', 'checkout', branch], subFolder)
                    runFromSubFolder(['git', 'reset', '--hard', revision], subFolder)
                else:
                    print(
                        f"Skip cloning {{url}} since target '{{targetFolder}}' already exists."
                    )


            def applyPatch(subFolder, patch):
                sfPath = os.path.join(top, subFolder)
                patchPath = os.path.join(sfPath, 'tmp.patch')
                with open(patchPath, 'w') as patchFile:
                    patchFile.write(patch)
                runFromSubFolder(['git', 'apply', 'tmp.patch'], subFolder)
                os.remove(patchPath)
        """
            )
        )

    def writeInstallation(self, dependency):
        """Write installation part of the script"""
        self.ostream.write(
            "installModule("
            f'"{dependency["folder"]}", '
            f'"{dependency["remote"]}", '
            f'"{dependency["branch"]}", '
            f'"{dependency["revision"]}", '
            ")\n"
        )

    def writePatchApplication(self, folder, patchContent):
        """Write patch application part of the script"""
        self.ostream.write('patch = """\n')
        for line in patchContent.rstrip("\n").split("\n"):
            line = getRawString(line)
            self.ostream.write(escapeCharacters(line, ['"']))
            self.ostream.write("\n")
        self.ostream.write('"""\n')

        self.ostream.write(f'applyPatch("{folder}", patch)\n')

    def writeConfiguration(self, opts):
        """Write configure part of the script"""
        self.ostream.write(
            "runFromSubFolder(\n"
            f"    ['./dune-common/bin/dunecontrol', '--opts={opts}', 'all'],\n"
            "    '.'\n"
            ")\n"
        )
