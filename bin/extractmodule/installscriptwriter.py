
# Language-specific backends for install script generation

import os
import sys

try:
    path = os.path.split(os.path.abspath(__file__))[0]
    sys.path.append(os.path.join(path, '../bin/util'))

    from common import addPrefixToLines
except Exception:
    sys.exit("Could not import required modules")


class InstallScriptWriterInterface():
    def __init__(self):
        self.ostream = None

    def setOutputStream(self, stream):
        self.ostream = stream

    def writeSheBang(self):
        raise NotImplementedError("Derived class must provide 'writeSheBang'")

    def writeComment(self, comment):
        raise NotImplementedError("Derived class must provide 'writeComment'")

    def writeMessageOutput(self, message):
        raise NotImplementedError(
            "Derived class must provide 'writeMessageOutput'"
        )

    def writeRuntimeArgsParse(self):
        raise NotImplementedError(
            "Derived class must provide 'writeRuntimeArgsParse'"
        )

    def writePreamble(self, topFolderName=None):
        pass

    def writeInstallation(self, dependency):
        raise NotImplementedError(
            "Derived class must implement 'writeInstallation'"
        )

    def writePatchApplication(self, folder, patchName):
        raise NotImplementedError(
            "Derived class must implement 'writePatch'"
        )

    def writeConfiguration(self, optsFile):
        raise NotImplementedError(
            "Derived class must implement 'writeConfiguration'"
        )


class InstallScriptWriterBash(InstallScriptWriterInterface):
    def __init__(self):
        super().__init__()

    def writeSheBang(self):
        self.ostream.write('#!/bin/bash\n')

    def writeComment(self, comment):
        comment = addPrefixToLines('#', comment)
        self.ostream.write(comment)

    def writeMessageOutput(self, message):
        self.ostream.write(f'echo "{message}"\n')

    def writePreamble(self, topFolderName=None):
        self.ostream.write("""
installModule()
{
    FOLDER=$1
    URL=$2
    BRANCH=$3
    REVISION=$4

    git clone $URL
    pushd $FOLDER
        git checkout $BRANCH
        git reset --hard $REVISION
    popd
}

applyPatch()
{
    FOLDER=$1
    PATCH=$2
    if ! test -f $PATCH; then
        echo "Patch $PATCH does not exist"
        exit 1
    fi

    pushd $FOLDER
        git apply $PATCH
    popd
}
""")
        top = topFolderName if topFolderName else "."
        self.ostream.write('TOP="{}"\n'.format(top))
        self.ostream.write('mkdir -p $TOP\n')
        self.ostream.write('cd $TOP\n')

    def writeRuntimeArgsParse(self):
        self.ostream.write("""
if [[ -z "$1" ]]; then
    PATCHFOLDER="$(dirname $(realpath -s $0))"
else
    PATCHFOLDER="$(realpath -s $1)"
fi
""")

    def writeInstallation(self, dependency):
        self.writeCommandWithErrorCheck_('installModule {} {} {} {}'
                                         .format(dependency['folder'],
                                                 dependency['remote'],
                                                 dependency['branch'],
                                                 dependency['revision']))

    def writePatchApplication(self, folder, patchName):
        self.writeCommandWithErrorCheck_(
            f'applyPatch {folder} $PATCHFOLDER/{patchName}'
        )

    def writeConfiguration(self, opts):
        self.writeCommandWithErrorCheck_(
            f'./dune-common/bin/dunecontrol --opts={opts} all'
        )

    def writeCommandWithErrorCheck_(self, cmd):
        self.ostream.write(f'if ! {cmd}; then exit 1; fi\n')


class InstallScriptWriterPython(InstallScriptWriterInterface):
    def __init__(self):
        super().__init__()

    def writeSheBang(self):
        self.ostream.write('#!/usr/bin/env python3\n')

    def writeComment(self, comment):
        comment = addPrefixToLines('#', comment)
        self.ostream.write(comment)

    def writeMessageOutput(self, message):
        self.ostream.write(f'print("{message}")\n')

    def writePreamble(self, topFolderName=None):
        top = topFolderName if topFolderName else "."
        self.ostream.write(f"""
import os
import sys
import subprocess

cwd = os.getcwd()
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
    runFromSubFolder(['git', 'clone', url], '.')
    runFromSubFolder(['git', 'checkout', branch], subFolder)
    runFromSubFolder(['git', 'reset', '--hard', revision], subFolder)

def applyPatch(subFolder, patch):
    patchPatch = os.path.abspath(patch)
    runFromSubFolder(['git', 'apply', patchPatch], subFolder)
""")

    def writeRuntimeArgsParse(self):
        self.ostream.write("""
patchFolder = "."
if len(sys.argv) > 1:
    patchFolder = sys.argv[1]
patchFolder = os.path.abspath(patchFolder)
""")

    def writeInstallation(self, dependency):
        self.ostream.write('installModule("{}", "{}", "{}", "{}")\n'
                           .format(dependency['folder'],
                                   dependency['remote'],
                                   dependency['branch'],
                                   dependency['revision']))

    def writePatchApplication(self, folder, patchName):
        self.ostream.write(
            f'applyPatch("{folder}", os.path.join(patchFolder, "{patchName}"))'
        )

    def writeConfiguration(self, opts):
        self.ostream.write(
            "runFromSubFolder(\n"
            f"    ['./dune-common/bin/dunecontrol', '--opts={opts}', 'all'],\n"
            "    '.'\n"
            "\n)"
        )
