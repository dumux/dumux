
# Language-specific backends for install script generation

import os
import sys

try:
    path = os.path.split(os.path.abspath(__file__))[0]
    sys.path.append(os.path.join(path, '../bin/util'))

    from common import addPrefixToLines, indent
except Exception:
    sys.exit("Could not import commoon module")


def preparePatchForWrite(patch):
    patch = patch.rstrip('\n')
    patch = patch.replace('\"', '\\"')
    patch = patch.replace('"', '\"')
    patch = r'{}\n'.format(patch)
    return patch


class InstallScriptWriteInterface():
    def writeSheBang(self, file):
        raise NotImplementedError(
            "Derived class must implement 'writeSheBang'"
        )

    def writeImports(self, file):
        raise NotImplementedError(
            "Derived class must implement 'writeImports'"
        )

    def writeComment(self, comment, file):
        raise NotImplementedError(
            "Derived class must implement 'writeComment'"
        )

    def writeInstallation(self, dependencies, file):
        raise NotImplementedError(
            "Derived class must implement 'writeInstallation'"
        )

    def writeConfiguration(self, dependencies, file):
        raise NotImplementedError(
            "Derived class must implement 'writeConfiguration'"
        )


class InstallScriptWriterBash(InstallScriptWriteInterface):
    def __init__(self):
        super().__init__()
        self.installFunction = """
installModule()
{
    FOLDER=$1
    URL=$2
    BRANCH=$3
    REVISION=$4

    echo "Pulling module code from $URL"

    git clone $URL
    pushd $FOLDER
        git checkout $BRANCH
        git reset --hard $REVISION
    popd
}

"""

    def writeSheBang(self, file):
        file.write("#!/bin/bash\n")

    def writeImports(self, file):
        pass

    def writeComment(self, comment, file):
        comment = addPrefixToLines('#', comment)
        file.write(comment)

    def writeInstallation(self, dependencies, file):
        file.write(f"{self.installFunction}\n")
        for dep in dependencies:
            name = dep['name']
            file.write(f'echo "Installing {name}"\n')
            file.write('installModule {} {} {} {}\n'
                       .format(dep['folder'],
                               dep['remote'],
                               dep['branch'],
                               dep['revision']))
            file.write('\n')

        file.write('echo "Applying patches"\n')
        for dep in dependencies:
            def writePatch(folder, patch):
                file.write(f'echo "Applying patch in {folder}"\n')

                patch = preparePatchForWrite(patch)
                file.write('pushd {}\n'.format(folder))
                file.write('    echo "{}\n" > tmp.patch\n'.format(patch))
                file.write('    git apply tmp.patch\n')
                file.write('    rm tmp.patch\n')
                file.write('popd\n')
                file.write('\n')

            if dep['uncommitted'] is not None:
                writePatch(dep['folder'], dep['uncommitted'])
            if dep['unpublished'] is not None:
                writePatch(dep['folder'], dep['unpublished'])

    def writeConfiguration(self, opts, file):
        file.write('echo "configuring project"\n')
        file.write(
            f'./dune-common/bin/dunecontrol --opts={opts} all'
        )


class InstallScriptWriterPython(InstallScriptWriteInterface):
    def __init__(self, version=3, indent=' '*4):
        super().__init__()
        self.indentation = indent
        self.version = str(int(version))
        self.packages = ['subprocess']
        self.installFunction = """
def callFromFolder(folder, cmd):
    subprocess.call(cmd, cwd=folder)

def installModule(folder, url, branch, revision):
    print("Pulling module code from {}".format(url))

    callFromFolder(".", ['git', 'clone', url])
    callFromFolder(folder, ['git', 'checkout', branch])
    callFromFolder(folder, ['git', 'reset', '--hard', revision])

def applyPatch(folder, patch):
    print("Applying patch in {}".format(folder))

    callFromFolder(folder, ['git', 'apply', patch])

"""

    def writeSheBang(self, file):
        file.write("#!/usr/bin/env python{}\n".format(self.version))

    def writeImports(self, file):
        for pkg in self.packages:
            file.write(f"import {pkg}\n")

    def writeComment(self, comment, file):
        comment = addPrefixToLines('#', comment)
        file.write(comment)

    def writeInstallation(self, dependencies, file):
        file.write(f"\n{self.installFunction}\n")
        file.write('\n')
        file.write('if __name__ == "__main__":\n')

        for dep in dependencies:
            name = dep['name']
            self.writeIndented_(f'print("Installing {name}")', file)
            self.writeIndented_('installModule("{}", "{}", "{}", "{}")'
                                .format(dep['folder'],
                                        dep['remote'],
                                        dep['branch'],
                                        dep['revision']),
                                file)
            file.write('\n')

        self.writeIndented_('print("Applying patches")', file)
        for dep in dependencies:
            def writePatch(folder, patch):
                self.writeIndented_(
                    f'print("Applying patch in {folder}")', file
                )

                patch = preparePatchForWrite(patch)
                self.writeIndented_('patch = """', file)
                file.write('{}\n'.format(patch))
                file.write('"""\n')
                self.writeIndented_(
                    'applyPatch("{}", patch)'.format(folder), file
                )
                file.write('\n')

            if dep['uncommitted'] is not None:
                writePatch(dep['folder'], dep['uncommitted'])
            if dep['unpublished'] is not None:
                writePatch(dep['folder'], dep['unpublished'])

    def writeConfiguration(self, opts, file):
        self.writeIndented_('print("configuring project")', file)
        self.writeIndented_(
            'subprocess.call('
            f'["./dune-common/bin/dunecontrol", "--opts={opts}", "all"]'
            ')',
            file
        )

    def writeIndented_(self, text, file):
        text.rstrip('\n')
        file.write(indent(text, self.indentation))
        file.write('\n')
