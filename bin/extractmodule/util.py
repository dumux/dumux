#!/usr/bin/env python3

import os
import sys

try:
    path = os.path.split(os.path.abspath(__file__))[0]
    sys.path.append(os.path.join(path, '../util'))

    from common import callFromPath, runCommand
except Exception:
    sys.exit('Could not import common module')


def isGitRepository(pathToRepo='.'):
    try:
        run = callFromPath(pathToRepo)(runCommand)
        run('git status')
        return True
    except Exception:
        return False


def getRemote(pathToRepo='.'):
    run = callFromPath(pathToRepo)(runCommand)
    return run('git ls-remote --get-url').strip('\n')


def fetchRepo(remote, pathToRepo='.'):
    run = callFromPath(pathToRepo)(runCommand)
    run('git fetch {}'.format(remote))


def hasUntrackedFiles(pathToRepo='.'):
    run = callFromPath(pathToRepo)(runCommand)
    return run('git ls-files --others --exclude-standard') != ''


def isPersistentBranch(branchName):
    if branchName == 'origin/master':
        return True
    if branchName.startswith('origin/releases/'):
        return True
    return False


# get the most recent commit that also exists on remote master/release branch
# may be used to find a commit we can use as basis for a pub module
def mostRecentCommonCommitWithRemote(modFolderPath,
                                     branchFilter=isPersistentBranch):
    run = callFromPath(modFolderPath)(runCommand)

    def findBranches(sha):
        candidates = run('git branch -r --contains {}'.format(sha)).split('\n')
        candidates = [branch.strip().split(' ->')[0] for branch in candidates]
        return list(filter(branchFilter, candidates))

    revList = run('git rev-list HEAD').split('\n')
    for rev in revList:
        branches = findBranches(rev)
        if branches:
            return branches[0], rev

    raise RuntimeError('Could not find suitable ancestor commit'
                       ' on a branch that matches the given filter')


# function to extract persistent, remotely available git versions for all
def getPersistentVersions(modFolderPaths, ignoreUntracked=False):

    result = {}
    for modFolderPath in modFolderPaths:

        if not isGitRepository(modFolderPath):
            raise Exception('Folder is not a git repository')

        if hasUntrackedFiles(modFolderPath) and not ignoreUntracked:
            raise Exception("Found untracked files in '{}'."
                            "Please commit, stash, or remove them."
                            .format(modFolderPath))

        result[modFolderPath] = {}
        result[modFolderPath]['remote'] = getRemote(modFolderPath)

        # update remote to make sure we find all upstream commits
        fetchRepo(result[modFolderPath]['remote'], modFolderPath)

        branch, rev = mostRecentCommonCommitWithRemote(modFolderPath)
        run = callFromPath(modFolderPath)(runCommand)

        result[modFolderPath]['revision'] = rev
        result[modFolderPath]['date'] = run(
            'git log -n 1 --format=%ai {}'.format(rev)
        ).strip('\n')
        result[modFolderPath]['author'] = run(
            'git log -n 1 --format=%an {}'.format(rev)
        ).strip('\n')

        # this may return HEAD if we are on some detached HEAD tree
        result[modFolderPath]['branch'] = branch

    return result


def getPatches(persistentVersions):
    result = {}
    for path, gitInfo in persistentVersions.items():
        run = callFromPath(path)(runCommand)

        unCommPatch = run('git diff')
        unpubPatch = run(
            'git format-patch --stdout {}'.format(gitInfo['revision'])
        )

        if unpubPatch or unCommPatch:
            result[path] = {}
            if unpubPatch:
                result[path]['unpublished'] = unpubPatch
            if unCommPatch:
                result[path]['uncommitted'] = unCommPatch
    return result


def versionTable(versions):
    return "| {:<50} | {:<50} | {:<50} | {:<30} |\n".format(
        'module folder', 'branch', 'commit hash', 'commit date'
    ) + "| {:<50} | {:<50} | {:<50} | {:<30} |\n".format(
        '-'*50, '-'*50, '-'*50, '-'*30
    ) + "\n".join(
        ["| {:<50} | {:<50} | {:<50} | {:<30} |".format(
            folder, versionInfo['branch'], versionInfo['revision'], versionInfo['date']
        ) for folder, versionInfo in versions.items()]
    ) + "\n"


def printVersionTable(versions):
    print(versionTable(versions=versions))


def writeShellInstallScript(instFileName,
                            modName, modFolder,
                            folders, versions,
                            patches, patchModule, patchRelPath,
                            topFolderName, optsRelPath="dumux/cmake.opts"):
    """
    function to write the content into the generated shell install script

    Keyword Arguments:

        - instFileName -- The name of the generated shell install script

        - modName -- The name of the module to be installed
        - modFolder -- The folder containing the module to be installed

        - folders -- The list containing the folders of the module to be
            installed and all its dependencies
        - versions -- The persistent, remotely available git versions for
            the module to be installed and all its dependencies

        - patches -- The patches for unpublished commits and uncommitted changes
        - patchModule -- The paths for the modules which has unpublished commits and uncommited changes
        - patchRelPath -- The realative paths of the generated patch files

        - topFolderName -- The of the folder that the install script creates upon execution to install the module in.
            If an empty string is passed, no folder will be created and installation happens in place.
        - optsRelPath -- The relative path of custom opts file called by 'dunecontrol'

    """

    with open(instFileName, 'w') as installFile:
        installFile.write("#!/bin/bash\n\n")
        installFile.write("#"*80 + "\n")
        installFile.write("# This script installs the module '{}'"
                          " together with all dependencies.\n"
                          "\n".format(modName))

        # write function to install the new module and its dependencies into install script
        installFile.write(
            '# defines function to clone modules, set branches and commits\n'
            'installModule()\n'
            '{\n'
            '    URLS=$1\n'
            '    DEPFOLDERS=$2\n'
            '    DEPBRANCHES=$3\n'
            '    DEPSHAS=$4\n\n'
            '    for url in ${URLS[@]}; do\n'
            '        if ! git clone $url; then\n'
            '            echo "--Error: failed to clone $url"\n'
            '        fi\n'
            '    done\n\n'
            '    for index in ${!DEPFOLDERS[@]}; do\n'
            '        if ! cd ${DEPFOLDERS[index]}; then\n'
            '            echo "Error: could not enter folder ${DEPFOLDERS[index]}"\n'
            '        fi\n'
            '        if ! git checkout ${DEPBRANCHES[index]}; then\n'
            '            echo "-- Error: failed to check out branch  ${DEPBRANCHES[index]} in ${DEPFOLDERS[index]}."\n'
            '        fi\n'
            '        if ! git reset --hard ${DEPSHAS[index]}; then\n'
            '            echo "-- Error: failed to check out commit ${DEPSHAS[index]} in module ${DEPFOLDERS[index]}."\n'
            '        fi\n'
            '        cd ..\n'
            '    done\n'
            '}\n\n'

            '# defines function to apply patches\n'
            'applyPatch()\n'
            '{\n'
            '    PATCHES=$1\n'
            '    PATCHFOLDERS=$2\n\n'
            '    for i in ${!PATCHES[@]}; do\n'
            '        cd ${PATCHFOLDERS[i]}\n'
            '        if ! git apply ${PATCHES[i]}; then\n'
            '            echo "--Error: failed to apply patch $patch"\n'
            '        fi\n'
            '        cd ..\n'
            '    done\n'
            '}\n\n'
        )

        # write urls, depfolders/branches/shas and patches into install script
        installFile.write('# Information about url, folder/branch/sha and patch\n')
        installFile.write('URLS=(\n')
        for dep in folders:
            installFile.write(versions[dep]['remote'] + '\n')
        installFile.write(')\n\n')

        installFile.write('DEPFOLDERS=(\n')
        installFile.write('\n'.join(folders))
        installFile.write('\n)\n\n')

        installFile.write('DEPBRANCHES=(\n')
        for dep in folders:
            installFile.write(versions[dep]['branch'] + '\n')
        installFile.write(')\n\n')

        installFile.write('DEPSHAS=(\n')
        for dep in folders:
            installFile.write(versions[dep]['revision'] + '\n')
        installFile.write(')\n\n')

        # write patch information into install script
        installFile.write('PATCHFOLDERS=(\n')
        installFile.write('\n'.join(patchModule))
        installFile.write('\n)\n\n')

        installFile.write('PATCHES=(\n')
        installFile.write('\n'.join(patchRelPath))
        installFile.write('\n)\n\n')

        def writeCommandWithErrorCheck(command, errorMessage, indentationLevel=0):
            unitIndentation = ' '*4
            indent = unitIndentation*indentationLevel
            nextIndent = indent + unitIndentation
            installFile.write(
                f"{indent}if ! {command}; then\n"
                f"{nextIndent}echo \"{errorMessage}\"\n"
                "fi\n"
            )

        if topFolderName:
            installFile.write(
                '# Everything will be installed into a new folder "{}".\n'
                .format(topFolderName)
            )
            installFile.write(
                'echo "Creating the folder {} to install the modules in"\n'
                .format(topFolderName)
            )

            writeCommandWithErrorCheck(
                'mkdir -p {}'.format(topFolderName),
                '--Error: could not create top folder {}'.format(topFolderName)
            )

            writeCommandWithErrorCheck(
                'cd {}'.format(topFolderName),
                '--Error: could not enter top folder {}'.format(topFolderName))

            installFile.write('\n')
        else:
            installFile.write(
                '# Everything will be installed inside the folder from which the'
                ' script is executed.\n\n'
            )

        installFile.write('installModule $URLS $DEPFOLDERS $DEPBRANCHES $DEPSHAS\n\n')
        for depModPath in set(patchModule):
            if 'unpublished' in patches[depModPath]:
                installFile.write(f"cat >> {depModPath}/unpublished.patch <<'EOF'\n"
                                  + patches[depModPath]['unpublished']
                                  + "EOF\n")
            if 'uncommitted' in patches[depModPath]:
                installFile.write(f"cat >> {depModPath}/uncommitted.patch <<'EOF'\n"
                                  + patches[depModPath]['uncommitted']
                                  + "EOF\n")
        installFile.write('applyPatch $PATCHES $PATCHFOLDERS\n\n')

        # write configure command
        installFile.write('echo "-- All modules haven been cloned successfully. '
                          'Configuring project..."\n')
        writeCommandWithErrorCheck(
            './dune-common/bin/dunecontrol --opts={} all'.format(optsRelPath),
            '--Error: could not configure project'
        )

        # build the tests of the module
        installFile.write(
            '\n'
            'echo "-- Configuring successful. Compiling applications..."\n'
        )
        writeCommandWithErrorCheck(
            'cd {}/build-cmake'.format(modFolder),
            '--Error: could not enter build directory at {}/build-cmake'
            .format(modFolder)
        )
        writeCommandWithErrorCheck(
            'make build_tests',
            '--Error: applications could not be compiled. '
            'Please try to compile them manually.'
        )


def writePythonInstallScript(instFileName,
                             modName, modFolder,
                             folders, versions,
                             patches, patchModule, patchRelPath,
                             topFolderName, optsRelPath="dumux/cmake.opts"):
    """
    function to write the content into the generated python install script

    Keyword Arguments:

        - instFileName -- The name of the generated python install script

        - modName -- The name of the module to be installed
        - modFolder -- The folder containing the module to be installed

        - folders -- The list containing the folders of the module to be
            installed and all its dependencies
        - versions -- The persistent, remotely available git versions for
            the module to be installed and all its dependencies

        - patches -- The patches for unpublished commits and uncommitted changes
        - patchModule -- The paths for the modules which has unpublished commits and uncommited changes
        - patchRelPath -- The realative paths of the generated patch files

        - topFolderName -- The of the folder that the install script creates upon execution to install the module in.
            If an empty string is passed, no folder will be created and installation happens in place.
        - optsRelPath -- The relative path of custom opts file called by 'dunecontrol'

    """

    with open(instFileName, 'w') as installFile:
        # env and intro about script
        installFile.write("##!/usr/bin/env python3\n\n")
        installFile.write("'''\n")
        installFile.write("This script installs the module {} together with all dependencies.\n".format(modName))
        installFile.write("'''\n")

        # import libraies
        installFile.write("""
import os
import subprocess

""")

        # define funtions (git utilites, installModule, print messages)
        installFile.write("""
def show_message(message):
    print("*" * 120)
    print(message)
    print("*" * 120 + "\\n")


# execute a command in a folder and retrieve the output
def runCommandFromPath(command, path="."):
    curPath = os.getcwd()
    os.chdir(path)
    subprocess.run(command)
    os.chdir(curPath)


def git_clone(url):
    clone = ["git", "clone"]
    runCommandFromPath(command=[*clone, url])


def git_setbranch(branch):
    checkout = ["git", "checkout", branch]
    runCommandFromPath(command=checkout)


def git_setrevision(sha):
    reset = ["git", "reset", "--hard", sha]
    runCommandFromPath(command=reset)


def git_apply_patch(folder, patch):
    apply = ["git", "apply", patch]
    runCommandFromPath(command=apply, path=folder)


def installModule(urls, deps):
    for url in urls:
        git_clone(url)
    for folder, branch, sha in zip(deps["folders"], deps["branches"], deps["shas"]):
        os.chdir(folder)
        git_setbranch(branch)
        git_setrevision(sha)
        os.chdir("..")


def applyPatch(patches):
    for folder, patch in zip(patches["folders"], patches["files"]):
        git_apply_patch(folder, patch)\n
""")

        # write main function and step1: clone repositories, set branches and commits
        indent = ' '*4

        installFile.write("""
if __name__ == '__main__':
    ###########################################################################
    # (1/3) Clone repositories, set branches and commits
    ###########################################################################
    show_message("(1/3) Cloning repositories, setting branches and commits...")
""")
        if topFolderName:
            installFile.write('    os.makedirs("./DUMUX", exist_ok=True)\n'
                              '    os.chdir("DUMUX")\n\n')
        installFile.write(indent + 'urls = [\n')
        for dep in folders:
            installFile.write(indent*2 + "\"" + (versions[dep]['remote']) + "\"," + '\n')
        installFile.write(indent*2 + ']\n')

        installFile.write(indent + r'deps = {}' + '\n')

        installFile.write(indent + 'deps["folders"] = [\n')
        for dep in folders:
            installFile.write(indent*2 + "\"" + dep + "\"," + '\n')
        installFile.write(indent*2 + ']\n')

        installFile.write(indent + 'deps["branches"] = [\n')
        for dep in folders:
            installFile.write(indent*2 + "\"" + versions[dep]['branch'] + "\"," + '\n')
        installFile.write(indent*2 + ']\n')

        installFile.write(indent + 'deps["shas"] = [\n')
        for dep in folders:
            installFile.write(indent*2 + "\"" + versions[dep]['revision'] + "\"," + '\n')
        installFile.write(indent*2 + ']\n')
        installFile.write(indent + 'installModule(urls, deps)\n')
        installFile.write(indent + 'show_message("(1/3) Repositories are cloned and set properly.")')

        # step 2: generate and apply patches
        installFile.write("""
    ###########################################################################
    # (2/3) Generate and apply patches
    ###########################################################################
    show_message("(2/3) Creating patches for unpublished commits and uncommitted changes...")
""")
        for depModPath in set(patchModule):
            if 'unpublished' in patches[depModPath]:
                patchString = str(patches[depModPath]['unpublished']).replace('\\', r'\\').replace("'", "\\'").replace('"', '\\"')
                installFile.write(indent + "with open('{}/unpublished.patch', 'w') as patchFile:\n".format(depModPath))
                installFile.write(
                    indent*2 + "patchFile.write(\"\"\""
                    + patchString
                    + '\"\"\" )' + "\n"
                )
            if 'uncommitted' in patches[depModPath]:
                patchString = str(patches[depModPath]['uncommitted']).replace('\\', r'\\').replace("'", "\\'").replace('"', '\\"')
                installFile.write(
                    indent + "with open('{}/uncommitted.patch', 'w') as patchFile:\n".format(depModPath))
                installFile.write(
                    indent*2 + "patchFile.write(\"\"\""
                    + patchString
                    + '\"\"\" )' + "\n"
                    )

        installFile.write(indent + 'show_message("(2/3) Applying patches for unpublished commits and uncommitted changes...")\n')
        installFile.write(indent + r'patches = {}' + '\n')
        installFile.write(indent + 'patches["folders"] = [\n')
        for patchfolder in patchModule:
            installFile.write(indent*2 + "\"" + patchfolder + "\"," + '\n')
        installFile.write(indent*2 + ']\n')

        installFile.write(indent + 'patches["files"] = [\n')
        for patchPath in patchRelPath:
            installFile.write(indent*2 + "\"" + patchPath + "\"," + '\n')
        installFile.write(indent*2 + ']\n')
        installFile.write(indent + 'applyPatch(patches)\n')
        installFile.write("\n" + indent + 'show_message("(2/3) Step completed. All patch files are generated and applied.")')

        # step3: configure with dunecontrol and build tests
        installFile.write("""
    ###########################################################################
    # (3/3) Configure and build
    ###########################################################################
    show_message("(3/3) Configure and build dune modules and dumux using dunecontrol....")
    runCommandFromPath(command=["./dune-common/bin/dunecontrol", "--opts={0}", "all"])
    os.chdir("{1}/build-cmake")
    runCommandFromPath(command=["make", "build_tests"])
    show_message("(3/3) Step completed. Succesfully configured and built tests.")
""".format(*[optsRelPath, modFolder]))
