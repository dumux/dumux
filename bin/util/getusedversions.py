#!/usr/bin/env python3

import os
import argparse
from common import runCommand
from common import callFromPath

# print warning message for scanned folders that are not git repositories
def printNoGitRepoWarning(folderPath):
    print("Folder " + folderPath + " does not seem to be the top level" \
          "of a git repository and will be skipped. Make sure not to call " \
          "this script from a sub-directory of a git repository.")

# raise error due to untracked files present in the given module folder
def raiseUntrackedFilesError(folderPath):
    raise RuntimeError('Found untracked files in module folder: "' + folderPath + '". ' \
                       'Please commit, stash, or remove them.')

# returns true if the given folder is a git repository
def isGitRepository(modFolderPath):
    return os.path.exists(os.path.join(modFolderPath, '.git'))


# returns true if a module contains untracked files
def hasUntrackedFiles(modFolderPath):
    run = callFromPath(modFolderPath)(runCommand)
    return run('git ls-files --others --exclude-standard') != ''


# get the most recent commit that also exists on the remote master/release branch
# maybe used to find a commit we can use as basis for a pub module
def mostRecentCommonCommitWithRemote(modFolderPath, match=['/master', '/release']):
    run = callFromPath(modFolderPath)(runCommand)
    revList = run('git rev-list HEAD').split('\n')
    for rev in revList:
        remoteBranches = run('git branch -r --contains {}'.format(rev))
        for m in match:
            if m in remoteBranches:
                return rev
    
    raise RuntimeError('Could not suitable find ancestor commit'
                       ' that is synced with master on remote')


# function to extract git version information for modules
# returns a dictionary containing module information for each given module folder
def getUsedVersions(modFolderPaths, ignoreUntracked=False):

    result = {}
    for modFolderPath in modFolderPaths:

        # make sure this is the top level of a git repository
        if not isGitRepository(modFolderPath):
            printNoGitRepoWarning(modFolderPath)
        else:
            if not ignoreUntracked and hasUntrackedFiles(modFolderPath):
                raiseUntrackedFilesError(modFolderPath)

            run = callFromPath(modFolderPath)(runCommand)
            result[modFolderPath] = {}
            result[modFolderPath]['remote'] = run('git ls-remote --get-url').strip('\n')
            # update remote to make sure we find all upstream commits
            run('git fetch {}'.format(result[modFolderPath]['remote']))
            rev = mostRecentCommonCommitWithRemote(modFolderPaths)
            result[modFolderPath]['revision'] = rev
            result[modFolderPath]['date'] = run('git log -n 1 --format=%ai {}'.format(rev)).strip('\n')
            result[modFolderPath]['author'] = run('git log -n 1 --format=%an {}'.format(rev)).strip('\n')
            # this may return HEAD if we are on some detached HEAD tree
            result[modFolderPath]['branch'] = run('git rev-parse --abbrev-ref HEAD').strip('\n')

    return result

# create patches for unpublished commits and uncommitted changes in modules
def getPatches(modFolderPaths, ignoreUntracked=False):

    result = {}
    for modFolderPath in modFolderPaths:

        # make sure this is the top level of a git repository
        if not isGitRepository(modFolderPath):
            printNoGitRepoWarning(modFolderPath)
        else:
            if not ignoreUntracked and hasUntrackedFiles(modFolderPath):
                raiseUntrackedFilesError(modFolderPath)

            run = callFromPath(modFolderPath)(runCommand)
            unpubPatch = run('git format-patch --stdout @{upstream}')
            unCommPatch = run('git diff')
            if unpubPatch != '' or unCommPatch != '': result[modFolderPath] = {}
            if unpubPatch != '': result[modFolderPath]['unpublished'] = unpubPatch
            if unCommPatch != '': result[modFolderPath]['uncommitted'] = unCommPatch

    return result

# prints the detected versions as table
def printVersionTable(versions):
    print("\t| {:^50} | {:^50} | {:^50} | {:^30} |".format('module folder', 'branch', 'commit hash', 'commit date'))
    print("\t" + 193*'-')
    for folder, versionInfo in versions.items():
        print("\t| {:^50} | {:^50} | {:^50} | {:^30} |".format(folder, versionInfo['branch'], versionInfo['revision'], versionInfo['date']))


# For standalone execution
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='This script extracts the used dune/dumux versions.')
    parser.add_argument('-p', '--path', required=False, help='the path to the top folder containing your dune/dumux modules')
    parser.add_argument('-i', '--ignoreuntracked', required=False, action='store_true', help='use this flag to ignore untracked files present in the modules')
    parser.add_argument('-s', '--skipfolders', required=False, nargs='*', help='a list of module folders to be skipped')
    cmdArgs = vars(parser.parse_args())

    modulesPath = os.getcwd() if not cmdArgs['path'] else os.path.join(os.getcwd(), cmdArgs['path'])
    print('\nDetermining the versions of all dune modules in the folder: ' + modulesPath)

    def getPath(modFolder):
        return os.path.join(modulesPath, modFolder)

    modFolderPaths = [getPath(dir) for dir in os.listdir(modulesPath) if os.path.isdir(getPath(dir))]
    if cmdArgs['skipfolders']:
        cmdArgs['skipfolders'] = [f.strip('/') for f in cmdArgs['skipfolders']]
        modFolderPaths = [d for d in modFolderPaths if os.path.basename(d.strip('/')) not in cmdArgs['skipfolders']]
    versions = getUsedVersions(modFolderPaths, True)

    print("\nDetected the following versions:")
    printVersionTable(versions)

    # maybe check untracked files
    if not cmdArgs['ignoreuntracked']:
        modsWithUntracked = [f for f in versions if hasUntrackedFiles(f)]
        if modsWithUntracked:
            print('\n')
            print('#'*56)
            print('WARNING: Found untracked files in the following modules:\n\n')
            print('\n'.join(modsWithUntracked))
            print('\nPlease make sure that these are not required for your purposes.')
            print('If not, you can run this script with the option -i/--ignoreuntracked to suppress this warning.')
            print('#'*56)
