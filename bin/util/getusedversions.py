#!/usr/bin/env python3

import os
import subprocess

# print warning message for scanned folders that are not git repositories
def printNoGitRepoWarning():
    print("Folder " + modFolder + " does not seem to be the top level" \
          "of a git repository and will be skipped. Make sure not to call " \
          "this script from a sub-directory of a git repository.")

# function to run a shell command and retrieve the output
def runCommand(command):
    return subprocess.run(command, shell=True, check = True,
                                   text=True, capture_output=True).stdout

# check for untracked files
def checkUntrackedFiles(modFolder):
    if runCommand('git ls-files --others --exclude-standard') != '':
        raise RuntimeError('Found untracked files in module folder: "' + modFolder + '". ' \
                           'Please commit, stash, or remove them.')

# function to extract Git version information for modules
def getUsedVersions(modFolders=None):
    curPath = os.getcwd()
    result = {}

    if modFolders == None:
        modFolders = [dir for dir in os.listdir(curPath) if os.path.isdir(dir)]

    for modFolder in modFolders:
        os.chdir(os.path.join(curPath, modFolder))

        # make sure this is the top level of a git repository
        if not os.path.exists( os.path.join(os.getcwd(), '.git') ):
            printNoGitRepoWarning()
        else:
            # checkUntrackedFiles(modFolder)

            result[modFolder] = {}
            result[modFolder]['remote'] = runCommand('git ls-remote --get-url').strip('\n')
            result[modFolder]['revision'] = runCommand('git log -n 1 --format=%H @{upstream}').strip('\n')
            result[modFolder]['date'] = runCommand('git log -n 1 --format=%ai @{upstream}').strip('\n')
            result[modFolder]['author'] = runCommand('git log -n 1 --format=%an @{upstream}').strip('\n')
            result[modFolder]['branch'] = runCommand('git rev-parse --abbrev-ref HEAD').strip('\n')

    os.chdir(curPath)
    return result

# function to create patches for unpublished commits and uncommitted changes in modules
def getPatches(modFolders=None):
    curPath = os.getcwd()
    result = {}

    if modFolders == None:
        modFolders = [dir for dir in os.listdir(curPath) if os.path.isdir(dir)]

    for modFolder in modFolders:
        os.chdir(os.path.join(curPath, modFolder))

        # make sure this is the top level of a git repository
        if not os.path.exists( os.path.join(os.getcwd(), '.git') ):
            printNoGitRepoWarning()
        else:
            # checkUntrackedFiles(modFolder)

            unpubPatch = runCommand('git format-patch --stdout @{upstream}')
            unCommPatch = runCommand('git diff')
            if unpubPatch != '' or unCommPatch != '': result[modFolder] = {}
            if unpubPatch != '': result[modFolder]['unpublished'] = unpubPatch
            if unCommPatch != '': result[modFolder]['uncommitted'] = unCommPatch

    os.chdir(curPath)
    return result
