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
