#!/usr/bin/env python3

"""
Get the names of the files that differ between two git trees
"""

import os
import subprocess
from argparse import ArgumentParser


def getCommandOutput(command, cwd=None):
    return subprocess.check_output(command, encoding='ascii', cwd=cwd)


# get the files that differ between two trees in a git repo
def getChangedFiles(gitFolder, sourceTree, targetTree):

    gitFolder = os.path.abspath(gitFolder)
    root = getCommandOutput(
        command=['git', 'rev-parse', '--show-toplevel'],
        cwd=gitFolder
    ).strip('\n')
    changedFiles = getCommandOutput(
        command=["git", "diff-tree", "-r", "--name-only", sourceTree, targetTree],
        cwd=gitFolder
    ).splitlines()

    return [os.path.join(root, file) for file in changedFiles]


if __name__ == '__main__':

    # parse input arguments
    parser = ArgumentParser(
        description='Get the files that differ between two git-trees'
    )
    parser.add_argument('-f', '--folder',
                        required=False, default='.',
                        help='The path to a folder within the git repository')
    parser.add_argument('-s', '--source-tree',
                        required=False, default='HEAD',
                        help='The source tree (default: `HEAD`)')
    parser.add_argument('-t', '--target-tree',
                        required=False, default='master',
                        help='The tree to compare against (default: `master`)')
    parser.add_argument('-o', '--outfile',
                        required=False, default='changedfiles.txt',
                        help='The file in which to write the changed files')
    args = vars(parser.parse_args())

    changedFiles = getChangedFiles(args['folder'],
                                   args['source_tree'],
                                   args['target_tree'])

    with open(args['outfile'], 'w') as outFile:
        for file in changedFiles:
            outFile.write(f"{os.path.abspath(file)}\n")
