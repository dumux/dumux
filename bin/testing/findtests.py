#!/usr/bin/env python3

"""
Find those tests that are affected by changes
Run this in the build directory
Warning: This runs 'make clean' on the build directory
Note: This only works with compilers that support the following flags:
      - "-MM": write included user headers (no system headers) to stdout
      - "-H": Show header includes and nesting depth
"""

import json
import subprocess
from argparse import ArgumentParser
from glob import glob
from subprocess import PIPE
import os
from multiprocessing import Pool
from functools import partial


# Check if the set a contains a member of list b
def hasCommonMember(myset, mylist):
    return not myset.isdisjoint(mylist)


# make dry run and return the compilation command
def getCompileCommand(testConfig):
    lines = subprocess.check_output(["make", "--dry-run",
                                     testConfig["target"]],
                                    encoding='ascii').splitlines()

    def hasCppCommand(line):
        return any(cpp in line for cpp in ['g++', 'clang++'])

    commands = list(filter(lambda line: hasCppCommand(line), lines))
    assert len(commands) <= 1
    return commands[0] if commands else None


# get the command and folder to compile the given test
def buildCommandAndDir(testConfig, cache):
    compCommand = getCompileCommand(testConfig)
    if compCommand is None:
        with open(cache) as c:
            data = json.load(c)
            return data["command"], data["dir"]
    else:
        (_, dir), command = [comm.split() for comm in compCommand.split("&&")]
        with open(cache, "w") as c:
            json.dump({"command": command, "dir": dir}, c)
        return command, dir


# check if a test is affected by changes in the given files
def isAffectedTest(testConfigFile, changedFiles):
    with open(testConfigFile) as configFile:
        testConfig = json.load(configFile)

    cacheFile = "TestTargets/" + testConfig["target"] + ".json"
    command, dir = buildCommandAndDir(testConfig, cacheFile)
    mainFile = command[-1]

    # detect headers included in this test
    # -MM skips headers from system directories
    # -H  prints the name(+path) of each used header
    # for some reason g++ writes to stderr
    headers = subprocess.run(command + ["-MM", "-H"],
                             stderr=PIPE, stdout=PIPE, cwd=dir,
                             encoding='ascii').stderr.splitlines()

    # filter only headers from this project and turn them into relative paths
    projectDir = os.path.abspath(os.getcwd().rstrip("build-cmake"))

    def isProjectHeader(headerPath):
        return projectDir in headerPath

    testFiles = [os.path.relpath(mainFile.lstrip(". "), projectDir)]
    testFiles.extend([os.path.relpath(header.lstrip(". "), projectDir)
                      for header in filter(isProjectHeader, headers)])
    testFiles = set(testFiles)

    if hasCommonMember(changedFiles, testFiles):
        return True, testConfig["name"], testConfig["target"]

    return False, testConfig["name"], testConfig["target"]


if __name__ == '__main__':

    # parse input arguments
    parser = ArgumentParser(description='Find tests affected by changes')
    parser.add_argument('-s', '--source',
                        required=False, default='HEAD',
                        help='The source tree (default: `HEAD`)')
    parser.add_argument('-t', '--target',
                        required=False, default='master',
                        help='The tree to compare against (default: `master`)')
    parser.add_argument('-np', '--num-processes',
                        required=False, type=int, default=4,
                        help='Number of processes (default: 4)')
    parser.add_argument('-f', '--outfile',
                        required=False, default='affectedtests.json',
                        help='The file in which to write the affected tests')
    args = vars(parser.parse_args())

    # find the changes files
    changedFiles = subprocess.check_output(
        ["git", "diff-tree", "-r", "--name-only", args['source'],  args['target']],
        encoding='ascii'
    ).splitlines()
    changedFiles = set(changedFiles)

    # clean build directory
    subprocess.run(["make", "clean"])
    subprocess.run(["make"])

    # create cache folder
    os.makedirs("TestTargets", exist_ok=True)

    # detect affected tests
    print("Detecting affected tests:")
    affectedTests = {}
    tests = glob("TestMetaData/*json")

    numProcesses = max(1, args['num_processes'])
    findAffectedTest = partial(isAffectedTest, changedFiles=changedFiles)
    with Pool(processes=numProcesses) as p:
        for affected, name, target in p.imap_unordered(findAffectedTest, tests, chunksize=4):
            if affected:
                affectedTests[name] = {'target': target}
                print('\t- {} (target: {})'.format(name, target))

    print("Detected {} affected tests".format(len(affectedTests)))

    with open(args['outfile'], 'w') as jsonFile:
        json.dump(affectedTests, jsonFile)
