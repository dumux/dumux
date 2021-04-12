#!/usr/bin/env python3

"""
Find those tests that are affected by changes
Run this in the build directory
Warning: This runs make clean on the build directory
"""

import json
import subprocess
from argparse import ArgumentParser
from glob import glob
from subprocess import PIPE
import os


# Whether the two lists a and b have a common member
def has_common_member(myset, mylist):
    return not myset.isdisjoint(mylist)


# make dry run and return the compilation command
def get_compile_command(testConfig):
    lines = subprocess.check_output(["make", "--dry-run",
                                     testConfig["target"]],
                                    encoding='ascii').splitlines()
    commands = list(filter(lambda comm: 'g++' in comm, lines))
    assert len(commands) <= 1
    return commands[0] if commands else None


# get the command and folder to compile the given test
def build_command_and_dir(testConfig, cache):
    compCommand = get_compile_command(testConfig)
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
def is_affected_test(testConfigFile, changed_files):
    with open(testConfigFile) as configFile:
        testConfig = json.load(configFile)

    cacheFile = "TestTargets/" + testConfig["target"] + ".json"
    command, dir = build_command_and_dir(testConfig, cacheFile)

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

    test_files = set([os.path.relpath(header.lstrip(". "), projectDir)
                      for header in filter(isProjectHeader, headers)])

    if has_common_member(changed_files, test_files):
        return True, testConfig["name"], testConfig["target"]

    return False, testConfig["name"], testConfig["target"]


if __name__ == '__main__':

    # parse input arguments
    parser = ArgumentParser(description='Find tests affected by changes')
    parser.add_argument('-s', '--source', required=False, default='HEAD',
                        help='The source tree (default: `HEAD`)')
    parser.add_argument('-t', '--target', required=False, default='master',
                        help='The tree to compare against (default: `master`)')
    args = vars(parser.parse_args())

    # find the changes files
    changed_files = subprocess.check_output(["git", "diff-tree",
                                             "-r", "--name-only",
                                             args['source'],  args['target']],
                                            encoding='ascii').splitlines()
    changed_files = set(changed_files)

    # clean build directory
    subprocess.run(["make", "clean"])
    subprocess.run(["make"])

    # create cache folder
    os.makedirs("TestTargets", exist_ok=True)

    # detect affected tests
    print("Detecting affected tests:")
    count = 0
    affectedTests = {}
    for test in glob("TestMetaData/*json"):
        affected, name, target = is_affected_test(test, changed_files)
        if affected:
            print("\t- {}".format(name))
            affectedTests[name] = {'target': target}
            count += 1
    print("Detected {} affected tests".format(count))

    with open('affectedtests.json', 'w') as jsonFile:
        json.dump(affectedTests, jsonFile)
