#!/usr/bin/env python3
# SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later


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


def hasCommonMember(myset, mylist):
    """Check if the set a contains a member of list b"""
    return not myset.isdisjoint(mylist)


def getCompileCommand(testConfig, buildTreeRoot="."):
    """make dry run and return the compilation command"""
    target = testConfig["target"]
    lines = subprocess.check_output(
        ["make", "-B", "--dry-run", target], encoding="ascii", cwd=buildTreeRoot
    ).splitlines()

    def hasCppCommand(line):
        return any(cpp in line for cpp in ["g++", "clang++"])

    # there may be library build commands first, last one is the actual target
    commands = list(filter(hasCppCommand, lines))
    return commands[-1] if commands else None


def buildCommandAndDir(testConfig, buildTreeRoot="."):
    """get the command and folder to compile the given test"""
    compCommand = getCompileCommand(testConfig, buildTreeRoot)
    if compCommand is None:
        raise Exception(f"Could not determine compile command for {testConfig}")

    (_, directory), command = [comm.split() for comm in compCommand.split("&&")]
    return command, directory


def isAffectedTest(testConfigFile, changedFiles, buildTreeRoot="."):
    """check if a test is affected by changes in the given files"""
    with open(testConfigFile) as configFile:
        testConfig = json.load(configFile)

    # first we check if the CMakeLists.txt file of the test is changed
    # then we always mark the test as affected
    if f"{testConfig['source_dir']}/CMakeLists.txt" in changedFiles:
        return True, testConfig["name"], testConfig["target"]

    # next we use the compiler to detect changed header files
    command, directory = buildCommandAndDir(testConfig, buildTreeRoot)
    mainFile = command[-1]

    # detect headers included in this test
    # -MM skips headers from system directories
    # -H  prints the name(+path) of each used header
    # for some reason g++ writes to stderr
    headers = subprocess.run(
        command + ["-MM", "-H"],
        stderr=PIPE,
        stdout=PIPE,
        cwd=directory,
        encoding="ascii",
        check=False,
    ).stderr.splitlines()
    headers = [h.lstrip(". ") for h in headers]
    headers.append(mainFile)

    if hasCommonMember(changedFiles, headers):
        return True, testConfig["name"], testConfig["target"]

    return False, testConfig["name"], testConfig["target"]


if __name__ == "__main__":

    # parse input arguments
    parser = ArgumentParser(description="Find tests affected by changes")
    parser.add_argument(
        "-l", "--file-list", required=True, help="A file containing a list of files that changed"
    )
    parser.add_argument(
        "-np",
        "--num-processes",
        required=False,
        type=int,
        default=8,
        help="Number of processes (default: 8)",
    )
    parser.add_argument(
        "-o",
        "--outfile",
        required=False,
        default="affectedtests.json",
        help="The file in which to write the affected tests",
    )
    parser.add_argument(
        "-b",
        "--build-dir",
        required=False,
        default=".",
        help="The path to the top-level build directory of the project to be checked",
    )
    args = vars(parser.parse_args())

    buildDir = os.path.abspath(args["build_dir"])
    targetFile = os.path.abspath(args["outfile"])
    with open(args["file_list"]) as files:
        changedFileList = set(line.strip("\n") for line in files.readlines())

    # clean build directory
    subprocess.run(["make", "clean"], cwd=buildDir, check=False)
    subprocess.run(["make", "all"], cwd=buildDir, check=False)

    # detect affected tests
    print("Detecting affected tests:")
    affectedTests = {}
    tests = glob(os.path.join(buildDir, "TestMetaData") + "/*json")

    numProcesses = max(1, args["num_processes"])
    findAffectedTest = partial(isAffectedTest, changedFiles=changedFileList, buildTreeRoot=buildDir)
    with Pool(processes=numProcesses) as pool:
        for affected, name, cmakeTarget in pool.imap_unordered(
            findAffectedTest, tests, chunksize=4
        ):
            if affected:
                affectedTests[name] = {"target": cmakeTarget}
                print(f"\t- {name} (target: {cmakeTarget})")

    print(f"Detected {len(affectedTests)} affected tests")

    with open(targetFile, "w") as jsonFile:
        json.dump(affectedTests, jsonFile)
