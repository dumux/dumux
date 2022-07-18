#!/usr/bin/env python3

"""
Find the headers that are not included by our test suite
Run this in the build directory
Warning: This runs 'make clean' on the build directory
Note: This only works with compilers that support the following flags:
      - "-MM": write included user headers (no system headers) to stdout
      - "-H": Show header includes and nesting depth
"""

import json
import os
import subprocess
from argparse import ArgumentParser
from glob import glob
from findtests import getTestHeaders, buildCommandAndDir

if __name__ == "__main__":

    # parse input arguments
    parser = ArgumentParser(description="Find tests affected by changes")
    parser.add_argument(
        "-o",
        "--outfile",
        required=False,
        default="nonincluded_headers.txt",
        help="The file in which to write the non-included headers",
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

    # clean build directory
    subprocess.run(["make", "clean"], cwd=buildDir, check=False)
    subprocess.run(["make", "all"], cwd=buildDir, check=False)

    def isProjectHeader(path):
        return os.path.abspath(os.path.join(buildDir, "..")) in path

    print("Detecting headers included by the test suite:")
    allTestedHeaders = set()
    for testConfig in glob(os.path.join(buildDir, "TestMetaData") + "/*json"):
        with open(testConfig) as config:
            testConfig = json.load(config)
        cmd, dir = buildCommandAndDir(testConfig, buildTreeRoot=buildDir)
        for header in filter(lambda n: isProjectHeader(n), getTestHeaders(cmd, dir)):
            allTestedHeaders.add(header)
    print(f" ... found {len(allTestedHeaders)} headers")

    print("Detecting all headers of the project")
    allProjectHeaders = set()
    dumuxRoot = os.path.join(buildDir, "..")
    for root, _, files in os.walk(dumuxRoot):
        for f in filter(lambda n: os.path.splitext(n)[1] == ".hh", files):
            allProjectHeaders.add(os.path.abspath(
                os.path.join(os.path.join(dumuxRoot, root), f))
            )
    print(f" ... found {len(allProjectHeaders)} headers")

    def isPythonHeader(path):
        return "python" in os.path.dirname(path)

    print("Removing headers of the python bindings")
    for header in list(filter(isPythonHeader, allTestedHeaders)):
        allTestedHeaders.remove(header)
    for header in list(filter(isPythonHeader, allProjectHeaders)):
        allProjectHeaders.remove(header)

    print("The following headers are not included by the test suite:")
    nonTestedHeaders = allProjectHeaders.difference(allTestedHeaders)
    print("\n".join(str(test) for test in nonTestedHeaders))

    with open(args["outfile"], "w") as outFile:
        outFile.write("\n".join(str(test) for test in nonTestedHeaders))
