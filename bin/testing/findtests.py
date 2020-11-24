#!/usr/bin/env python3

"""
Find those tests that are affected by changes
Run this in the build directory
Warning: This runs make clean on the build directory
"""

import multiprocessing as mp
import json
import subprocess
from glob import glob
from subprocess import PIPE
import os

import sys

# Whether the two lists a and b have a common member
def has_common_member(myset, mylist):
    return not myset.isdisjoint(mylist)

def make_dryrun(config):
    lines = subprocess.check_output(["make", "--dry-run", config["target"]], encoding='ascii').splitlines()
    return [l for l in lines if "g++" in l]

def build_command_and_dir(config, cache):
    lines = make_dryrun(config)
    if len(lines) == 0:
        with open(cache) as c:
            data = json.load(c)
            return data["command"], data["dir"]
    else:
        (_, dir), command = [l.split() for l in lines[0].split("&&")]
        with open(cache, "w") as c:
            json.dump({"command":command, "dir":dir}, c)

    return command, dir

# find the changes files
changed_files = subprocess.check_output(["git", "diff-tree", "-r", "--name-only", "HEAD",  "master"], encoding='ascii').splitlines()
changed_files = set(changed_files)

def is_affected_tests(test):
    with open(test) as config_file:
        config = json.load(config_file)

    command, dir = build_command_and_dir(config, "TestTargets/"+config["target"]+".json")
    # for some reason g++ writes to stderr
    lines = subprocess.run(command + ["-MM", "-H"], stderr=PIPE, stdout=PIPE, cwd=dir, encoding='ascii').stderr.splitlines()

    # filter lines
    project_dir = os.path.abspath(os.getcwd().rstrip("build-cmake"))
    test_files = set([os.path.relpath(l.lstrip(". "), project_dir) for l in lines if project_dir in l])

    if has_common_member(changed_files, test_files):
        return True, config["name"], config["target"]

    return False, config["name"], config["target"]

if __name__ == '__main__':

    subprocess.run(["make", "clean"])
    subprocess.run(["make"])

    # create cache folder
    os.makedirs("TestTargets", exist_ok=True)

    for test in glob("TestMetaData/*json"):
        print(is_affected_tests(test))
