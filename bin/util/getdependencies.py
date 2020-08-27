#!/usr/bin/env python3

import os
import subprocess

# function to determine module name from module file
def getModuleName(modulePath):
    modFile = os.path.join(modulePath, 'dune.module')
    if not os.path.exists(modFile):
        return None

    with open(modFile, 'r') as modFile:
        for line in modFile.readlines():
            line = line.strip('\n').split(' ')
            if line[0] == 'Module:':
                return line[1]
    return None

# get the dependencies of a dune module located in the given directory
def getDependencies(modulePath):
    parentPath = os.path.join(modulePath, '../')
    duneControlPath = os.path.join(parentPath, 'dune-common/bin/dunecontrol')
    if not os.path.exists(duneControlPath):
        raise RuntimeError('Could not find dunecontrol, expected to be in {}'.format(duneControlPath))

    curPath = os.getcwd()
    os.chdir(parentPath)
    output = subprocess.run(['./dune-common/bin/dunecontrol', '--module={}'.format(getModuleName(modulePath))],
                            check=True, text=True, capture_output=True)
    for line in output.stdout.split('\n'):
        if "going to build" in line:
            line = line.replace('going to build', '').replace('---', '').replace('done', '')
            line = line.strip('\n').strip()
            line = line.split(' ')
            deps = line

    result = []
    for dir in [d for d in os.listdir(os.getcwd()) if os.path.isdir(d)]:
        depModName = getModuleName(os.path.join(os.getcwd(), dir))
        if depModName in deps:
            result.append({})
            result[-1]['name'] = depModName
            result[-1]['folder'] = dir

    if len(result) != len(deps):
        raise RuntimeError("Could not find the sub-folders of all dependencies")

    os.chdir(curPath)
    return result
