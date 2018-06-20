#!/usr/bin/env python

"""
Create files CMakeLists.txt for the given folder and all subfolders,
including the add_subdirectory(...) and install(...) commands.
Defaults to the folder `dumux` that contains the header files,
if no folder was specified.
"""

import os
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('folder', type=str, nargs='?', help='the folder to create CMakeLists.txt\'s for', default=None)
args = vars(parser.parse_args())

# default to the dumux folder (relative path to the location of this script)
if args['folder'] is None:
    rootDir = os.path.dirname(os.path.abspath(__file__)) + "/../../dumux"
else:
    rootDir = args['folder']

for folderName, subFolders, files in os.walk(rootDir):
    subFolders = sorted(subFolders)
    files = sorted(files)

    with open(folderName + "/CMakeLists.txt", "w") as cmakelists:

        for subFolder in subFolders:
            cmakelists.write("add_subdirectory({})\n".format(subFolder))

        headersExist = False
        for fileName in files:
            if fileName != "CMakeLists.txt":
                headersExist = True
                break

        if headersExist:
            if subFolders:
                cmakelists.write("\n")

            cmakelists.write("install(FILES\n")

            for fileName in files:
                if fileName != "CMakeLists.txt":
                    cmakelists.write("{}\n".format(fileName))

            cmakelists.write("DESTINATION ${{CMAKE_INSTALL_INCLUDEDIR}}/dumux/{})\n".format(folderName.replace(rootDir + '/', '')))
