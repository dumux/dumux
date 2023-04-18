#!/usr/bin/env python3
# SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later


"""
Create files CMakeLists.txt for the given folder and all subfolders,
including the add_subdirectory(...) and install(...) commands.
Defaults to the folder `dumux` that contains the header files,
if no folder was specified.
"""

import os
import argparse

SPDX_HEADER = """\
# SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

"""


def createCMakeLists():
    """Create the CMakeLists.txt files"""

    parser = argparse.ArgumentParser()
    parser.add_argument(
        "folder",
        type=str,
        nargs="?",
        help="the folder to create CMakeLists.txt's for",
        default=None,
    )
    args = vars(parser.parse_args())

    # default to the dumux folder (relative path to the location of this script)
    if args["folder"] is None:
        rootDir = os.path.dirname(os.path.abspath(__file__)) + "/../dumux"
    else:
        rootDir = args["folder"]

    ignoreFolders = ["", "io/format/fmt", "io/xml"]
    extensions = [".hh", ".inc"]
    for fullFolderName, subFolders, files in os.walk(rootDir):
        # alphabetically sort
        subFolders = sorted(subFolders)
        files = sorted(files)
        # get folder name relative to dumux
        folderName = fullFolderName.replace(rootDir + "/", "").replace(rootDir, "")
        if folderName not in ignoreFolders:
            with open(fullFolderName + "/CMakeLists.txt", "w") as cmakeLists:
                cmakeLists.write(SPDX_HEADER)
                # add subfolders
                for subFolder in subFolders:
                    cmakeLists.write(f"add_subdirectory({subFolder})\n")

                headersExist = False
                for fileName in files:
                    ext = os.path.splitext(fileName)[1]
                    if ext in extensions:
                        headersExist = True
                        break

                if headersExist:
                    if subFolders:
                        cmakeLists.write("\n")
                    # collect all files to be installed in a CMake variable
                    headerGuard = "DUMUX_" + folderName.upper().replace("/", "_") + "_HEADERS"
                    cmakeLists.write(f"file(GLOB {headerGuard}{' *'.join([''] + extensions)})\n")
                    cmakeLists.write(f"install(FILES ${{{headerGuard}}}\n")
                    cmakeLists.write(
                        f"        DESTINATION ${{CMAKE_INSTALL_INCLUDEDIR}}/dumux/{folderName})\n"
                    )


if __name__ == "__main__":
    createCMakeLists()
