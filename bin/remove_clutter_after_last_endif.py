#!/usr/bin/env python3
# SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later


""""
Remove clutter after the last #endif (header guard)
in C++ header files
"""

import os


def clearAfterLastEndIf(fileName):
    """Clear a single headerfile with name fileName"""
    with open(fileName, "r") as header:
        split = header.read().split("#endif")
        split[-1] = "\n"
    with open(fileName, "w") as header:
        header.write("#endif".join(split))


def run():
    """Main driver: go through all header in directory recursively"""
    for root, _, files in os.walk(os.getcwd()):
        for file in files:
            if file.endswith(".hh"):
                clearAfterLastEndIf(os.path.join(root, file))


if __name__ == "__main__":
    run()
