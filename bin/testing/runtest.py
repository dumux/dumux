#!/usr/bin/env python3
# SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

"""
Helper script to run tests: use dumux_runtest.py instead
"""


from dumux_runtest import runTest

if __name__ == "__main__":
    print(">>>>> THIS SCRIPT IS DEPRECATED <<<<<<<<< Please use dumux_runtest.py instead.")
    runTest()
