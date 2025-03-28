#!/usr/bin/env python3
# SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

""" A module for fuzzy comparing data files.

This module provides methods to compare two data files.
Applicable for all style formats like e.g. csv files.
Fuzzy compares numbers by using absolute and/or relative difference comparison.

"""
import argparse
import json
import sys
from dumux_fuzzycompare_legacy import compareData

# Note: these issues can be improved on by factoring out functions
# but we ignore it for know ("legacy code")
# pylint: disable=too-many-arguments,too-many-locals,too-many-branches,too-many-statements

# main program if called as script return appropriate error codes
if __name__ == "__main__":
    # handle arguments and print help message
    parser = argparse.ArgumentParser(
        description="Fuzzy compare of two data files (e.g csv). \
        The files are accepted if for every value the difference is below the absolute error \
        or below the relative error or below both."
    )
    parser.add_argument("data_file_1", type=str, help="first file to compare")
    parser.add_argument("data_file_2", type=str, help="second file to compare")
    parser.add_argument("delimiter", type=str, help="second file to compare")
    parser.add_argument(
        "-r", "--relative", type=float, default=1e-2, help="maximum relative error (default=1e-2)"
    )
    parser.add_argument(
        "-a",
        "--absolute",
        type=float,
        default=1.5e-7,
        help="maximum absolute error (default=1.5e-7)",
    )
    parser.add_argument("-v", "--verbose", type=bool, default=True, help="verbosity of the script")
    parser.add_argument(
        "-z",
        "--zeroThreshold",
        type=json.loads,
        default="{}",
        help=(
            "Thresholds for treating numbers as zero for a parameter as a python dict "
            'e.g. {"vel":1e-7,"delP":1.0}'
        ),
    )
    parser.add_argument(
        "-i",
        "--ignore",
        nargs="+",
        help=("Space separated list of fields to ignore in the comparison"),
    )
    args = vars(parser.parse_args())

    sys.exit(
        compareData(
            args["data_file_1"],
            args["data_file_2"],
            args["delimiter"],
            args["absolute"],
            args["relative"],
            args["zeroThreshold"],
            args["verbose"],
            args["ignore"],
        )
    )
