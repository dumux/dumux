#!/usr/bin/env python3
# SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

""" A module for fuzzy comparing VTK files.

This module provides methods to compare two VTK files. Applicable
for all VTK style formats like VTK files. Fuzzy compares numbers by
using absolute and/or relative difference comparison.

"""
import argparse
import json
import sys

from dumux_fuzzycompare_legacy import compareVTK

# Note: these issues can be improved on by factoring out functions
# but we ignore it for know ("legacy code")
# pylint: disable=too-many-arguments,too-many-locals,too-many-branches,too-many-statements

# main program if called as script return appropriate error codes
if __name__ == "__main__":
    # handle arguments and print help message
    parser = argparse.ArgumentParser(
        description="Fuzzy compare of two VTK\
        (Visualization Toolkit) files. The files are accepted if for every\
        value the difference is below the absolute error or below the\
        relative error or below both.  If a pvd file is given instead, the\
        corresponding possibly parallel vtk file(s) have to be present and\
        will be converted to a (series of) sequential vtk file(s). The last\
        one in the natural ordering of these files will be taken for\
        comparison."
    )
    parser.add_argument("vtk_file_1", type=str, help="first file to compare")
    parser.add_argument("vtk_file_2", type=str, help="second file to compare")
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
    parser.add_argument(
        "-z",
        "--zeroThreshold",
        type=json.loads,
        default="{}",
        help=(
            "Thresholds for treating numbers as zero for a parameter as a python dict"
            'e.g. {"vel":1e-7,"delP":1.0}'
        ),
    )
    parser.add_argument(
        "-i",
        "--ignore",
        nargs="+",
        help=("Space separated list of fields to ignore in the comparison"),
    )
    parser.add_argument("-v", "--verbose", dest="verbose", action="store_true")
    parser.add_argument("--no-verbose", dest="verbose", action="store_false")
    parser.set_defaults(verbose=True)
    args = vars(parser.parse_args())

    sys.exit(
        compareVTK(
            args["vtk_file_1"],
            args["vtk_file_2"],
            args["absolute"],
            args["relative"],
            args["zeroThreshold"],
            args["verbose"],
            args["ignore"],
        )
    )
