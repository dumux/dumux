# SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

""" A module for fuzzy comparing data files.

This module provides methods to compare two data files.
Applicable for all style formats like e.g. csv files.
Fuzzy compares numbers by using absolute and/or relative difference comparison.

"""
import argparse
import csv
import json
import sys
from fuzzycomparevtu import isFuzzyEqualText

# Note: these issues can be improved on by factoring out functions
# but we ignore it for know ("legacy code")
# pylint: disable=too-many-arguments,too-many-locals,too-many-branches,too-many-statements


def compareData(
    dataFile1,
    dataFile2,
    delimiter,
    absolute=1.5e-7,
    relative=1e-2,
    zeroValueThreshold=None,
    ignoreFields=None,
    verbose=True,
):
    """take two data files and compare them. Returns an exit key as returnvalue.

    Arguments:
    ----------
    dataFile1, dataFile2 : string
        The filenames of the data files to compare
    delimiter: string
        The delimiter for the columns

    Keyword Arguments:
    ------------------
    absolute : float
        The epsilon used for comparing numbers with an absolute criterion
    relative: float
        The epsilon used for comparing numbers with an relative criterion
    zeroValueThreshold: dict
        A dictionary of parameter value pairs that set the threshold under
        which a number is treated as zero for a certain parameter. Use this parameter if
        you have to avoid comparisons of very small numbers for a certain parameter.
    ignoreFields: list
        A list of field names to be ignored in the comparison
    verbose : bool
        If the script should produce informative output. Enabled by default as the details
        give the tester a lot more information on why tests fail.
    """

    if verbose:
        print(
            f"Comparing {dataFile1} and {dataFile2}\n"
            f"... with a maximum relative error of {relative} and "
            f"a maximum absolute error of {absolute}*max_abs_parameter_value."
        )

    zeroValueThreshold = zeroValueThreshold or {}

    # implement ignoring fields by setting a very high threshold
    if ignoreFields is not None:
        for field in ignoreFields:
            zeroValueThreshold[field] = 1e100

    # construct element tree from data files
    with open(dataFile1, "r") as data1:
        data1 = list(csv.reader(data1, delimiter=delimiter))
    with open(dataFile1, "r") as data2:
        data2 = list(csv.reader(data2, delimiter=delimiter))

    if len(data1) != len(data2):
        print(
            "Length of data1 and data2 not equal: ref=",
            len(data1),
            ",new=",
            len(data2),
            ". Aborting!",
        )
        sys.exit(3)

    isEqual = True
    for i in range(0, len(data1[0])):
        valueA = data1[0][i]
        valueB = data2[0][i]
        for j in range(1, len(data1)):
            valueA += f" {data1[j][i]}"
            valueB += f" {data2[j][i]}"

        if not isFuzzyEqualText(
            valueA,
            valueB,
            f"row {i}",
            len(data1),
            absolute,
            relative,
            zeroValueThreshold,
            verbose,
        ):
            if verbose:
                isEqual = False
            else:
                return False

    return 0 if isEqual else 1


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
