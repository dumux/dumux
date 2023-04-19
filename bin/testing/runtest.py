#!/usr/bin/env python3
# SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later


"""
Helper script to run tests in DuMux and enable regression tests
by data and vtu comparisons
"""

import argparse
import shlex
import os
import sys
import subprocess
import json


try:
    from fieldcompare import FieldDataComparator, protocols, DefaultFieldComparisonCallback
    from fieldcompare.mesh import MeshFieldsComparator
    from fieldcompare.predicates import DefaultEquality, ScaledTolerance
    from fieldcompare.io import read_as

    # pylint: disable=too-many-arguments
    def makePredicateSelector(
        relThreshold,
        absThreshold,
        zeroValueThreshold,
        sourceFieldNameTransform=lambda name: name,
    ):
        """Create a predicate selector for fieldcompare emulates the Dumux behaviour"""

        def _selector(sourceField: protocols.Field, _: protocols.Field) -> protocols.Predicate:
            sourceFieldName = sourceFieldNameTransform(sourceField.name)
            absTol = zeroValueThreshold.get(
                sourceFieldName,
                ScaledTolerance(base_tolerance=absThreshold),
            )
            return DefaultEquality(abs_tol=absTol, rel_tol=relThreshold)

        return _selector

    def fieldcompareMeshData(
        source,
        ref,
        absThreshold=0.0,
        relThreshold=1e-7,
        zeroValueThreshold=None,
        ignoreFields=None,
    ):
        """Compares mesh data with the fieldcompare library"""

        print(f"-- Comparing {source} and {ref}")

        zeroValueThreshold = zeroValueThreshold or {}
        if zeroValueThreshold:
            print(f"-- Using the following absolute thresholds: {zeroValueThreshold}")

        # read the files
        sourceFields = read_as("mesh", source)
        referenceFields = read_as("mesh", ref)

        # hard-code some values for the mesh comparisons (as for Dumux legacy backend)
        sourceFields.domain.set_tolerances(abs_tol=ScaledTolerance(1e-6), rel_tol=1.5e-7)
        referenceFields.domain.set_tolerances(abs_tol=ScaledTolerance(1e-6), rel_tol=1.5e-7)

        ignoreFields = ignoreFields or []
        compare = MeshFieldsComparator(
            source=sourceFields,
            reference=referenceFields,
            field_exclusion_filter=lambda name: name in ignoreFields,
        )
        result = compare(
            predicate_selector=makePredicateSelector(
                relThreshold=relThreshold,
                absThreshold=absThreshold,
                zeroValueThreshold=zeroValueThreshold,
            ),
            fieldcomp_callback=DefaultFieldComparisonCallback(verbosity=1),
            reordering_callback=lambda msg: print(f"-- {msg}"),
        )

        print(f"-- Summary: {result.status} ({result.report})\n")

        if not result:
            return 1
        return 0

    def fieldcompareCSVData(
        source,
        ref,
        delimiter,
        absThreshold=0.0,
        relThreshold=1e-7,
        zeroValueThreshold=None,
        ignoreFields=None,
    ):
        """Compares CSV data with the fieldcompare library"""

        print(f"-- Comparing {source} and {ref}")

        zeroValueThreshold = zeroValueThreshold or {}
        if zeroValueThreshold:
            print(f"-- Using the following absolute thresholds: {zeroValueThreshold}")

        sourceFields = read_as("dsv", source, delimiter=delimiter, use_names=False)
        referenceFields = read_as("dsv", ref, delimiter=delimiter, use_names=False)

        ignoreFields = ignoreFields or []
        compare = FieldDataComparator(
            source=sourceFields,
            reference=referenceFields,
            field_exclusion_filter=lambda name: name in ignoreFields,
        )
        result = compare(
            predicate_selector=makePredicateSelector(
                relThreshold=relThreshold,
                absThreshold=absThreshold,
                zeroValueThreshold=zeroValueThreshold,
                sourceFieldNameTransform=lambda name: f"row {float(name.strip('field_'))}",
            ),
            fieldcomp_callback=DefaultFieldComparisonCallback(verbosity=1),
        )

        print(f"-- Summary: {result.status} ({result.report})\n")

        if not result:
            return 1
        return 0

    BACKEND = "fieldcompare"


# fall back to Dumux legacy backend if we don't have fieldcompare
except ImportError:
    from fuzzycomparevtu import compareVTK as fieldcompareMeshData
    from fuzzycomparedata import compareData as fieldcompareCSVData

    BACKEND = "legacy"


def readCmdParameters():
    """Read the command line parameters"""
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-c",
        "--command",
        nargs=1,
        help="The executable and optional arguments as a single string",
        required=True,
    )
    parser.add_argument(
        "-s",
        "--script",
        nargs=1,
        help=(
            "The comparison script. [fuzzy, fuzzyData, exact, <path_to_script>]"
            " where the script takes two files as arguments."
        ),
    )
    parser.add_argument(
        "-f",
        "--files",
        nargs="+",
        help=(
            "Pairs of file names (first reference, then current). "
            "Usage: '[-f ref1 cur1 [[ref2] [cur2] ...]]'"
        ),
    )
    parser.add_argument(
        "-d", "--delimiter", type=str, default=",", help="Column delimiter for data files"
    )
    parser.add_argument(
        "-r",
        "--relative",
        type=float,
        default=1e-2,
        help="maximum relative error (default=1e-2) when using fuzzy comparison",
    )
    parser.add_argument(
        "-a",
        "--absolute",
        type=float,
        default=1.5e-7,
        help="maximum absolute error (default=1.5e-7) when using fuzzy comparison",
    )
    parser.add_argument(
        "-z",
        "--zeroThreshold",
        type=json.loads,
        default="{}",
        help=(
            "Thresholds for treating numbers as zero for"
            ' a parameter as a python dict e.g. {"vel":1e-7,"delP":1.0}'
        ),
    )
    parser.add_argument(
        "-i",
        "--ignore",
        nargs="+",
        help=("Space separated list of fields to ignore in the comparison"),
    )
    args = vars(parser.parse_args())

    # check parameters
    if args["script"]:
        if len(args["files"]) % 2 != 0 or not args["files"]:
            sys.stderr.write(
                "The files have to be pairs of reference and current solution files."
                " Usage '-f [ref1] [cur1] [[ref2] [cur2] ...]'"
            )
            parser.print_help()
            sys.exit(1)
        for i in range(0, len(args["files"]) // 2):
            # delete the vtu files to compare
            referenceDirectory = (
                os.path.dirname(os.path.abspath(__file__)).rstrip("bin") + "test/references"
            )
            if os.path.dirname(args["files"][(i * 2) + 1]) == referenceDirectory:
                sys.stderr.write(
                    "Tried to delete a reference solution. "
                    "Specify reference file first, then the current solution. "
                    "Usage: '[-f ref1 cur1 [[ref2] [cur2] ...]]'"
                )
                sys.exit(1)
            subprocess.call(["rm", "-fv", args["files"][(i * 2) + 1]])

    return args


def _exactComparison(args):
    """Exact comparison driver"""
    returnCode = 0
    for i in range(0, len(args["files"]) // 2):
        print("\nExact comparison...")
        result = subprocess.call(["diff", args["files"][i * 2], args["files"][(i * 2) + 1]])
        if result:
            returnCode = 1
    return returnCode


def _fuzzyMeshComparison(args):
    """Fuzzy mesh comparison driver"""
    numFailed = 0
    for i in range(0, len(args["files"]) // 2):
        print(f"\nFuzzy data comparison with {BACKEND} backend")
        source, ref = args["files"][i * 2], args["files"][(i * 2) + 1]
        if "reference" in source and "reference" not in ref:
            source, ref = ref, source
        relThreshold = args["relative"]
        absThreshold = args["absolute"]
        zeroValueThreshold = args["zeroThreshold"]
        ignoreFields = args["ignore"]
        numFailed += fieldcompareMeshData(
            source,
            ref,
            absThreshold,
            relThreshold,
            zeroValueThreshold,
            ignoreFields,
        )

    return int(numFailed > 0)


def _fuzzyDataComparison(args):
    """Fuzzy data comparison driver"""
    numFailed = 0
    for i in range(0, len(args["files"]) // 2):
        print(f"\nFuzzy data comparison with {BACKEND} backend")
        source, ref = args["files"][i * 2], args["files"][(i * 2) + 1]
        if "reference" in source and "reference" not in ref:
            source, ref = ref, source
        delimiter = args["delimiter"]
        relThreshold = args["relative"]
        absThreshold = args["absolute"]
        zeroValueThreshold = args["zeroThreshold"]
        ignoreFields = args["ignore"]
        numFailed += fieldcompareCSVData(
            source,
            ref,
            delimiter,
            absThreshold,
            relThreshold,
            zeroValueThreshold,
            ignoreFields,
        )

    return int(numFailed > 0)


def _scriptComparison(args):
    """Script comparison driver"""
    returnCode = 0
    for i in range(0, len(args["files"]) // 2):
        print(f"\n{args['script']} comparison")
        result = subprocess.call(args["script"], args["files"][i * 2], args["files"][(i * 2) + 1])
        if result:
            returnCode = 1
    return returnCode


def runRegressionTest(args):
    """Run regression test scripts against reference data"""

    # exact comparison?
    if args["script"] == ["exact"]:
        sys.exit(_exactComparison(args))

    # fuzzy mesh comparison?
    elif args["script"] == ["fuzzy"] or args["script"] == [
        os.path.dirname(os.path.abspath(__file__)) + "/fuzzycomparevtu.py"
    ]:
        sys.exit(_fuzzyMeshComparison(args))

    # fuzzy comparison of CSV-like data sets?
    elif args["script"] == ["fuzzyData"]:
        sys.exit(_fuzzyDataComparison(args))

    # other script?
    else:
        sys.exit(_scriptComparison(args))


def runTest():
    """DuMux test driver"""

    args = readCmdParameters()

    # run the test
    res = 1
    try:
        res = subprocess.call(shlex.split(args["command"][0]))
    except OSError:
        print(args["command"][0].split())
        print("OSError: Command not found. Most likely the executable specified doesn't exist.")
        sys.exit(1)
    if res:
        sys.exit(res)

    # run the comparison
    if args["script"]:
        runRegressionTest(args=args)

    # everything is fine
    sys.exit(0)


if __name__ == "__main__":
    runTest()
