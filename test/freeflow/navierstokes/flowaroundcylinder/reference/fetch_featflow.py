#!/usr/bin/env python3
# SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

"""
Download the FeatFlow time series of the DFG cases 2D-2 and 2D-3 and extract the reference values

The archives of the FeatFlow benchmark pages
(https://wwwold.mathematik.tu-dortmund.de/~featflow/en/benchmarks/cfdbenchmarking/flow.html)
are downloaded into featflow/<dataset>/ next to this script and unpacked there (about 60 MB,
excluded from git). Archives that exist are not downloaded again (unless --force is given).

| dataset    | content                                                                         |
| :--------- | :------------------------------------------------------------------------------ |
| 2dcase2new | 2D-2, Q2/P1disc, levels 2-6, dt = 1/200, 1/400, 1/1600, t in [9, 10] s          |
|            | ("new reference results", higher-order line integration of the forces)          |
| 2dcase2old | 2D-2, Q1~/Q0 and Q2/P1disc, levels 2/3-6, dt = 1/200 to 1/1600, t in [0, 10] s  |
| 2dcase3    | 2D-3, Q2/P1disc, levels 1-6, dt = 1/1600, t in [0, 8] s                         |

The force files (bdforces_*) contain the time (column 1) and the drag and lift coefficients
(columns 3, 4), the point value files (pointvalues_*) the pressure at A = (0.15, 0.2) (column 6)
and B = (0.25, 0.2) (column 11). All values are given at the half steps of the Crank-Nicolson
method. The Q2/P1disc results of 2D-2 are taken from 2dcase2new, the Q1~/Q0 results from
2dcase2old (the old Q2/P1disc results are superseded by the new ones).

The quantities are evaluated with the same functions as the DuMux results (visualize.py) and
written to featflow.json, sorted by the discretization ("q2p1": Q2/P1disc, "q1p0": Q1~/Q0),
the time step size (coarse to fine) and the number of unknowns:

    {"<case>": {"<scheme>": {"1/<steps per second>": [{"level", "numElements", "numDofs", ...}]}}}

numDofs is the total number of scalar unknowns (velocity components and pressure), as counted
in the DuMux runs. (The tables of the FeatFlow pages list #Q2 nodes + #elements instead.)
The values of case 2D-1 are taken from the FeatFlow benchmark page. The bounds of the paper
(Schäfer & Turek, 1996) are in 1996.json.

Usage:
    python3 fetch_featflow.py
"""

import argparse
import glob
import json
import os
import re
import sys
import urllib.request
import zipfile

import numpy as np

REFERENCE_DIR = os.path.dirname(os.path.realpath(__file__))
FEATFLOW_DIR = os.path.join(REFERENCE_DIR, "featflow")
JSON_FILE = os.path.join(REFERENCE_DIR, "featflow.json")

MEDIA_URL = "https://wwwold.mathematik.tu-dortmund.de/~featflow/media"
DATASETS = {
    "2dcase2new": ["dfg_bench2new_2d/20141124_dfg_bench2_official.zip"],
    "2dcase2old": [
        "dfg_bench2_2d/draglift_q1t_cn_lv3-6_dt1-4.zip",
        "dfg_bench2_2d/draglift_q2_cn_lv3-6_dt1-4.zip",
        "dfg_bench2_2d/pressure_q1t_cn_lv3-6_dt1-4.zip",
        "dfg_bench2_2d/pressure_q2_cn_lv3-6_dt1-4.zip",
    ],
    "2dcase3": [
        "dfg_bench3_2d/draglift_q2_cn_lv1-6_dt4.zip",
        "dfg_bench3_2d/pressure_q2_cn_lv1-6_dt4.zip",
    ],
}

# force files of the evaluated runs per case and scheme (the point value files are found by
# replacing bdforces by pointvalues)
FORCE_FILES = {
    2: {
        "q2p1": os.path.join(
            "2dcase2new", "logs", "results_cc2d_bench2_referencenew", "bdforces_q2_lv*_dt*"
        ),
        "q1p0": os.path.join("2dcase2old", "bdforces_q1t_org_lv*_dt*"),
    },
    3: {"q2p1": os.path.join("2dcase3", "bdforces_lv*")},
}

# old data of case 2D-2 on the same level and time step size for both discretizations: level 6
# with dt = 1/800 s, the finest complete old Q2/P1disc run (the one with dt = 1/1600 s stops at
# t = 0.025 s). Only plotted for comparison, the old Q2/P1disc results are superseded.
OLD_SERIES = {
    "q2p1": os.path.join("2dcase2old", "bdforces_q2_lv6_dt3"),
    "q1p0": os.path.join("2dcase2old", "bdforces_q1t_org_lv6_dt3"),
}

# vertices and elements of the FeatFlow meshes per level (benchmark pages)
MESHES = {
    1: (156, 130),
    2: (572, 520),
    3: (2184, 2080),
    4: (8528, 8320),
    5: (33696, 33280),
    6: (133952, 133120),
}

# FeatFlow reference values of case 2D-1 (benchmark page, Q2/P1disc on level 6)
FEATFLOW_2D1 = {"CD": 5.57953523384, "CL": 0.010618948146, "dP": 0.11752016697}


def numDofs(scheme, level):
    """
    Total number of scalar unknowns of a FeatFlow run. The channel with the cylinder has one
    hole, so the Euler characteristic gives #edges = #vertices + #elements.
    """
    vertices, elements = MESHES[level]
    edges = vertices + elements
    if scheme == "q2p1":
        return 2 * (vertices + edges + elements) + 3 * elements
    return 2 * edges + elements


def unpack(archive, folder):
    """Unpack a zip archive into folder, refusing members outside of it"""
    root = os.path.realpath(folder)
    with zipfile.ZipFile(archive) as zipFile:
        for member in zipFile.namelist():
            target = os.path.realpath(os.path.join(root, member))
            if os.path.commonpath([root, target]) != root:
                raise RuntimeError(f"Unsafe path {member} in {archive}")
        zipFile.extractall(root)


def fetch(featflowDir, force):
    """Download the missing archives into featflow/<dataset>/ and unpack all archives there"""
    for dataset, paths in DATASETS.items():
        folder = os.path.join(featflowDir, dataset)
        os.makedirs(folder, exist_ok=True)
        for path in paths:
            archive = os.path.join(folder, os.path.basename(path))
            if force or not os.path.exists(archive):
                url = f"{MEDIA_URL}/{path}"
                print(f"Downloading {url}")
                with urllib.request.urlopen(url, timeout=60) as response:
                    data = response.read()
                with open(archive + ".part", "wb") as archiveFile:
                    archiveFile.write(data)
                os.replace(archive + ".part", archive)
            print(f"Unpacking {archive}")
            unpack(archive, folder)


def evaluation():
    """The evaluation functions of visualize.py (in the parent folder)"""
    sys.path.insert(0, os.path.dirname(REFERENCE_DIR))
    import visualize  # pylint: disable=import-outside-toplevel

    return visualize


def readTable(fileName):
    """Read a whitespace-separated FeatFlow log file"""
    with open(fileName, encoding="utf-8") as file:
        rows = [line for line in file if line.strip() and not line.startswith("#")]
    return np.loadtxt(rows, ndmin=2)


def readSeries(forcesFile):
    """Time series of the indicators (as written by the DuMux simulation) of a FeatFlow run"""
    forces = readTable(forcesFile)
    pointsFile = os.path.join(
        os.path.dirname(forcesFile),
        os.path.basename(forcesFile).replace("bdforces", "pointvalues"),
    )
    points = readTable(pointsFile)
    length = min(len(forces), len(points))
    return {
        "time": forces[:length, 1],
        "cDrag": forces[:length, 3],
        "cLift": forces[:length, 4],
        "pDiff": points[:length, 6] - points[:length, 11],
    }


def finestSeries(case, scheme="q2p1", featflowDir=FEATFLOW_DIR):
    """
    Time series of the finest evaluated run (level 6, dt = 1/1600 s) of a scheme for case 2D-2
    or 2D-3, None if the scheme is not available for the case or the data has not been downloaded
    """
    if scheme not in FORCE_FILES[case]:
        return None
    pattern = FORCE_FILES[case][scheme].replace("lv*", "lv6").replace("dt*", "dt4")
    files = glob.glob(os.path.join(featflowDir, pattern))
    return readSeries(files[0]) if files else None


def oldSeries(scheme, featflowDir=FEATFLOW_DIR):
    """Time series of the old data of case 2D-2 (level 6, dt = 1/800 s), None if not downloaded"""
    forcesFile = os.path.join(featflowDir, OLD_SERIES[scheme])
    return readSeries(forcesFile) if os.path.exists(forcesFile) else None


def stepsPerSecond(time):
    """Number of time steps per second, from the spacing of the half steps"""
    return int(round(1.0 / np.median(np.diff(time[1:]))))


def evaluateRun(case, forcesFile):
    """Quantities of one FeatFlow run, None if it cannot be evaluated"""
    series = readSeries(forcesFile)
    if case == 2:
        result = evaluation().evaluatePeriodic(series, tStart=0.0)
        if result is not None:
            del result["t0"]
    else:
        result = evaluation().evaluateTransient(series)
    if result is None:
        print(f"{forcesFile} cannot be evaluated, skipped")
        return None, None
    return stepsPerSecond(series["time"]), {
        key: float(f"{value:.10g}") for key, value in result.items()
    }


def extract(case, featflowDir):
    """Quantities of all FeatFlow runs of a case as {scheme: {dt: [run sorted by numDofs]}}"""
    references = {}
    for scheme, pattern in FORCE_FILES[case].items():
        runs = {}
        for forcesFile in glob.glob(os.path.join(featflowDir, pattern)):
            level = int(re.search(r"lv(\d+)", os.path.basename(forcesFile)).group(1))
            steps, result = evaluateRun(case, forcesFile)
            if result is not None:
                run = {"level": level, "numElements": MESHES[level][1]}
                run["numDofs"] = numDofs(scheme, level)
                runs.setdefault(steps, []).append({**run, **result})
        if not runs:
            sys.exit(f"No FeatFlow runs found for {pattern}, check the downloaded data")
        references[scheme] = {
            f"1/{steps}": sorted(runs[steps], key=lambda run: run["numDofs"])
            for steps in sorted(runs)
        }
    return references


def main():
    """Download the FeatFlow data, extract the reference values and write featflow.json"""
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument(
        "--featflow-dir", default=FEATFLOW_DIR, help="folder for the FeatFlow time series"
    )
    parser.add_argument("--force", action="store_true", help="download the archives again")
    parser.add_argument(
        "--no-download", action="store_true", help="only extract from the existing data"
    )
    args = parser.parse_args()

    if not args.no_download:
        fetch(args.featflow_dir, args.force)

    references = {
        "1": FEATFLOW_2D1,
        "2": extract(2, args.featflow_dir),
        "3": extract(3, args.featflow_dir),
    }
    with open(JSON_FILE, "w", encoding="utf-8") as jsonFile:
        json.dump(references, jsonFile, indent=4)
        jsonFile.write("\n")
    print(f"Wrote {JSON_FILE}")


if __name__ == "__main__":
    main()
