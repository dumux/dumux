#!/usr/bin/env python3
# SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

"""
DFG flow around a cylinder benchmark (Schäfer & Turek, 1996): grid convergence of case 2D-1

Runs the stationary case 2D-1 (Re = 20) for all discretization schemes on a sequence of
uniformly refined grids. The coarse meshes meshes/cylinder_channel_coarse.msh (triangles) and
meshes/cylinder_channel_quad_coarse.msh (quadrilaterals) are refined with Grid.Refinement.
For comparison, each scheme is also run once on the default mesh of the regression tests
(meshes/cylinder_channel.msh or meshes/cylinder_channel_quad.msh, run <scheme>_default), which
resolves the cylinder twice as fine as the coarse meshes.
The advective term uses central differences (Flux.UpwindWeight = 0.5).

Each run gets its own folder runs/2d1/<name>/ (<name> = <scheme>_r<level> or <scheme>_default)
with the log, the benchmark indicators (<name>.json) and the VTK output (vtk/). A failed run is
reported and skipped.
The errors and the experimental orders of convergence are written to results/convergence_2d1.md,
and visualize.py plots them.

Usage (in the build directory test/freeflow/navierstokes/flow_around_cylinder):
    python3 benchmark.py
    python3 benchmark.py --schemes pq1bubble_box --levels 0 1 2
    python3 benchmark.py -- -Newton.MaxRelativeShift 1e-10   # further simulation parameters
    python3 visualize.py
"""

import argparse
import json
import math
import os
import resource
import subprocess

SCRIPT_DIR = os.path.dirname(os.path.realpath(__file__))
RUNS_DIR = os.path.join("runs", "2d1")
RESULTS_DIR = os.path.join(SCRIPT_DIR, "results")

# scheme: (executable, mesh)
SCHEMES = {
    "diamond_tpfa": (
        "test_ff_navierstokes_dfg_benchmark_stationary_diamond_tpfa",
        "meshes/cylinder_channel_coarse.msh",
    ),
    "pq1bubble_box": (
        "test_ff_navierstokes_dfg_benchmark_stationary_pq1bubble_box",
        "meshes/cylinder_channel_coarse.msh",
    ),
    "pq1bubble_diamond": (
        "test_ff_navierstokes_dfg_benchmark_stationary_pq1bubble_diamond",
        "meshes/cylinder_channel_quad_coarse.msh",
    ),
}

# default meshes of the regression tests (mesh size 0.01 in the channel, 0.002 on the cylinder)
DEFAULT_MESHES = {
    "diamond_tpfa": "meshes/cylinder_channel.msh",
    "pq1bubble_box": "meshes/cylinder_channel.msh",
    "pq1bubble_diamond": "meshes/cylinder_channel_quad.msh",
}

# FeatFlow reference values of case 2D-1
REFERENCE = {"cDrag": 5.57953523384, "cLift": 0.010618948146, "pDiff": 0.11752016697}


def runName(scheme, level=None):
    """Name of a run (folder, problem name and JSON file), level None for the default mesh"""
    return f"{scheme}_default" if level is None else f"{scheme}_r{level}"


def resultFile(name):
    """Path of the JSON file with the benchmark indicators of a run"""
    return os.path.join(RUNS_DIR, name, f"{name}.json")


def limitMemory(gigabytes):
    """Return a function limiting the address space of the child process"""

    def setLimit():
        limit = int(gigabytes * 1024**3)
        hardLimit = resource.getrlimit(resource.RLIMIT_AS)[1]
        if hardLimit != resource.RLIM_INFINITY:
            limit = min(limit, hardLimit)
        resource.setrlimit(resource.RLIMIT_AS, (limit, hardLimit))

    return setLimit if gigabytes > 0 else None


def runSimulation(scheme, level, args, extraArgs):
    """
    Run case 2D-1 for one scheme on the given level of the coarse mesh (level None: default
    mesh of the regression tests), return True on success
    """
    executable, mesh = SCHEMES[scheme]
    if level is None:
        mesh = DEFAULT_MESHES[scheme]
    name = runName(scheme, level)
    folder = os.path.join(RUNS_DIR, name)
    if os.path.exists(resultFile(name)) and not args.overwrite:
        print(f"Skipping {name} (results exist, use --overwrite to rerun)")
        return True

    vtkFolder = os.path.join(folder, "vtk")
    os.makedirs(vtkFolder, exist_ok=True)
    cmd = [
        os.path.abspath(executable),
        os.path.abspath("params_1.input"),
        "-Grid.File",
        os.path.abspath(mesh),
        "-Grid.Refinement",
        str(level or 0),
        "-Flux.UpwindWeight",
        "0.5",
        "-Problem.EnableInertiaTerms",
        "true",
        "-Problem.Name",
        name,
        "-Problem.IndicatorFile",
        os.path.join("..", name),
    ] + extraArgs

    print(f"Running {name}: " + " ".join(cmd))
    with open(os.path.join(folder, f"{name}.log"), "w", encoding="utf-8") as logFile:
        result = subprocess.run(
            cmd,
            cwd=vtkFolder,
            stdout=logFile,
            stderr=subprocess.STDOUT,
            preexec_fn=limitMemory(args.memory_limit),
            check=False,
        )
    if result.returncode != 0:
        print(f"  {name} failed (exit code {result.returncode}), see {logFile.name}")
        return False
    return True


def readResult(name):
    """The benchmark indicators of a finished run, None if the run does not exist or failed"""
    fileName = resultFile(name)
    if not os.path.exists(fileName):
        return None
    with open(fileName, encoding="utf-8") as jsonFile:
        return json.load(jsonFile)


def readResults(schemes, levels):
    """
    The benchmark indicators of all finished runs as {scheme: ([result per level], default)},
    where default is the result on the default mesh (or None)
    """
    results = {}
    for scheme in schemes:
        runs = [readResult(runName(scheme, level)) for level in levels]
        results[scheme] = ([run for run in runs if run], readResult(runName(scheme)))
    return results


def orderOfConvergence(errorCoarse, errorFine, dofsCoarse, dofsFine):
    """Experimental order of convergence with respect to h ~ N^(-1/2) in 2D"""
    if errorCoarse <= 0.0 or errorFine <= 0.0:
        return float("nan")
    return -2.0 * math.log(errorFine / errorCoarse) / math.log(dofsFine / dofsCoarse)


def tableRow(run, previous, label):
    """Table row of a run with the errors and the orders of convergence (w.r.t. the previous run)"""
    row = [label, str(run["numElements"]), str(run["numDofs"]["total"])]
    for key, reference in REFERENCE.items():
        error = abs(run[key] - reference)
        eoc = ""
        if previous is not None:
            previousError = abs(previous[key] - reference)
            dofs = (previous["numDofs"]["total"], run["numDofs"]["total"])
            eoc = f"{orderOfConvergence(previousError, error, *dofs):.2f}"
        row += [f"{run[key]:.8g}", f"{error:.3e}", eoc]
    return "| " + " | ".join(row) + " |"


def writeTable(results, fileName):
    """Write the errors and the orders of convergence of all schemes as Markdown tables"""
    references = ", ".join(f"{key} = {value}" for key, value in REFERENCE.items())
    lines = [
        "# Grid convergence of case 2D-1 (Re = 20, Flux.UpwindWeight = 0.5)",
        "",
        f"Errors with respect to the FeatFlow reference values {references}.",
        "The experimental order of convergence (EOC) is computed with respect to "
        "h ~ N^(-1/2), where N is the total number of unknowns.",
    ]
    for scheme, (runs, default) in results.items():
        lines += [
            "",
            f"## {scheme}",
            "",
            "| level | elements | unknowns | C_D | error C_D | EOC | C_L | error C_L | EOC "
            "| Δp | error Δp | EOC |",
            "| ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |",
        ]
        lines += [
            tableRow(run, runs[i - 1] if i > 0 else None, str(run["refinement"]))
            for i, run in enumerate(runs)
        ]
        if default is not None:
            lines.append(tableRow(default, None, "default mesh"))
    os.makedirs(os.path.dirname(fileName), exist_ok=True)
    with open(fileName, "w", encoding="utf-8") as tableFile:
        tableFile.write("\n".join(lines) + "\n")
    print(f"Wrote {fileName}")


def main():
    """Run the grid convergence study of case 2D-1"""
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument(
        "--schemes", nargs="+", choices=list(SCHEMES), default=list(SCHEMES), help="schemes"
    )
    parser.add_argument(
        "--levels", nargs="+", type=int, default=[0, 1, 2, 3], help="values of Grid.Refinement"
    )
    parser.add_argument(
        "--memory-limit",
        type=float,
        default=6.0,
        help="address space limit per simulation in GB (0: no limit)",
    )
    parser.add_argument("--overwrite", action="store_true", help="rerun existing runs")
    parser.add_argument(
        "--results-dir", default=RESULTS_DIR, help="output folder of the convergence table"
    )
    args, extraArgs = parser.parse_known_args()
    if extraArgs and extraArgs[0] == "--":
        extraArgs = extraArgs[1:]

    for scheme in args.schemes:
        executable = SCHEMES[scheme][0]
        if not os.path.exists(executable):
            print(f"Executable {executable} not found, build it first. Skipping {scheme}.")
            continue
        for level in args.levels + [None]:
            runSimulation(scheme, level, args, extraArgs)

    writeTable(
        readResults(args.schemes, args.levels),
        os.path.join(args.results_dir, "convergence_2d1.md"),
    )


if __name__ == "__main__":
    main()
