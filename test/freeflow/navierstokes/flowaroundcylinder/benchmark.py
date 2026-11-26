#!/usr/bin/env python3
# SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

"""
DFG flow around a cylinder benchmark (Schäfer & Turek, 1996): simulation runs

Case 2D-1 (--case 1, default): grid convergence study of the stationary case (Re = 20) for all
discretization schemes. The coarse meshes meshes/cylinder_channel_coarse.msh (triangles) and
meshes/cylinder_channel_quad_coarse.msh (quadrilaterals) are refined with Grid.Refinement.
For comparison, each scheme is also run once on the default mesh of the regression tests
(meshes/cylinder_channel.msh or meshes/cylinder_channel_quad.msh, run <scheme>_default), which
resolves the cylinder twice as fine as the coarse meshes.
The errors with respect to the FeatFlow reference values (reference/featflow.json) and the
experimental orders of convergence are written to results/convergence_2d1.md.

Cases 2D-2 and 2D-3 (--case 2, --case 3): one time-dependent run (Re = 100) with the PQ1Bubble/Box
or the hybrid PQ1Bubble/Box scheme (--scheme) and the DIRK2 method (--time-scheme). These runs take
hours, so they are started in the background, detached from the terminal (--foreground waits for
them). The process id is written to <name>.pid in the run folder, and the run is stopped with
kill $(cat <run folder>/<name>.pid).

All cases use central differences for the advective term (Flux.UpwindWeight = 0.5). Each run gets
its own folder runs/2d<case>/<name>/ with the log, the benchmark indicators (<name>.json) and the
VTK output (vtk/). Existing results are never overwritten unless --overwrite is given, and a failed
run of the convergence study is reported and skipped. visualize.py evaluates and plots the results.

Usage (in the build directory test/freeflow/navierstokes/flowaroundcylinder):
    python3 benchmark.py
    python3 benchmark.py --schemes pq1bubble_box --levels 0 1 2
    python3 benchmark.py --case 2
    python3 benchmark.py --case 3 --dt 0.00125
    python3 benchmark.py --case 2 --scheme hybrid_pq1bubble_box
    python3 benchmark.py -- -Newton.MaxRelativeShift 1e-10   # further simulation parameters
    python3 visualize.py --case 1
"""

import argparse
import json
import math
import os
import resource
import subprocess
import sys

SCRIPT_DIR = os.path.dirname(os.path.realpath(__file__))
RESULTS_DIR = os.path.join(SCRIPT_DIR, "results")

# scheme: (executable, coarse mesh of the convergence study)
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
    "hybrid_pq1bubble_box": (
        "test_ff_navierstokes_dfg_benchmark_stationary_hybrid_pq1bubble_box",
        "meshes/cylinder_channel_coarse.msh",
    ),
}

# default meshes of the regression tests (mesh size 0.01 in the channel, 0.002 on the cylinder)
DEFAULT_MESHES = {
    "diamond_tpfa": "meshes/cylinder_channel.msh",
    "pq1bubble_box": "meshes/cylinder_channel.msh",
    "pq1bubble_diamond": "meshes/cylinder_channel_quad.msh",
    "hybrid_pq1bubble_box": "meshes/cylinder_channel.msh",
}

# the time-dependent cases are only implemented for the schemes with the new problem interface
TRANSIENT_SCHEMES = ["pq1bubble_box", "hybrid_pq1bubble_box"]


def readReferences():
    """FeatFlow reference values of case 2D-1 (reference/featflow.json)"""
    fileName = os.path.join(SCRIPT_DIR, "reference", "featflow.json")
    with open(fileName, encoding="utf-8") as jsonFile:
        featflow = json.load(jsonFile)["1"]
    return {"cDrag": featflow["CD"], "cLift": featflow["CL"], "pDiff": featflow["dP"]}


def runsDir(case):
    """Folder with the runs of a case"""
    return os.path.join("runs", f"2d{case}")


def runName(scheme, level=None, dt=None, timeScheme=None):
    """Name of a run (folder, problem name and JSON file), level None for the default mesh"""
    name = f"{scheme}_default" if level is None else f"{scheme}_r{level}"
    if dt is not None:
        name += f"_dt{dt:g}"
    return name if timeScheme is None else f"{name}_{timeScheme.lower()}"


def resultFile(case, name):
    """Path of the JSON file with the benchmark indicators of a run"""
    return os.path.join(runsDir(case), name, f"{name}.json")


def limitMemory(gigabytes):
    """Return a function limiting the address space of the child process"""

    def setLimit():
        limit = int(gigabytes * 1024**3)
        hardLimit = resource.getrlimit(resource.RLIMIT_AS)[1]
        if hardLimit != resource.RLIM_INFINITY:
            limit = min(limit, hardLimit)
        resource.setrlimit(resource.RLIMIT_AS, (limit, hardLimit))

    return setLimit if gigabytes > 0 else None


def simulationCommand(executable, params, name, extraArgs):
    """Command line of a run, executed in the vtk/ folder of the run"""
    return [
        os.path.abspath(executable),
        os.path.abspath(params),
        "-Flux.UpwindWeight",
        "0.5",
        "-Problem.EnableInertiaTerms",
        "true",
        "-Problem.Name",
        name,
        "-Problem.IndicatorFile",
        os.path.join("..", name),
    ] + extraArgs


###############################################################################
# Case 2D-1: grid convergence study
###############################################################################


def runStationary(scheme, level, args, extraArgs):
    """
    Run case 2D-1 for one scheme on the given level of the coarse mesh (level None: default
    mesh of the regression tests), return True on success
    """
    executable, mesh = SCHEMES[scheme]
    if level is None:
        mesh = DEFAULT_MESHES[scheme]
    name = runName(scheme, level)
    folder = os.path.join(runsDir(1), name)
    if os.path.exists(resultFile(1, name)) and not args.overwrite:
        print(f"Skipping {name} (results exist, use --overwrite to rerun)")
        return True

    vtkFolder = os.path.join(folder, "vtk")
    os.makedirs(vtkFolder, exist_ok=True)
    gridArgs = ["-Grid.File", os.path.abspath(mesh), "-Grid.Refinement", str(level or 0)]
    cmd = simulationCommand(executable, "params_1.input", name, gridArgs + extraArgs)

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
    """The benchmark indicators of a finished run of case 2D-1, None if it does not exist"""
    fileName = resultFile(1, name)
    if not os.path.exists(fileName):
        return None
    with open(fileName, encoding="utf-8") as jsonFile:
        return json.load(jsonFile)


def readStationaryResults(schemes, levels):
    """
    The benchmark indicators of all finished runs of case 2D-1 as
    {scheme: ([result per level], default)}, where default is the result on the default mesh
    (or None)
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


def tableRow(run, previous, label, references):
    """Table row of a run with the errors and the orders of convergence (w.r.t. the previous run)"""
    row = [label, str(run["numElements"]), str(run["numDofs"]["total"])]
    for key, reference in references.items():
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
    references = readReferences()
    referenceText = ", ".join(f"{key} = {value}" for key, value in references.items())
    lines = [
        "# Grid convergence of case 2D-1 (Re = 20, Flux.UpwindWeight = 0.5)",
        "",
        f"Errors with respect to the FeatFlow reference values {referenceText}.",
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
            tableRow(run, runs[i - 1] if i > 0 else None, str(run["refinement"]), references)
            for i, run in enumerate(runs)
        ]
        if default is not None:
            lines.append(tableRow(default, None, "default mesh", references))
    os.makedirs(os.path.dirname(fileName), exist_ok=True)
    with open(fileName, "w", encoding="utf-8") as tableFile:
        tableFile.write("\n".join(lines) + "\n")
    print(f"Wrote {fileName}")


def runConvergenceStudy(args, extraArgs):
    """Run case 2D-1 for all schemes and levels and write the convergence table"""
    for scheme in args.schemes:
        executable = SCHEMES[scheme][0]
        if not os.path.exists(executable):
            print(f"Executable {executable} not found, build it first. Skipping {scheme}.")
            continue
        for level in args.levels + [None]:
            runStationary(scheme, level, args, extraArgs)

    writeTable(
        readStationaryResults(args.schemes, args.levels),
        os.path.join(args.results_dir, "convergence_2d1.md"),
    )


###############################################################################
# Cases 2D-2 and 2D-3: time-dependent runs
###############################################################################


def runningProcess(pidFile):
    """Process id stored in the pid file if the process is still running, else None"""
    if not os.path.exists(pidFile):
        return None
    with open(pidFile, encoding="utf-8") as file:
        pid = int(file.read().strip())
    try:
        os.kill(pid, 0)
    except OSError:
        return None
    return pid


def transientArgs(case, args):
    """Mesh, grid level, time step size and end time of a time-dependent run"""
    mesh = os.path.abspath(os.path.join("meshes", "cylinder_channel.msh"))
    result = ["-Grid.File", mesh, "-Grid.Refinement", str(args.level)]
    result += ["-TimeLoop.Scheme", args.time_scheme]
    if args.dt is not None:
        result += ["-TimeLoop.Dt2" if case == 2 else "-TimeLoop.Dt", str(args.dt)]
    if args.tend is not None:
        result += ["-TimeLoop.T2" if case == 2 else "-TimeLoop.TEnd", str(args.tend)]
    return result


def runTransient(case, args, extraArgs):
    """Run case 2D-2 or 2D-3 (in the background unless --foreground is given)"""
    executable = SCHEMES[args.scheme][0]
    if not os.path.exists(executable):
        sys.exit(f"Executable {executable} not found, build it first.")

    name = runName(args.scheme, args.level, args.dt, args.time_scheme)
    folder = os.path.join(runsDir(case), name)
    pidFile = os.path.join(folder, f"{name}.pid")
    pid = runningProcess(pidFile)
    if pid is not None:
        sys.exit(f"{name} is still running (pid {pid}). Stop it first (kill {pid}).")
    if os.path.exists(os.path.join(folder, f"{name}.log")) and not args.overwrite:
        sys.exit(f"{folder}/ already contains results. Move them away or pass --overwrite.")

    vtkFolder = os.path.join(folder, "vtk")
    os.makedirs(vtkFolder, exist_ok=True)
    cmd = simulationCommand(executable, f"params_{case}.input", name, transientArgs(case, args))
    cmd += extraArgs
    logName = os.path.join(folder, f"{name}.log")
    print(f"Running {name}: " + " ".join(cmd))
    with open(logName, "w", encoding="utf-8") as logFile:
        options = {
            "cwd": vtkFolder,
            "stdin": subprocess.DEVNULL,
            "stdout": logFile,
            "stderr": subprocess.STDOUT,
            "preexec_fn": limitMemory(args.memory_limit),
        }
        if args.foreground:
            subprocess.run(cmd, check=False, **options)
            return

        # a new session detaches the process from the terminal, so it survives closing it
        process = subprocess.Popen(  # pylint: disable=consider-using-with
            cmd, start_new_session=True, **options
        )
    with open(pidFile, "w", encoding="utf-8") as file:
        file.write(f"{process.pid}\n")
    print(
        f"Started {name} in the background (pid {process.pid}).\n"
        f"  log:  tail -f {logName}\n"
        f"  stop: kill {process.pid}"
    )


###############################################################################
# Main
###############################################################################


def main():
    """Run the benchmark cases"""
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("--case", type=int, choices=[1, 2, 3], default=1, help="benchmark case")
    parser.add_argument(
        "--schemes",
        nargs="+",
        choices=list(SCHEMES),
        default=list(SCHEMES),
        help="2D-1: schemes",
    )
    parser.add_argument(
        "--levels", nargs="+", type=int, default=[0, 1, 2, 3], help="2D-1: Grid.Refinement values"
    )
    parser.add_argument(
        "--scheme",
        choices=TRANSIENT_SCHEMES,
        default=TRANSIENT_SCHEMES[0],
        help="2D-2, 2D-3: scheme",
    )
    parser.add_argument(
        "--level", type=int, default=0, help="2D-2, 2D-3: Grid.Refinement of the default mesh"
    )
    parser.add_argument(
        "--dt", type=float, default=None, help="2D-2 (fine step), 2D-3: fixed time step size"
    )
    parser.add_argument("--tend", type=float, default=None, help="2D-2, 2D-3: end time")
    # the choices of TimeLoop.Scheme in main.cc
    parser.add_argument(
        "--time-scheme",
        choices=["DIRK2", "DIRK3", "ImplicitEuler", "CrankNicolson"],
        default="DIRK2",
        help="2D-2, 2D-3: time stepping scheme (TimeLoop.Scheme)",
    )
    parser.add_argument(
        "--foreground", action="store_true", help="2D-2, 2D-3: wait for the run to finish"
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

    if args.case == 1:
        runConvergenceStudy(args, extraArgs)
    else:
        runTransient(args.case, args, extraArgs)


if __name__ == "__main__":
    main()
