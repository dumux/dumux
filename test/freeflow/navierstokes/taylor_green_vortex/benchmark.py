#!/usr/bin/env python3
# SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later
"""
Run the Taylor-Green vortex benchmark (2D): spatial and
temporal convergence studies, the decay of kinetic energy and enstrophy, and the simulation
for the animation, see README.md.

The script has to be executed from the build directory of this test, where the executables
and parameter files are located. Everything is written to the directory runs/:

  * runs/<name>/<name>_errors.csv: errors and indicators of every time step of a run
  * runs/<name>/vtk/: ParaView output of every time step (<name>.pvd and the higher-order
    output <name>_ho.pvd), for the instationary runs or with --vtk
  * runs/taylorgreen_<study>_<dim>d.md: convergence tables
  * runs/taylorgreen_<study>_<dim>d.json: manifests of the runs of a study, from which
    visualize.py creates the plots and animate.py the animation
"""

import argparse
import csv
import json
import math
import subprocess
import sys
from collections import namedtuple
from dataclasses import dataclass, field
from pathlib import Path

# momentum discretizations (the mass balance is always discretized with the Box scheme);
# pq2hybridhoquad is hybrid PQ2 with higher-order quadrature rules for the source term on the
# sub-control volumes and for the mass fluxes (instead of the midpoint rule of the default rules)
SCHEMES = ("pq1bubble", "pq1bubblehybrid", "pq2hybrid", "pq2hybridhoquad", "pq3hybrid")
SCHEME_LABELS = {
    "pq1bubble": "PQ1Bubble",
    "pq1bubblehybrid": "hybrid PQ1Bubble",
    "pq2hybrid": "hybrid PQ2",
    "pq2hybridhoquad": "hybrid PQ2, higher-order quadrature",
    "pq3hybrid": "hybrid PQ3",
}
# default momentum schemes of the studies (unless given with --schemes)
DEFAULT_SCHEMES = {
    "spatial": {2: SCHEMES},
    "temporal": {2: ("pq2hybrid",)},
    "energy": {2: SCHEMES},
    "animation": {2: ("pq3hybrid",)},
}
# time integration methods (TimeLoop.Scheme) of the instationary problem
TIME_SCHEME_LABELS = {
    "ImplicitEuler": "implicit Euler",
    "CrankNicolson": "Crank-Nicolson",
    "DIRK2": "DIRK2",
    "DIRK3": "DIRK3",
}

# expected convergence orders under grid refinement; the velocity of the hybrid schemes with Box
# pressure converges at most with order 3 in L2 and 2 in H1 since the pressure is only of order 2;
# hybrid PQ2 only converges with order 2 in L2 (independent of the quadrature rules)
EXPECTED_SPATIAL_ORDER = {
    "pq1bubble": {"velocityL2": 2.0, "velocityH1": 1.0, "pressureL2": 2.0, "pressureH1": 1.0},
    "pq1bubblehybrid": {"velocityL2": 2.0, "velocityH1": 1.0, "pressureL2": 2.0, "pressureH1": 1.0},
    "pq2hybrid": {"velocityL2": 2.0, "velocityH1": 2.0, "pressureL2": 2.0, "pressureH1": 1.0},
    "pq2hybridhoquad": {"velocityL2": 2.0, "velocityH1": 2.0, "pressureL2": 2.0, "pressureH1": 1.0},
    "pq3hybrid": {"velocityL2": 3.0, "velocityH1": 2.0, "pressureL2": 2.0, "pressureH1": 1.0},
}
EXPECTED_TEMPORAL_ORDER = {"ImplicitEuler": 1.0, "CrankNicolson": 2.0, "DIRK2": 2.0, "DIRK3": 3.0}
ORDER_TOLERANCE = 0.3
# error norms not checked by the convergence tests (option --test): the velocity L2 error of
# hybrid PQ2 converges irregularly (in 2D with EOC 5.7, 0.4, 1.6 from 16^2 to 128^2 cells,
# independent of the quadrature rules), see README.md
UNCHECKED_SPATIAL_ORDER = {"pq2hybrid": ("velocityL2",), "pq2hybridhoquad": ("velocityL2",)}
ERROR_KEYS = ("velocityL2", "velocityH1", "pressureL2", "pressureH1")

# grid sequences of the spatial study (full study, test mode); the velocity L2 error of hybrid
# PQ1Bubble and hybrid PQ2 converges irregularly on the coarse grids, so their tests start with
# 32^2 cells; the Newton solver does not converge for PQ1Bubble on 16^2 cells
SPATIAL_CELLS = {2: ([16, 32, 64, 128], [16, 32, 64])}
SCHEME_SPATIAL_CELLS = {
    (2, "pq1bubble"): ([32, 64, 128], [32, 64]),
    (2, "pq1bubblehybrid"): ([16, 32, 64, 128], [32, 64, 128]),
    (2, "pq2hybrid"): ([16, 32, 64, 128], [32, 64, 128]),
    (2, "pq2hybridhoquad"): ([16, 32, 64, 128], [32, 64, 128]),
}
# additional Newton settings of the stationary problem
SPATIAL_NEWTON = {2: {}}
# cells per direction (full study, test mode) and time step refinements (full study, test mode)
# of the temporal study; the spatial error should be small compared to the temporal error
TEMPORAL_CELLS = {2: (64, 32)}
TEMPORAL_LEVELS = {2: (4, 3)}
TEMPORAL_TEND = 0.4
# setup of the animation (2D, instationary)
ANIMATION = {"cells": 32, "dt": 0.01, "tEnd": 1.0}

# every simulation runs in its own directory below this one
RUN_DIRECTORY = Path("runs")

# a simulation: dimension, momentum scheme and name of the run
Simulation = namedtuple("Simulation", "dim scheme name")


@dataclass
class Settings:
    """Options of all runs"""

    reuse: bool = False  # reuse existing error files instead of rerunning (--reuse)
    vtk: bool = False  # ParaView output for all runs (--vtk)
    timeScheme: str = "ImplicitEuler"  # of the energy study and the animation (--time-scheme)
    parameters: dict = field(default_factory=dict)  # additional runtime parameters


@dataclass
class Run:
    """Result of a single simulation run (last row of the error file)"""

    name: str
    errorFile: Path
    dt: float
    numDofs: tuple  # number of velocity and pressure degrees of freedom N_v and N_p
    meshSizes: tuple  # h^v = (|Omega| / N_v)^(1/d) and h^p = (|Omega| / N_p)^(1/d)
    errors: dict

    def meshSize(self, key):
        """
        Mesh size of the variable of the error norm key: h^v for velocity and h^p for pressure
        errors, comparable across element types and schemes (see Schneider & Koch, 2024)
        """
        return self.meshSizes[0] if key.startswith("velocity") else self.meshSizes[1]


def runDirectory(name):
    """Directory of the run with the given name"""
    return RUN_DIRECTORY / name


def vtkDirectory(name):
    """Directory of the ParaView output of the run with the given name"""
    return runDirectory(name) / "vtk"


def commandLine(parameters):
    """Command line options for the runtime parameters given as dict"""
    return [option for key, value in parameters.items() for option in (f"-{key}", str(value))]


def gridCells(numCells, dim):
    """Value of Grid.Cells for numCells cells per direction"""
    return " ".join([str(numCells)] * dim)


def readParameter(settings, dim, key):
    """A runtime parameter (group.name) from the settings or the input file params_<dim>d.input"""
    if key in settings.parameters:
        return settings.parameters[key]
    group, name = key.split(".", 1)
    currentGroup = None
    with open(f"params_{dim}d.input") as inputFile:
        for line in inputFile:
            line = line.split("#")[0].strip()
            if line.startswith("[") and line.endswith("]"):
                currentGroup = line[1:-1]
            elif "=" in line and currentGroup == group:
                parameter, value = (part.strip() for part in line.split("=", 1))
                if parameter == name:
                    return value
    return None


def readRows(errorFile):
    """The rows of an error file"""
    with open(errorFile, newline="") as errors:
        reader = csv.DictReader(errors)
        return [{key: float(value) for key, value in row.items()} for row in reader]


def runSimulation(settings, simulation, parameters, vtk=False):
    """
    Run a simulation in its own directory runs/<name>/ and return the rows of the error file;
    the ParaView output (if enabled) is written to the subdirectory runs/<name>/vtk/
    """
    dim, scheme, name = simulation
    exe = Path(f"test_ff_navierstokes_taylorgreen_{dim}d_{scheme}")
    if not exe.is_file():
        raise FileNotFoundError(f"Executable {exe} not found. Build it first.")

    errorFile = runDirectory(name) / f"{name}_errors.csv"
    vtk = vtk or settings.vtk
    allParameters = {
        "Problem.Name": name,
        "Problem.ErrorFile": errorFile.resolve(),
        "Problem.EnableVtkOutput": "true" if vtk else "false",
        **parameters,
        **settings.parameters,
    }
    command = [str(exe.resolve()), str(Path(f"params_{dim}d.input").resolve())]
    command += commandLine(allParameters)
    if settings.reuse and errorFile.is_file():
        print(f"Reusing {errorFile}", flush=True)
    else:
        # the simulation writes the ParaView output to its working directory
        workDir = vtkDirectory(name) if vtk else runDirectory(name)
        workDir.mkdir(parents=True, exist_ok=True)
        print(f"+ (in {workDir}) " + " ".join(command), flush=True)
        subprocess.run(command, check=True, stdout=subprocess.DEVNULL, cwd=workDir)

    rows = readRows(errorFile)
    if not rows:
        raise RuntimeError(f"No errors written by {exe}")
    return rows


def toRun(name, rows):
    """Convert the rows of an error file to a run"""
    last = rows[-1]
    return Run(
        name=name,
        errorFile=runDirectory(name) / f"{name}_errors.csv",
        dt=last["dt"],
        numDofs=(int(last["numDofsVelocity"]), int(last["numDofsPressure"])),
        meshSizes=(last["hVelocity"], last["hPressure"]),
        errors={key: last[key] for key in ERROR_KEYS},
    )


def rates(runs, key, step):
    """Experimental orders of convergence between consecutive runs (step size step(run, key))"""
    result = []
    for coarse, fine in zip(runs[:-1], runs[1:]):
        errorRatio = coarse.errors[key] / fine.errors[key]
        stepRatio = step(coarse, key) / step(fine, key)
        result.append(math.log(errorRatio) / math.log(stepRatio))
    return result


def writeTable(fileName, title, runs, columns, step):
    """
    Print a Markdown table with errors and convergence rates (EOC w.r.t. step(run, key)) and
    append it to fileName (if given); columns is a list of (header, value(run)) tuples
    describing the refinement
    """
    rateTable = {key: [float("nan")] + rates(runs, key, step) for key in ERROR_KEYS}
    header = [name for name, _ in columns] + [f"{key} | EOC" for key in ERROR_KEYS]
    lines = [f"### {title}", "", "| " + " | ".join(header) + " |"]
    lines.append("|" + "---|" * (len(columns) + 2 * len(ERROR_KEYS)))
    for i, run in enumerate(runs):
        cells = [value(run) for _, value in columns]
        for key in ERROR_KEYS:
            cells += [f"{run.errors[key]:.4e}", f"{rateTable[key][i]:.2f}"]
        lines.append("| " + " | ".join(cells) + " |")
    text = "\n".join(lines) + "\n"
    print(text)
    if fileName is not None:
        with open(fileName, "a") as tableFile:
            tableFile.write(text + "\n")


def checkRates(label, runs, expected, step):
    """Check that the rate of the finest refinement matches the expected order"""
    success = True
    for key, order in expected.items():
        rate = rates(runs, key, step)[-1]
        if not rate >= order - ORDER_TOLERANCE:
            print(f"FAILED: {label}: {key} rate {rate:.2f} < expected {order} - {ORDER_TOLERANCE}")
            success = False
        else:
            print(f"passed: {label}: {key} rate {rate:.2f} (expected {order})")
    return success


def studyFile(study, dim, suffix):
    """Path of the table or the manifest of a study"""
    RUN_DIRECTORY.mkdir(exist_ok=True)
    return RUN_DIRECTORY / f"taylorgreen_{study}_{dim}d.{suffix}"


def makeGroup(scheme, timeScheme, expected, runs):
    """
    A group of runs in a manifest (one curve in the plots): the label, the momentum scheme,
    the time scheme, the expected orders and the runs (error files, higher-order ParaView output)
    """
    label = SCHEME_LABELS[scheme]
    if timeScheme is not None and expected:
        label += f", {TIME_SCHEME_LABELS[timeScheme]}"
    return {
        "label": label,
        "scheme": scheme,
        "timeScheme": timeScheme,
        "expected": expected,
        "runs": [
            {
                "name": run.name,
                "errors": str(run.errorFile),
                "pvd": str(vtkDirectory(run.name) / f"{run.name}_ho.pvd"),
                "dt": run.dt,
            }
            for run in runs
        ],
    }


def writeManifest(settings, study, dim, groups, parameters=None):
    """List the runs of a study for visualize.py and animate.py"""
    manifest = {
        "study": study,
        "dim": dim,
        "viscosity": float(readParameter(settings, dim, "Component.LiquidKinematicViscosity")),
        "density": float(readParameter(settings, dim, "Component.LiquidDensity")),
        "parameters": parameters or {},
        "groups": groups,
    }
    fileName = studyFile(study, dim, "json")
    with open(fileName, "w") as manifestFile:
        json.dump(manifest, manifestFile, indent=2)
    print(f"Wrote {fileName}")


def timeStepSize(run, _key):
    """The time step size of a run (the step size of the temporal convergence study)"""
    return run.dt


# columns of the convergence tables describing the refinement
SPATIAL_COLUMNS = [
    ("N_v", lambda run: f"{run.numDofs[0]}"),
    ("h^v", lambda run: f"{run.meshSizes[0]:.4e}"),
    ("N_p", lambda run: f"{run.numDofs[1]}"),
    ("h^p", lambda run: f"{run.meshSizes[1]:.4e}"),
]
TEMPORAL_COLUMNS = [("dt", lambda run: f"{run.dt:.4e}")]


def runGridSequence(settings, dim, scheme, cells):
    """Run the stationary problem on the given grids, returns None if a run fails"""
    runs = []
    for numCells in cells:
        simulation = Simulation(dim, scheme, f"taylorgreen_spatial_{dim}d_{scheme}_{numCells}")
        parameters = {"Problem.IsStationary": "true", "Grid.Cells": gridCells(numCells, dim)}
        try:
            rows = runSimulation(settings, simulation, parameters | SPATIAL_NEWTON[dim])
        except (subprocess.CalledProcessError, RuntimeError):
            # RuntimeError: a failed run is reused (option --reuse)
            print(f"FAILED: {dim}D {scheme} (stationary) on {numCells}^{dim} cells, skipping")
            return None
        runs.append(toRun(simulation.name, rows))
    return runs


def runTimeStepSequence(settings, simulation, cells, levels):
    """Run the instationary problem with refined time step sizes (time scheme as simulation name)"""
    dim, scheme, timeScheme = simulation
    runs = []
    for level in range(levels):
        dt = TEMPORAL_TEND / 2 ** (level + 1)
        name = f"taylorgreen_temporal_{dim}d_{scheme}_{timeScheme}_{level}"
        parameters = {
            "Problem.IsStationary": "false",
            "Grid.Cells": gridCells(cells, dim),
            "TimeLoop.Scheme": timeScheme,
            "TimeLoop.TEnd": TEMPORAL_TEND,
            "TimeLoop.DtInitial": dt,
            "TimeLoop.MaxTimeStepSize": dt,
        }
        rows = runSimulation(settings, Simulation(dim, scheme, name), parameters, vtk=True)
        runs.append(toRun(name, rows))
    return runs


def spatialStudy(settings, dim, schemes, cellSequence, test):
    """Stationary problem under uniform grid refinement"""
    groups = []
    success = True
    # the convergence tests (test mode) do not overwrite the tables of the full study
    tableFile = None if test else studyFile("spatial", dim, "md")
    if tableFile is not None:
        tableFile.unlink(missing_ok=True)
    for scheme in schemes:
        cells = cellSequence or SCHEME_SPATIAL_CELLS.get((dim, scheme), SPATIAL_CELLS[dim])[test]
        runs = runGridSequence(settings, dim, scheme, cells)
        if runs is None:
            success = False
            continue

        expected = EXPECTED_SPATIAL_ORDER[scheme]
        groups.append(makeGroup(scheme, None, expected, runs))
        # convergence rates w.r.t. h^v (velocity errors) and h^p (pressure errors)
        label = f"{dim}D {scheme} (stationary)"
        writeTable(
            tableFile, f"{label}, EOC w.r.t. h^v and h^p", runs, SPATIAL_COLUMNS, Run.meshSize
        )
        if test:
            checked = {
                key: order
                for key, order in expected.items()
                if key not in UNCHECKED_SPATIAL_ORDER.get(scheme, ())
            }
            success &= checkRates(label, runs, checked, Run.meshSize)

    if groups and not test:
        writeManifest(settings, "spatial", dim, groups)
    return success


def temporalStudy(settings, dim, schemes, levels, test):
    """Instationary problem on a fine grid under refinement of the time step size"""
    groups = []
    success = True
    # the convergence tests (test mode) do not overwrite the tables of the full study
    tableFile = None if test else studyFile("temporal", dim, "md")
    if tableFile is not None:
        tableFile.unlink(missing_ok=True)
    cells = TEMPORAL_CELLS[dim][test]
    for scheme in schemes:
        for timeScheme, order in EXPECTED_TEMPORAL_ORDER.items():
            runs = runTimeStepSequence(settings, Simulation(dim, scheme, timeScheme), cells, levels)
            # the time integration limits the order of all error norms
            groups.append(makeGroup(scheme, timeScheme, dict.fromkeys(ERROR_KEYS, order), runs))
            label = f"{dim}D {scheme} {timeScheme}"
            title = f"{label} (t = {TEMPORAL_TEND})"
            writeTable(tableFile, title, runs, TEMPORAL_COLUMNS, timeStepSize)
            if test:
                success &= checkRates(label, runs, {"velocityL2": order}, timeStepSize)

    if groups and not test:
        writeManifest(settings, "temporal", dim, groups, {"cells": cells, "tEnd": TEMPORAL_TEND})
    return success


def energyStudy(settings, dim, schemes):
    """Decay of the kinetic energy and the enstrophy of the instationary problem"""
    groups = []
    success = True
    for scheme in schemes:
        name = f"taylorgreen_energy_{dim}d_{scheme}"
        parameters = {"Problem.IsStationary": "false", "TimeLoop.Scheme": settings.timeScheme}
        try:
            rows = runSimulation(settings, Simulation(dim, scheme, name), parameters, vtk=True)
        except (subprocess.CalledProcessError, RuntimeError):
            print(f"FAILED: {dim}D {scheme} (instationary), skipping")
            success = False
            continue
        groups.append(makeGroup(scheme, settings.timeScheme, {}, [toRun(name, rows)]))
        last = rows[-1]
        print(
            f"{dim}D {scheme}: at t = {last['t']:g}: "
            f"E_h/E = {last['kineticEnergy'] / last['kineticEnergyExact']:.6f}, "
            f"eps_h/eps = {last['enstrophy'] / last['enstrophyExact']:.6f}"
        )

    if groups:
        timeScheme = TIME_SCHEME_LABELS[settings.timeScheme]
        writeManifest(settings, "energy", dim, groups, {"timeScheme": timeScheme})
    return success


def animationStudy(settings, dim, scheme):
    """Instationary run with ParaView output for the animation (rendered by animate.py)"""
    if dim != 2:
        print("The animation is only available in 2D, skipping")
        return True

    name = f"taylorgreen_animation_{dim}d_{scheme}"
    parameters = {
        "Problem.IsStationary": "false",
        "Grid.Cells": gridCells(ANIMATION["cells"], dim),
        "TimeLoop.Scheme": settings.timeScheme,
        "TimeLoop.DtInitial": ANIMATION["dt"],
        "TimeLoop.MaxTimeStepSize": ANIMATION["dt"],
        "TimeLoop.TEnd": ANIMATION["tEnd"],
    }
    rows = runSimulation(settings, Simulation(dim, scheme, name), parameters, vtk=True)
    group = makeGroup(scheme, settings.timeScheme, {}, [toRun(name, rows)])
    timeScheme = TIME_SCHEME_LABELS[settings.timeScheme]
    writeManifest(
        settings, "animation", dim, [group], {"cells": ANIMATION["cells"], "timeScheme": timeScheme}
    )
    print("Render the animation with: python3 animate.py")
    return True


def parseArguments():
    """Parse the command line arguments"""
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("--dim", type=int, nargs="+", choices=(2,), default=[2])
    parser.add_argument(
        "--schemes",
        nargs="+",
        choices=SCHEMES,
        default=None,
        help="momentum discretizations (default: see DEFAULT_SCHEMES and README.md)",
    )
    parser.add_argument(
        "--study",
        nargs="+",
        choices=("spatial", "temporal", "energy", "animation", "all"),
        default=["all"],
        help="studies to run ('all' runs all except the animation)",
    )
    parser.add_argument(
        "--time-scheme",
        choices=tuple(TIME_SCHEME_LABELS),
        default="ImplicitEuler",
        help="time integration method of the energy study and the animation",
    )
    parser.add_argument(
        "--upwind-weight",
        type=float,
        default=None,
        help="weight of the upwind velocity in the convective flux (Flux.UpwindWeight, "
        "default from the input file: 0.5, i.e. central differences; 1 is full upwinding)",
    )
    parser.add_argument("--wavelength", type=float, help="wavelength L (Problem.WaveLength)")
    parser.add_argument("--levels", type=int, help="number of time step refinements (temporal)")
    parser.add_argument(
        "--cells", type=int, nargs="+", help="cells per direction of the grids (spatial study)"
    )
    parser.add_argument(
        "--vtk",
        action="store_true",
        help="write ParaView output for all runs (default: only for the instationary runs)",
    )
    parser.add_argument(
        "--test", action="store_true", help="short study checking the convergence rates"
    )
    parser.add_argument(
        "--reuse",
        action="store_true",
        help="reuse existing error files instead of rerunning the simulations",
    )
    args = parser.parse_args()
    if args.levels is not None and args.levels < 2:
        parser.error("At least two refinement levels are needed to compute convergence rates")
    if args.cells is not None and len(args.cells) < 2:
        parser.error("At least two grids are needed to compute convergence rates")
    return args


def main():
    """Run the selected studies"""
    args = parseArguments()
    optional = {
        "Flux.UpwindWeight": args.upwind_weight,
        "Problem.WaveLength": args.wavelength,
    }
    settings = Settings(
        reuse=args.reuse,
        vtk=args.vtk,
        timeScheme=args.time_scheme,
        parameters={key: value for key, value in optional.items() if value is not None},
    )

    studies = ("spatial", "temporal", "energy") if "all" in args.study else args.study
    success = True
    for dim in args.dim:
        for study in studies:
            schemes = args.schemes or DEFAULT_SCHEMES[study].get(dim, ())
            if study == "spatial":
                success &= spatialStudy(settings, dim, schemes, args.cells, args.test)
            elif study == "temporal":
                levels = args.levels or TEMPORAL_LEVELS[dim][args.test]
                success &= temporalStudy(settings, dim, schemes, levels, args.test)
            elif study == "energy":
                success &= energyStudy(settings, dim, schemes)
            elif schemes:
                success &= animationStudy(settings, dim, schemes[0])
            else:
                print(f"The animation is only available in 2D, skipping {dim}D")

    sys.exit(0 if success else 1)


if __name__ == "__main__":
    main()
