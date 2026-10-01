#!/usr/bin/env python3
# SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later
"""Build, run, and visualize the 1pni heat convection benchmark.

The numerical solution is computed on a sequence of refined grids and compared
against two analytical references, each in the measure it is valid for:

- the retarded step front that the test problem itself writes out as
  ``temperatureExact`` neglects heat conduction, so it carries no information
  about the shape of the front, but its position x = v_T t is exact. It is
  therefore compared against the barycenter of the numerical profile only.
- the Ogata-Banks solution of the advection-diffusion equation is the solution
  of the equation the model actually solves, so the whole profile is compared
  against it in the L2 norm.
"""

import glob
import math
import subprocess
import xml.etree.ElementTree as ET
from dataclasses import dataclass
from pathlib import Path


TARGET = "test_1pni_convection_tpfa"
INPUT_FILE = "params_convection.input"
BASE_NAME = "1pni_1d_convection_benchmark"
TEMPERATURE_FIELD = "T"
STEP_FIELD = "temperatureExact"

# the thermal front reaches x ~ 4.25 m at the end of the simulation, so a
# shortened domain still keeps the outlet far behind the front while making the
# refined runs affordable
DOMAIN_LENGTH = 8.0
# cell counts of the refinement sequence; the time step is refined along with
# the grid so that the temporal numerical diffusion shrinks at the same rate
REFINEMENTS = (100, 400, 1600, 6400)
TIME_STEP_SIZES = (160.0, 40.0, 10.0, 2.5)

X_LIMITS = (2.5, 6.5)
TEMPERATURE_LIMITS = (289.9, 291.1)

# ordinal blue ramp (light -> dark) encoding the ordered refinement levels
PROFILE_COLORS = ("#a6c9f0", "#5b9ae4", "#2a6fc4", "#0f3f7a")


@dataclass(frozen=True)
class Case:
    cells: int
    max_dt: float
    name: str


@dataclass(frozen=True)
class Result:
    case: Case
    time: float
    vtp: Path
    front_velocity: float
    heat_capacity: float


def root_dir() -> Path:
    return Path(__file__).resolve().parents[4]


def case_build_dir() -> Path:
    return root_dir() / "build-cmake/test/porousmediumflow/1p/nonisothermal"


def run(command: list[str], cwd: Path | None = None) -> str:
    print(f"+ {' '.join(command)}")
    completed = subprocess.run(command, cwd=cwd, check=True,
                               capture_output=True, text=True)
    return completed.stdout


def last_reported(output: str, label: str) -> float:
    """Read the last value that the test problem reported under the given label.

    `OnePNIConvectionProblem::updateExactTemperature()` prints the material
    properties it derives from the current solution once per time step. Taking
    them from there keeps the analytical references consistent with the
    temperature- and pressure-dependent water properties of the simulation,
    which is more accurate than reading the front position off the discrete
    `temperatureExact` field, whose jump is quantized by the cell width.
    """
    values = [line.split(label, 1)[1] for line in output.splitlines() if label in line]
    if not values:
        raise RuntimeError(f"The test problem did not report '{label}'")

    return float(values[-1])


def remove_old_outputs(name: str) -> None:
    for pattern in (f"{name}-*.vtp", f"{name}.pvd", f"{name}-*.pvtp"):
        for file_name in glob.glob(str(case_build_dir() / pattern)):
            Path(file_name).unlink()


def final_output(case: Case, output: str) -> Result:
    """Read the last output time and the corresponding VTP file from the collection file."""
    collection = case_build_dir() / f"{case.name}.pvd"
    if not collection.exists():
        raise FileNotFoundError(f"No VTP output found for {case.name}")

    outputs = sorted(
        (float(data_set.attrib["timestep"]), data_set.attrib["file"])
        for data_set in ET.parse(collection).getroot().iter("DataSet")
    )
    time, file_name = outputs[-1]

    return Result(case=case, time=time, vtp=case_build_dir() / file_name,
                  front_velocity=last_reported(output, "retarded velocity:"),
                  heat_capacity=last_reported(output, "storage:"))


def build_and_run() -> list[Result]:
    build_dir = root_dir() / "build-cmake"
    run(["make", TARGET], cwd=build_dir)

    results = []
    for cells, max_dt in zip(REFINEMENTS, TIME_STEP_SIZES):
        case = Case(cells=cells, max_dt=max_dt, name=f"{BASE_NAME}_{cells}")
        remove_old_outputs(case.name)
        output = run([
            str(case_build_dir() / TARGET),
            INPUT_FILE,
            "-Problem.Name", case.name,
            "-Grid.UpperRight", str(DOMAIN_LENGTH),
            "-Grid.Cells", str(cells),
            "-TimeLoop.MaxTimeStepSize", str(max_dt),
            # only the final time step is needed for the comparison
            "-Problem.OutputInterval", "1000000",
        ], cwd=case_build_dir())
        results.append(final_output(case, output))

    return results


def require_plot_modules():
    try:
        import matplotlib.pyplot as plt
        import pyvista as pv
    except ImportError as error:
        raise SystemExit(
            "Post-processing requires PyVista and Matplotlib. "
            "Install them with: python3 -m pip install pyvista matplotlib"
        ) from error

    return pv, plt


def sorted_profile(mesh, field: str):
    """Return the x-coordinates and values of a cell or vertex field, sorted along x."""
    if field in mesh.cell_data:
        positions, values = mesh.cell_centers().points[:, 0], mesh.cell_data[field]
    elif field in mesh.point_data:
        positions, values = mesh.points[:, 0], mesh.point_data[field]
    else:
        raise KeyError(f"Missing field '{field}'")

    profile: dict[float, list[float]] = {}
    for position, value in zip(positions, values):
        profile.setdefault(round(float(position), 9), []).append(float(value))

    x = sorted(profile)
    return x, [sum(profile[position]) / len(profile[position]) for position in x]


def erfcx(z: float) -> float:
    """erfcx(z) = exp(z^2) erfc(z), for z >= 0.

    Used instead of erfc because at the Peclet number of this benchmark erfc
    underflows to zero while the exp(z^2) it is multiplied by still overflows.
    """
    if z < 25.0:
        return math.exp(z * z) * math.erfc(z)

    # asymptotic expansion, relative error below 1e-9 for z >= 25
    inverse = 1.0 / (2.0 * z * z)
    series = 1.0 + inverse * (-1.0 + inverse * (3.0 + inverse * (-15.0 + inverse * 105.0)))
    return series / (z * math.sqrt(math.pi))


class Analytical:
    """Analytical references for the retarded thermal front.

    The energy balance of the benchmark reduces to the advection-diffusion
    equation dT/dt + v dT/dx = D d2T/dx2, with the retarded front velocity
    v = v_D rho_w c_w / C_tot and the effective thermal diffusivity
    D = lambda_eff / C_tot.
    """

    def __init__(self, front_velocity: float, diffusivity: float,
                 t_low: float, t_high: float):
        self.v = front_velocity
        self.d = diffusivity
        self.t_low = t_low
        self.t_high = t_high

    def _scale(self, theta: float) -> float:
        return self.t_low + (self.t_high - self.t_low) * theta

    def ogata_banks(self, x: float, time: float) -> float:
        """Solution for a semi-infinite domain with a constant inlet temperature."""
        spread = 2.0 * math.sqrt(self.d * time)
        upstream = (x - self.v * time) / spread
        downstream = (x + self.v * time) / spread
        # exp(v x / D) erfc(downstream) == exp(-upstream^2) erfcx(downstream),
        # which avoids the overflow of exp(v x / D) at large Peclet numbers
        reflected = math.exp(-upstream**2) * erfcx(downstream)

        return self._scale(0.5 * (math.erfc(upstream) + reflected))

    def sharp_front(self, time: float) -> float:
        """Position of the pure-advection step front, x = v t.

        This is the field that the test problem writes out as
        ``temperatureExact``. It is compared against the barycenter of the
        numerical profile only, see `barycenter()`.
        """
        return self.v * time

    def step_profile(self, time: float, x_from: float, x_to: float):
        """The pure-advection step front as a profile, for plotting.

        Returned as the four corner points of the step rather than as samples of
        a discontinuous function, so that the riser stays vertical whatever the
        sampling of the axis.
        """
        front = self.sharp_front(time)
        return ([x_from, front, front, x_to],
                [self.t_high, self.t_high, self.t_low, self.t_low])

    def exact_barycenter(self, time: float) -> float:
        """Barycenter of the Ogata-Banks profile, which is exactly v t + D / v.

        The constant-temperature inlet of the type-I boundary condition lets
        heat in by conduction as well as by advection. That extra influx is
        largest while the front is still near the inlet and dies away
        afterwards, so it shifts the barycenter ahead of the step front by an
        amount that no longer depends on time.
        """
        return self.v * time + self.d / self.v


def front_parameters(results: list[Result]) -> Analytical:
    """Build the analytical references from the values the simulation reported.

    The front velocity v = v_D rho_w c_w / C_tot and the total volumetric heat
    capacity C_tot are taken from the test problem itself, so that the
    references use the same temperature- and pressure-dependent water
    properties as the simulation. Only the thermal conductivities, which the
    problem does not report, are taken from the benchmark parameters.
    """
    pv, _ = require_plot_modules()
    _, step = sorted_profile(pv.read(results[-1].vtp), STEP_FIELD)

    porosity = 0.4
    # thermal conductivity of water at the initial state, and of the solid as
    # set in params_convection.input; the model averages them by volume fraction
    lambda_fluid, lambda_solid = 0.598, 2.8
    lambda_eff = porosity * lambda_fluid + (1 - porosity) * lambda_solid

    return Analytical(front_velocity=results[-1].front_velocity,
                      diffusivity=lambda_eff / results[-1].heat_capacity,
                      t_low=min(step), t_high=max(step))


def barycenter(x: list[float], temperature: list[float], cell_width: float,
               analytical: Analytical) -> float:
    """Barycenter of the thermal front of a temperature profile.

    Defined as the first moment of the temperature gradient,

        x_bary = int x (-dT/dx) dx / int (-dT/dx) dx ,

    which, integrated by parts over a profile that runs from T_high at the
    inlet to T_low ahead of the front, is the equivalent step position

        x_bary = int (T - T_low) dx / (T_high - T_low) .

    It is the position of the sharp front that stores the same amount of heat
    as the smeared profile, so it is insensitive to how wide the front is and
    isolates how far the front has travelled.
    """
    excess = sum((value - analytical.t_low) * cell_width for value in temperature)
    return excess / (analytical.t_high - analytical.t_low)


def l2_error(x: list[float], numerical: list[float], reference, time: float,
             cell_width: float) -> float:
    """Discrete L2 norm of the difference, normalized by the domain length."""
    squared = sum((reference(position, time) - value)**2 * cell_width
                  for position, value in zip(x, numerical))
    return math.sqrt(squared / DOMAIN_LENGTH)


def convergence_table(results: list[Result], analytical: Analytical) -> None:
    """Report both comparisons: the front position and the front shape.

    The step front is compared against the barycenter only, because it is the
    correct reference for how far the front has travelled but not for how wide
    it is. The shape is compared against the Ogata-Banks solution in the L2 norm.
    """
    pv, _ = require_plot_modules()

    time = results[-1].time
    sharp_front = analytical.sharp_front(time)
    print("\nConvergence at t = %.0f s" % time)
    print("  sharp front position  v_T t   = %.6f m" % sharp_front)
    print("  Ogata-Banks barycenter        = %.6f m"
          % analytical.exact_barycenter(time))
    print()
    print("  cells  dx [m]  Pe_grid |  barycenter [m]  minus v_T t [m] |  L2 vs Ogata-Banks [K]   rate")

    previous = None
    for result in results:
        x, temperature = sorted_profile(pv.read(result.vtp), TEMPERATURE_FIELD)
        cell_width = DOMAIN_LENGTH / result.case.cells

        centre = barycenter(x, temperature, cell_width, analytical)
        shape_error = l2_error(x, temperature, analytical.ogata_banks, time, cell_width)

        rate = "     -"
        if previous is not None:
            rate = "  %5.2f" % math.log2(previous / shape_error)

        print("  %5d  %6.4f  %7.2f |  %12.6f  %+15.2e |  %.6f            %s"
              % (result.case.cells, cell_width, analytical.v * cell_width / analytical.d,
                 centre, centre - sharp_front, shape_error, rate))
        previous = shape_error


def create_line_plot(results: list[Result], analytical: Analytical, image_file: Path) -> None:
    pv, plt = require_plot_modules()

    fig, ax = plt.subplots(figsize=(8.2, 4.8), constrained_layout=True)
    time = results[-1].time

    # only profiles are drawn; the barycenters are reported as numbers by
    # convergence_table(), since their distance to v_T t is a few millimetres
    # and therefore invisible on the scale of the front
    for index, result in enumerate(results):
        x, temperature = sorted_profile(pv.read(result.vtp), TEMPERATURE_FIELD)
        ax.plot(x, temperature, linewidth=2, color=PROFILE_COLORS[index % len(PROFILE_COLORS)],
                label=rf"$T$ numerical ({result.case.cells} cells)")

    samples = [X_LIMITS[0] + (X_LIMITS[1] - X_LIMITS[0]) * index / 800 for index in range(801)]
    ax.plot(samples, [analytical.ogata_banks(position, time) for position in samples],
            linewidth=1.6, color="black", linestyle="--", label=r"$T$ exact (Ogata-Banks)")

    step_x, step_t = analytical.step_profile(time, *X_LIMITS)
    ax.plot(step_x, step_t, linewidth=1.4, color="#b03030", linestyle=":",
            label=r"$T$ retarded step front")

    ax.set_xlabel("x [m]")
    ax.set_ylabel(r"Temperature $T$ [K]")
    ax.set_title(rf"Thermal front at $t$ = {time / 3600:.2f} h")
    ax.set_xlim(*X_LIMITS)
    ax.set_ylim(*TEMPERATURE_LIMITS)
    ax.grid(True, alpha=0.3)
    ax.legend(loc="upper right", fontsize=9)

    fig.savefig(image_file, dpi=200)
    plt.close(fig)


def main() -> None:
    results = build_and_run()
    analytical = front_parameters(results)

    print(f"\nretarded front velocity: {analytical.v:.6e} m/s")
    print(f"thermal diffusivity:     {analytical.d:.6e} m^2/s")

    convergence_table(results, analytical)

    print("\nCreating line plot...")
    create_line_plot(results, analytical, case_build_dir() / f"{BASE_NAME}_lineplot.png")


if __name__ == "__main__":
    main()
