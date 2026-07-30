#!/usr/bin/env python3
# SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later
"""Build, run, and visualize the heatpipe benchmark.

The semi-analytical reference solution implements the steady-state ODE system
derived by Udell and Fitch (1985) for the heat-pipe effect, in the form given
by Huang, Kolditz and Shao (2015, Geothermal Energy 3:14,
https://doi.org/10.1186/s40517-015-0030-8), whose reference MATLAB
implementation is shipped as part of the OpenGeoSys benchmark suite. The
governing equations are four coupled first-order ODEs in space for the
effective wetting-phase saturation, gas pressure, gas-phase air mole
fraction, and temperature, integrated from the left (Dirichlet) boundary
towards the heat source with a Runge-Kutta scheme, here via
``scipy.integrate.solve_ivp``.
"""

import glob
import subprocess
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp

TARGET = "test_heatpipe_box"
BASE_NAME = "heatpipe"

# Grid resolutions compared in the convergence study; the default test (and
# CTest) always runs at 120 cells (grids/heatpipe.dgf) for a fast runtime.
# See the README for why the residual mismatch near the dry-out front shrinks
# with resolution.
RESOLUTIONS = {
    120: "grids/heatpipe.dgf",
    240: "grids/heatpipe_240.dgf",
    480: "grids/heatpipe_480.dgf",
}

SATURATION_FIELD = "S_liq"
PRESSURE_FIELD = "p_gas"
TEMPERATURE_FIELD = "T"
AIR_MOLEFRACTION_FIELD = "x^Air_gas"


def root_dir() -> Path:
    return Path(__file__).resolve().parents[4]


def case_source_dir() -> Path:
    return Path(__file__).resolve().parent


def case_build_dir() -> Path:
    return root_dir() / "build-cmake/test/porousmediumflow/2p2c/heatpipe"


def run(command: list[str], cwd: Path | None = None) -> None:
    print(f"+ {' '.join(command)}")
    subprocess.run(command, cwd=cwd, check=True)


def remove_old_outputs(name: str) -> None:
    for pattern in (f"{name}-*.vtu", f"{name}.pvd", f"{name}-*.pvtu"):
        for file_name in glob.glob(str(case_build_dir() / pattern)):
            Path(file_name).unlink()


def latest_vtu(name: str) -> Path:
    files = sorted(case_build_dir().glob(f"{name}-*.vtu"))
    if not files:
        raise FileNotFoundError(f"No VTU output found for {name}")
    return files[-1]


def build_and_run_all() -> dict[int, Path]:
    """Build once, then run the grid convergence study at each resolution.

    Returns a dict mapping cell count -> path of its final VTU output.
    """
    build_dir = root_dir() / "build-cmake"
    run(["cmake", "--build", ".", "--target", TARGET], cwd=build_dir)

    vtus = {}
    for cells, grid_file in RESOLUTIONS.items():
        name = f"{BASE_NAME}_{cells}"
        remove_old_outputs(name)
        run(
            [
                str(case_build_dir() / TARGET),
                "params.input",
                "-Problem.Name",
                name,
                "-Grid.File",
                str(case_source_dir() / grid_file),
            ],
            cwd=case_build_dir(),
        )
        vtus[cells] = latest_vtu(name)
    return vtus


# ---------------------------------------------------------------------------
# Semi-analytical reference solution (Udell & Fitch 1985 / Huang et al. 2015)
# ---------------------------------------------------------------------------

# physical parameters, matching spatialparams.hh
# (dumux/material/fluidmatrixinteractions/2p/heatpipelaw.hh), problem.hh and params.input
K = 1.0e-12  # m^2,  Problem.Permeability in params.input
PHI = 0.4  # -,    HeatPipeSpatialParams::porosity_
SWR = 0.15  # -,    Swr passed to HeatPipeLaw::EffToAbsParams in spatialparams.hh
RHO_W = 958.4  # kg/m^3, water density near the phase-change reference state
MU_W = 2.938e-4  # Pa*s
MU_G_A = 2.08e-5  # Pa*s, dynamic viscosity of air
MU_G_W = 1.2e-5  # Pa*s, dynamic viscosity of steam
GAMMA = 0.0588  # N/m, surface tension of water; matches the literal value used in
# HeatPipeSpatialParams (spatialparams.hh), not the 0.05878 quoted in older literature
MW = 0.01801518  # kg/mol, molar mass of water; Dumux::Components::H2O::molarMass (IAPWS::Common)
MA = 0.02896  # kg/mol, molar mass of air; Dumux::Components::Air::molarMass()
R_GAS = 8.314472  # J/(mol*K); Dumux::Constants<Scalar>::R
H_WG = 2.258e6  # J/kg, latent heat of vaporization of water

# Effective thermal conductivity: DuMux's TwoPTwoCNI model uses
# ThermalConductivitySomertonTwoP by default (dumux/porousmediumflow/2p2c/model.hh),
# NOT a fixed pair of wet/dry constants. It computes
#   lambda_wet = lambda_solid^(1-phi) * lambda_liquid_fluid^phi
#   lambda_dry = lambda_solid^(1-phi) * lambda_gas_fluid^phi
# (geometric mean, see dumux/material/fluidmatrixinteractions/2p/thermalconductivity/somerton.hh)
# with lambda_solid = 2.8 W/(m*K) (params.input, Component.SolidThermalConductivity),
# lambda_gas_fluid = 0.0255535 W/(m*K) (Dumux::Components::Air::gasThermalConductivity,
# a constant), and lambda_liquid_fluid the real IAPWS water thermal conductivity
# (Dumux::Components::H2O::liquidThermalConductivity), here evaluated at a
# representative T=370 K, p=1.1e5 Pa (~0.676 W/(m*K)). These values differ
# substantially from the historical Emmert/Udell table values (1.13/0.582 W/(m*K)),
# which do not apply to this particular solid conductivity choice.
LAMBDA_WET = 1.586  # W/(m*K), Somerton wet endpoint (see above)
LAMBDA_DRY = 0.4278  # W/(m*K), Somerton dry endpoint (see above)

# Binary diffusion coefficient: the default (unoverridden) H2OAir fluid system uses
# Dumux::BinaryCoeff::H2O_Air::gasDiffCoeff(T, p) (dumux/material/binarycoefficients/h2o_air.hh),
# not a fixed constant, and the default TwoPNC/TwoPTwoC effective diffusivity model is
# DiffusivityMillingtonQuirk (dumux/material/fluidmatrixinteractions/diffusivitymillingtonquirk.hh),
# not a simple porosity*saturation*tortuosity scaling.
_D_AW_THETA = 1.8
_D_AW_REF = 2.13e-5  # m^2/s, reference value at p0/T0 below
_D_AW_P0 = 1.0e5  # Pa
_D_AW_T0 = 273.15  # K

# boundary state at x=0, matching HeatPipeProblem::dirichletAtPos
PG_BC = 1.013e5  # Pa
SW_BC = 0.99  # -
T_BC = 341.75  # K
HEAT_FLUX = -100.0  # W/m^2, matching Problem.HeatFlux (sign as used in the ODE below)
DOMAIN_LENGTH = 2.4  # m


def _p0_leverett():
    return np.sqrt(PHI / K)


def _capillary_pressure(se: float) -> float:
    p0 = _p0_leverett()
    return p0 * GAMMA * (1.417 * (1 - se) - 2.120 * (1 - se) ** 2 + 1.263 * (1 - se) ** 3)


def _dpc_dse(se: float, eps: float = 1e-8) -> float:
    return (_capillary_pressure(se + eps) - _capillary_pressure(se - eps)) / (2 * eps)


def _krl(se: float) -> float:
    return se**3


def _krg(se: float) -> float:
    return (1 - se) ** 3


def _vapor_pressure(T: float) -> float:
    """Dumux::Components::H2O::vaporPressure(T) (dumux/material/components/iapws/region4.hh,
    Region4::saturationPressure), the IAPWS-97 Region 4 saturation-pressure correlation
    (backward equation for p_sat(T))."""
    n = [
        0.11670521452767e4, -0.72421316703206e6, -0.17073846940092e2,
        0.12020824702470e5, -0.32325550322333e7, 0.14915108613530e2,
        -0.48232657361591e4, 0.40511340542057e6, -0.23855557567849,
        0.65017534844798e3,
    ]
    sigma = T + n[8] / (T - n[9])
    A = (sigma + n[0]) * sigma + n[1]
    B = (n[2] * sigma + n[3]) * sigma + n[4]
    C = (n[5] * sigma + n[6]) * sigma + n[7]
    term = 2.0 * C / (np.sqrt(B * B - 4.0 * A * C) - B)
    return 1e6 * term**4


def _binary_diffusion_coefficient(T: float, p: float) -> float:
    """Dumux::BinaryCoeff::H2O_Air::gasDiffCoeff(T, p)."""
    return _D_AW_REF * (_D_AW_P0 / p) * (T / _D_AW_T0) ** _D_AW_THETA


def _millington_quirk(sg: float, d_binary: float) -> float:
    """Dumux::DiffusivityMillingtonQuirk::effectiveDiffusionCoefficient."""
    sg = max(sg, 0.0)
    return PHI * sg**3 * (PHI * sg) ** (1.0 / 3.0) * d_binary


def _rhs(x, z):
    """Right-hand side of the four coupled ODEs (Huang et al. 2015, eq. F1-F4).

    State vector z = [Se, pg, Xa, T], with Se the effective wetting-phase
    saturation (Se = (Sw - Swr)/(1 - Swr)). Absolute saturations are used
    wherever a bulk volume fraction (thermal conductivity, gas-phase pore
    volume available for diffusion) is required; relative permeability and
    capillary pressure are, by construction, functions of Se.
    """
    se, pg, xa, T = z
    sw = SWR + se * (1 - SWR)
    sg = 1 - sw

    d_pm = _millington_quirk(sg, _binary_diffusion_coefficient(T, pg))
    lam = LAMBDA_DRY + np.sqrt(sw) * (LAMBDA_WET - LAMBDA_DRY)
    pc = _capillary_pressure(se)
    dpc_dse = _dpc_dse(se)

    krg = _krg(se)
    krl = _krl(se)
    rho_g_a = MA * xa * pg / R_GAS / T
    rho_g_w = MW * (1 - xa) * pg / R_GAS / T
    rho_g = rho_g_a + rho_g_w
    mu_g = xa * MU_G_A + (1 - xa) * MU_G_W
    nu_g = mu_g / rho_g
    nu_w = MU_W / RHO_W
    beta = nu_w / nu_g
    alpha = 1 + pc / RHO_W / H_WG
    xi = (1 / krg) * (1 + (RHO_W * R_GAS * T) / pg / MW * (1 / (1 - xa))) + beta / krl
    delta = RHO_W * H_WG * H_WG * K * alpha / (lam * nu_g * T)
    zeta = (
        ((K * RHO_W * R_GAS * T) / (MW * rho_g * nu_g * d_pm))
        * (xa / (1 - xa))
        * (pg * MW / RHO_W / R_GAS / T + 1 / (1 - xa))
    )
    eta = delta / (delta + xi + zeta)

    q = HEAT_FLUX
    dse_dx = -(1 / (1 - xa) + beta * krg / krl) * eta * q * nu_g / K / H_WG / krg / dpc_dse
    dpg_dx = -(eta * q * nu_g / K / H_WG / krg) * (1 / (1 - xa))
    dxa_dx = eta * q * xa / H_WG / d_pm / rho_g / (1 - xa)
    dT_dx = -q * (1 - eta) / lam
    return [dse_dx, dpg_dx, dxa_dx, dT_dx]


def _dryout_event(x, z):
    return z[0] - 1e-6


_dryout_event.terminal = True
_dryout_event.direction = -1


def semianalytical_solution(x: np.ndarray) -> dict:
    """Evaluate the semi-analytical heat-pipe solution on the sample points x.

    Returns a dict with arrays "Sw", "p_gas", "x_a_gas", "T" of the same
    shape as x, covering the full domain: the two-phase heat-pipe zone
    (obtained by integrating the ODE system from the left boundary until the
    wetting phase dries out) followed by a single-phase, purely-conductive
    dry zone extending to the end of the domain.
    """
    se_bc = (SW_BC - SWR) / (1 - SWR)
    # Equilibrium air mole fraction at the boundary, from equating the liquid- and
    # gas-phase fugacities of water (Dumux::FluidSystems::H2OAir::fugacityCoefficient):
    # phi_liq * x_liq^w * p_liq = phi_gas * x_gas^w * p_gas, with phi_liq = p_vap(T)/p_liq
    # and phi_gas = 1 (ideal gas). The liquid-phase pressure p_liq cancels out of
    # phi_liq * p_liq, so capillary pressure does not enter here (default
    # useKelvinVaporPressure = false); dissolved air in the liquid is neglected
    # (x_liq^w ~= 1), giving simply xa_bc = 1 - p_vap(T_bc) / pg_bc.
    xa_bc = 1 - _vapor_pressure(T_BC) / PG_BC

    sol = solve_ivp(
        _rhs,
        [0, DOMAIN_LENGTH],
        [se_bc, PG_BC, xa_bc, T_BC],
        method="RK45",
        max_step=5e-3,
        events=_dryout_event,
        dense_output=True,
        rtol=1e-8,
        atol=1e-10,
    )
    x_dry = sol.t[-1]

    sw = np.empty_like(x)
    pg = np.empty_like(x)
    xa = np.empty_like(x)
    T = np.empty_like(x)

    wet = x <= x_dry
    se_wet, pg[wet], xa[wet], T[wet] = sol.sol(x[wet])
    sw[wet] = SWR + se_wet * (1 - SWR)

    # Beyond the two-phase zone the wetting phase has fully evaporated
    # (Se -> 0 is where kr_w -> 0, i.e. the wetting phase becomes immobile;
    # the numerical model then switches to a single, gas-only phase state
    # with Sw = 0 exactly, not Sw = Swr).
    dry = ~wet
    _, pg_dry, _, T_dry = sol.sol(x_dry)
    sw[dry] = 0.0
    xa[dry] = 0.0
    pg[dry] = pg_dry
    T[dry] = T_dry + (-HEAT_FLUX / LAMBDA_DRY) * (x[dry] - x_dry)

    return {"Sw": sw, "p_gas": pg, "x_a_gas": xa, "T": T, "heatpipe_length": x_dry}


# ---------------------------------------------------------------------------
# Post-processing / plotting
# ---------------------------------------------------------------------------


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


def point_data(mesh, field: str):
    # the box method (used by this test) writes results as point data
    if field not in mesh.point_data:
        raise KeyError(f"Missing point field '{field}'")
    x = mesh.points[:, 0]
    order = np.argsort(x)
    return x[order], mesh.point_data[field][order]


def _comparison_data(vtu: Path):
    pv, _ = require_plot_modules()
    mesh = pv.read(vtu)
    x_max = mesh.points[:, 0].max()
    x_analytic = np.linspace(0, x_max, 2000)
    analytic = semianalytical_solution(x_analytic)
    print(f"semi-analytical heat-pipe length: {analytic['heatpipe_length']:.4f} m")
    return mesh, x_max, x_analytic, analytic


def create_saturation_lineplot(vtu: Path, image_file: Path) -> None:
    """Single-panel wetting-phase saturation comparison (used as the benchmark thumbnail)."""
    _, plt = require_plot_modules()
    mesh, x_max, x_analytic, analytic = _comparison_data(vtu)

    fig, ax = plt.subplots(figsize=(8.2, 4.8), constrained_layout=True)
    x, sw = point_data(mesh, SATURATION_FIELD)
    ax.plot(x, sw, label="numerical", linewidth=2)
    ax.plot(x_analytic, analytic["Sw"], "k--", linewidth=2.4, label="semi-analytical")
    ax.set_xlabel("x [m]")
    ax.set_ylabel(r"Wetting-phase saturation $S_w$ [-]")
    ax.set_xlim(0, x_max)
    ax.set_ylim(-0.02, 1.02)
    ax.grid(True, alpha=0.3)
    ax.legend()
    fig.savefig(image_file, dpi=200)
    plt.close(fig)


def create_line_plot(vtu: Path, image_file: Path) -> None:
    _, plt = require_plot_modules()
    mesh, x_max, x_analytic, analytic = _comparison_data(vtu)

    fig, axes = plt.subplots(2, 2, figsize=(11, 7.5), constrained_layout=True)
    ax_sw, ax_T, ax_pg, ax_xa = axes[0, 0], axes[0, 1], axes[1, 0], axes[1, 1]

    x, sw = point_data(mesh, SATURATION_FIELD)
    ax_sw.plot(x, sw, label="numerical", linewidth=1.8)

    _, T = point_data(mesh, TEMPERATURE_FIELD)
    ax_T.plot(x, T, label="numerical", linewidth=1.8)

    _, pg = point_data(mesh, PRESSURE_FIELD)
    ax_pg.plot(x, pg, label="numerical", linewidth=1.8)

    _, xa = point_data(mesh, AIR_MOLEFRACTION_FIELD)
    ax_xa.plot(x, np.clip(xa, 0, None), label="numerical", linewidth=1.8)

    ax_sw.plot(x_analytic, analytic["Sw"], "k--", linewidth=2, label="semi-analytical")
    ax_T.plot(x_analytic, analytic["T"], "k--", linewidth=2, label="semi-analytical")
    ax_pg.plot(x_analytic, analytic["p_gas"], "k--", linewidth=2, label="semi-analytical")
    ax_xa.plot(x_analytic, analytic["x_a_gas"], "k--", linewidth=2, label="semi-analytical")

    ax_sw.set_ylabel(r"Wetting-phase saturation $S_w$ [-]")
    ax_T.set_ylabel("Temperature $T$ [K]")
    ax_pg.set_ylabel("Gas-phase pressure $p_g$ [Pa]")
    ax_xa.set_ylabel(r"Gas-phase air mole fraction $x_g^a$ [-]")
    for ax in (ax_sw, ax_T, ax_pg, ax_xa):
        ax.set_xlabel("x [m]")
        ax.set_xlim(0, x_max)
        ax.grid(True, alpha=0.3)
        ax.legend(fontsize=8)

    fig.savefig(image_file, dpi=200)
    plt.close(fig)


def create_saturation_image(vtu: Path, image_file: Path) -> None:
    pv, _ = require_plot_modules()
    mesh = pv.read(vtu)
    if SATURATION_FIELD not in mesh.array_names:
        raise KeyError(f"Missing field '{SATURATION_FIELD}' in {vtu}")

    plotter = pv.Plotter(
        off_screen=True, window_size=(1000, 260), border=False, border_color="white"
    )
    plotter.set_background("white")
    plotter.add_mesh(
        mesh,
        scalars=SATURATION_FIELD,
        show_edges=False,
        cmap="coolwarm",
        scalar_bar_args={
            "title": "S_w",
            "vertical": True,
            "position_x": 0.90,
            "position_y": 0.15,
            "width": 0.06,
            "height": 0.70,
        },
    )
    plotter.view_xy()
    plotter.camera.zoom(1.3)
    plotter.show(screenshot=str(image_file), auto_close=True)


def create_grid_convergence_plot(vtus: dict[int, Path], image_file: Path) -> None:
    """Wetting-phase saturation and temperature at each grid resolution, zoomed
    to the dry-out front where the resolution-dependent mismatch actually shows
    up (gas-phase pressure and air mole fraction barely change with resolution,
    so they are not repeated here)."""
    _, plt = require_plot_modules()
    pv, _ = require_plot_modules()

    finest = max(vtus)
    mesh_finest = pv.read(vtus[finest])
    x_max = mesh_finest.points[:, 0].max()
    x_analytic_full = np.linspace(0, x_max, 4000)
    x_zoom_min = max(0.0, semianalytical_solution(x_analytic_full)["heatpipe_length"] - 0.6)

    # Restrict all plotted data to the zoomed window up front, rather than plotting
    # the full domain and relying on `ax.set_xlim` to crop it: matplotlib's SVG
    # backend draws the un-cropped line past the axis bounds and hides the excess
    # via an SVG clip-path, which some SVG renderers (e.g. GitLab's markdown
    # sanitizer) strip, making the "hidden" part of the line visible again.
    x_analytic = np.linspace(x_zoom_min, x_max, 2000)
    analytic = semianalytical_solution(x_analytic)

    fig, (ax_sw, ax_T) = plt.subplots(1, 2, figsize=(11, 4.8), constrained_layout=True)
    colors = plt.cm.viridis(np.linspace(0.15, 0.85, len(vtus)))
    for color, cells in zip(colors, sorted(vtus)):
        mesh = pv.read(vtus[cells])
        x, sw = point_data(mesh, SATURATION_FIELD)
        _, T = point_data(mesh, TEMPERATURE_FIELD)
        zoom = x >= x_zoom_min
        ax_sw.plot(x[zoom], sw[zoom], color=color, linewidth=1.6, label=f"{cells} cells")
        ax_T.plot(x[zoom], T[zoom], color=color, linewidth=1.6, label=f"{cells} cells")

    ax_sw.plot(x_analytic, analytic["Sw"], "k--", linewidth=2, label="semi-analytical")
    ax_T.plot(x_analytic, analytic["T"], "k--", linewidth=2, label="semi-analytical")

    ax_sw.set_xlabel("x [m]")
    ax_sw.set_ylabel(r"Wetting-phase saturation $S_w$ [-]")
    ax_T.set_xlabel("x [m]")
    ax_T.set_ylabel("Temperature $T$ [K]")
    for ax in (ax_sw, ax_T):
        ax.set_xlim(x_zoom_min, x_max)
        ax.grid(True, alpha=0.3)
        ax.legend(fontsize=9)

    fig.savefig(image_file, dpi=200)
    plt.close(fig)


def main() -> None:
    vtus = build_and_run_all()
    finest = vtus[max(vtus)]
    out_dir = case_build_dir()

    print("Creating line plot (finest resolution)...")
    create_line_plot(finest, out_dir / f"{BASE_NAME}_lineplot_comparison.svg")
    print("Creating saturation line plot (finest resolution)...")
    create_saturation_lineplot(finest, out_dir / f"{BASE_NAME}_saturation_comparison.svg")
    print("Creating saturation field image (finest resolution)...")
    create_saturation_image(finest, out_dir / f"{BASE_NAME}_sw.png")
    print("Creating grid convergence plot...")
    create_grid_convergence_plot(vtus, out_dir / f"{BASE_NAME}_grid_convergence.svg")


if __name__ == "__main__":
    main()
