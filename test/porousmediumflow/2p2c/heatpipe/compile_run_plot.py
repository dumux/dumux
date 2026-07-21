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

SATURATION_FIELD = "S_liq"
PRESSURE_FIELD = "p_gas"
TEMPERATURE_FIELD = "T"
AIR_MOLEFRACTION_FIELD = "x^Air_gas"


def root_dir() -> Path:
    return Path(__file__).resolve().parents[4]


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


def build_and_run() -> Path:
    build_dir = root_dir() / "build-cmake"
    run(["cmake", "--build", ".", "--target", TARGET], cwd=build_dir)

    remove_old_outputs(BASE_NAME)
    run([str(case_build_dir() / TARGET), "heatpipe.input"], cwd=case_build_dir())
    return latest_vtu(BASE_NAME)


# ---------------------------------------------------------------------------
# Semi-analytical reference solution (Udell & Fitch 1985 / Huang et al. 2015)
# ---------------------------------------------------------------------------

# physical parameters, matching heatpipespatialparams.hh, krpcheatpipe.hh,
# heatpipeproblem.hh and heatpipe.input
K = 1.0e-12  # m^2,  Problem.Permeability in heatpipe.input
PHI = 0.4  # -,    HeatPipeSpatialParams::porosity_
SWR = 0.15  # -,    swr passed to KrPcHeatPipe::Params in heatpipespatialparams.hh
LAMBDA_DRY = 0.582  # W/(m*K)
LAMBDA_WET = 1.13  # W/(m*K)
RHO_W = 958.4  # kg/m^3, water density near the phase-change reference state
MU_W = 2.938e-4  # Pa*s
MU_G_A = 2.08e-5  # Pa*s, dynamic viscosity of air
MU_G_W = 1.2e-5  # Pa*s, dynamic viscosity of steam
D_BINARY = 2.6e-6  # m^2/s, binary diffusion coefficient air-water vapor
TAU = 0.5  # -, tortuosity
GAMMA = 0.05878  # N/m, surface tension of water at ~100.5°C
MW = 0.018016  # kg/mol, molar mass of water
MA = 0.02897  # kg/mol, molar mass of air
R_GAS = 8.3144621  # J/(mol*K)
H_WG = 2.258e6  # J/kg, latent heat of vaporization of water
P0_REF = 101325.0  # Pa, reference pressure for the Clausius-Clapeyron relation
T0_REF = 373.15  # K, reference (boiling) temperature for the Clausius-Clapeyron relation

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

    d_pm = PHI * sg * TAU * D_BINARY
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
    pc_bc = _capillary_pressure(se_bc)
    xa_bc = 1 - P0_REF / PG_BC * np.exp(
        (1 / T0_REF - 1 / T_BC) * H_WG * MW / R_GAS - pc_bc * MW / RHO_W / R_GAS / T_BC
    )

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

    se = np.empty_like(x)
    pg = np.empty_like(x)
    xa = np.empty_like(x)
    T = np.empty_like(x)

    wet = x <= x_dry
    wet_vals = sol.sol(x[wet])
    se[wet], pg[wet], xa[wet], T[wet] = wet_vals

    dry = ~wet
    se_dry, pg_dry, xa_dry, T_dry = sol.sol(x_dry)
    se[dry] = 0.0
    xa[dry] = 0.0
    pg[dry] = pg_dry
    T[dry] = T_dry + (-HEAT_FLUX / LAMBDA_DRY) * (x[dry] - x_dry)

    sw = SWR + se * (1 - SWR)
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


def create_line_plot(vtu: Path, image_file: Path) -> None:
    pv, plt = require_plot_modules()

    fig, axes = plt.subplots(2, 2, figsize=(11, 7.5), constrained_layout=True)
    ax_sw, ax_T, ax_pg, ax_xa = axes[0, 0], axes[0, 1], axes[1, 0], axes[1, 1]

    mesh = pv.read(vtu)
    x_max = mesh.points[:, 0].max()
    x_analytic = np.linspace(0, x_max, 2000)
    analytic = semianalytical_solution(x_analytic)
    print(f"semi-analytical heat-pipe length: {analytic['heatpipe_length']:.4f} m")

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


def main() -> None:
    vtu = build_and_run()
    out_dir = case_build_dir()

    print("Creating line plot...")
    create_line_plot(vtu, out_dir / f"{BASE_NAME}_lineplot_comparison.png")
    print("Creating saturation field image...")
    create_saturation_image(vtu, out_dir / f"{BASE_NAME}_sw.png")


if __name__ == "__main__":
    main()
