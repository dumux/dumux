# SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later
"""
Ramey (1962) analytical solution for the wellbore heat transport benchmark.

The overall heat transfer coefficient U and the characteristic length X are both
referred to the inner pipe radius, see README.md.  Only laminar pipe flow is
supported.

Usage::

    from ramey import RameySolution
    ramey = RameySolution(injection_temperature=20.0, formation_temperature=55.0,
                          injection_rate=2e-4, length=30.0,
                          r_inner=0.12913, r_outer=0.14, lambda_pg=0.962,
                          rho_f=1000.0, c_p_f=4190.0, mu_f=1.14e-3, lambda_f=0.59,
                          lambda_s=2.78018, rho_s=1800.0, c_p_s=1778.0)
    ramey(z, t)          # temperature along the borehole at time t
    ramey.outlet(times)  # outlet temperature over time
"""
import numpy as np

# ── Laminar pipe flow ─────────────────────────────────────────────────────────
RE_LAMINAR = 2300.0   # upper limit of laminar pipe flow
NU_LAMINAR = 4.364    # Nusselt number of laminar pipe flow with constant wall heat flux


def outlet_temp(T_s, T_i, z, X):
    """Fluid temperature at distance z along the borehole."""
    return T_s + (T_i - T_s) * np.exp(-z / X)


def coefficient_x(q, rho_f, c_p_f, lambda_s, r_ref, U, f_t):
    """Characteristic length X of the exponential temperature decay."""
    return (q * rho_f * c_p_f) * (lambda_s + r_ref * U * f_t) / (2 * np.pi * r_ref * U * lambda_s)


def dimensionless_time(lambda_s, delta_t, rho_s, c_p_s, r_b):
    return lambda_s * delta_t / (rho_s * c_p_s * r_b**2)


def time_function(t_D):
    if t_D > 1.5:
        return (0.4063 + 0.5 * np.log(t_D)) * (1 + 0.6 / t_D)
    return 1.1281 * np.sqrt(t_D) * (1 - 0.3 * np.sqrt(t_D))


class RameySolution:
    """Ramey solution T(z, t) in °C for a borehole with laminar pipe flow.

    injection_temperature -- inlet fluid temperature, °C
    formation_temperature -- undisturbed formation temperature, °C
    injection_rate        -- volumetric flow rate, m³/s
    length                -- borehole length, m
    r_inner, r_outer      -- inner pipe radius and borehole radius, m
    lambda_pg             -- combined thermal conductivity of pipe wall and grout, W/m/K
    rho_f, c_p_f, mu_f, lambda_f -- density (kg/m³), specific heat (J/kg/K), dynamic
                             viscosity (Pa·s) and thermal conductivity (W/m/K) of the fluid
    lambda_s, rho_s, c_p_s -- thermal conductivity (W/m/K), density (kg/m³)
                             and specific heat (J/kg/K) of the formation
    """

    def __init__(self, injection_temperature, formation_temperature, injection_rate, length,
                 r_inner, r_outer, lambda_pg, rho_f, c_p_f, mu_f, lambda_f, lambda_s, rho_s, c_p_s):
        self.T_i = injection_temperature
        self.T_s = formation_temperature
        self.q = injection_rate
        self.length = length
        self.r_pi = r_inner
        self.r_b = r_outer
        self.rho_f = rho_f
        self.c_p_f = c_p_f
        self.lambda_s = lambda_s
        self.rho_s = rho_s
        self.c_p_s = c_p_s

        v = self.q / (np.pi * self.r_pi**2)
        self.Pr = mu_f * c_p_f / lambda_f
        self.Re = rho_f * v * (2 * self.r_pi) / mu_f
        if self.Re >= RE_LAMINAR:
            raise ValueError(f"Re = {self.Re:.0f} is not laminar (Re < {RE_LAMINAR:.0f}), "
                             "only the laminar Nusselt number is implemented")
        self.Nu = NU_LAMINAR
        self.h = lambda_f * self.Nu / (2 * self.r_pi)

        # overall heat transfer coefficient, referred to the inner pipe radius
        self.U = 1 / (1 / self.h + self.r_pi * np.log(self.r_b / self.r_pi) / lambda_pg)

    def __call__(self, z, t):
        """Temperature in °C at distance z [m] along the borehole and time t [s]."""
        z = np.asarray(z, dtype=float)
        f_t = time_function(dimensionless_time(self.lambda_s, t, self.rho_s, self.c_p_s, self.r_b))
        X = coefficient_x(self.q, self.rho_f, self.c_p_f, self.lambda_s, self.r_pi, self.U, f_t)
        return outlet_temp(self.T_s, self.T_i, z, X)

    def outlet(self, times):
        """Outlet temperature in °C for an array of times [s]."""
        return np.array([float(self(self.length, t)) for t in times])
