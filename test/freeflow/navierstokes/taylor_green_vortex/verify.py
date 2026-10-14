#!/usr/bin/env python3
# SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later
"""
Verify symbolically (with sympy) the analytical solutions of the Taylor-Green vortex benchmark
and the formulas given in README.md:

  * 2D: Taylor-Green vortex (Taylor 1923) with amplitude U_0 and wavelength L
  * 3D: Beltrami (ABC) flow (Barbato, Berselli & Grisanti 2007) with amplitudes A, B, C
    and wavelength L

For each solution, the following is checked (with mu = rho nu and k = 2 pi / L):

  * the velocity is divergence free,
  * the instationary Navier-Stokes equations without source term are satisfied,
    rho du/dt + div(rho u (x) u) - div(mu grad u) + grad p = 0,
  * the identities behind the source terms of the stationary and the Stokes variants,
    rho du/dt = div(mu grad u) = -lambda mu u and div(rho u (x) u) = -grad p,
  * the vorticity (2D: the scalar dv/dx - du/dy, 3D: the vector curl u) is given by the stated
    formula,
  * the kinetic energy E = rho/2 int |u|^2 and the enstrophy eps = rho/2 int |omega|^2 over one
    wavelength in each direction are given by the stated closed forms,
  * the energy balance -dE/dt = 2 nu eps holds.
"""

import argparse
import sys
from collections import namedtuple
from dataclasses import dataclass

import sympy as sp

# the symbols shared by both solutions: time, kinematic viscosity, density and wavelength
Symbols = namedtuple("Symbols", "t nu rho wavelength")


@dataclass
class Solution:
    """An analytical solution together with the stated formulas derived from it"""

    name: str
    coords: list
    velocity: list
    pressure: sp.Expr
    vorticity: list  # the components of the vorticity (one in 2D)
    eigenvalue: sp.Expr  # lambda with laplace(u) = -lambda u
    indicators: tuple  # the closed forms of the kinetic energy and the enstrophy


def makeSymbols():
    """The symbols shared by both solutions"""
    return Symbols(*sp.symbols("t nu rho L", positive=True))


def report(label, expression):
    """Print whether the (simplified) expression vanishes and return the result"""
    passed = sp.simplify(expression) == 0
    print(f"  {'passed' if passed else 'FAILED'}: {label}")
    return passed


def divergence(vector, coords):
    """The divergence of a vector field"""
    return sum(sp.diff(component, x) for component, x in zip(vector, coords))


def laplacian(scalar, coords):
    """The Laplacian of a scalar field"""
    return sum(sp.diff(scalar, x, 2) for x in coords)


def curl(velocity, coords):
    """The components of the curl (one in 2D)"""
    if len(coords) == 2:
        return [sp.diff(velocity[1], coords[0]) - sp.diff(velocity[0], coords[1])]

    def derivative(component, direction):
        return sp.diff(velocity[component % 3], coords[direction % 3])

    return [derivative(i + 2, i + 1) - derivative(i + 1, i + 2) for i in range(3)]


def integrateOverPeriod(expression, coords, wavelength):
    """Integral over one wavelength [0, L] in each direction"""
    for x in coords:
        expression = sp.integrate(sp.expand(expression), (x, 0, wavelength))
    return expression


def checkMomentum(solution, symbols):
    """Check the mass and momentum balance and the identities behind the source terms"""
    coords, velocity = solution.coords, solution.velocity
    rho, mu = symbols.rho, symbols.rho * symbols.nu
    success = report("divergence free", divergence(velocity, coords))
    for i, (u, x) in enumerate(zip(velocity, coords)):
        storage = rho * sp.diff(u, symbols.t)
        viscous = mu * laplacian(u, coords)
        convection = rho * divergence([u * v for v in velocity], coords)
        pressureGradient = sp.diff(solution.pressure, x)
        success &= report(
            f"Navier-Stokes momentum balance (component {i})",
            storage + convection - viscous + pressureGradient,
        )
        success &= report(f"rho du/dt = div(mu grad u) (component {i})", storage - viscous)
        success &= report(
            f"stationary source -rho du/dt = lambda mu u (component {i})",
            -storage - solution.eigenvalue * mu * u,
        )
        success &= report(
            f"Stokes source -div(rho u (x) u) = grad p (component {i})",
            -convection - pressureGradient,
        )
    return success


def checkIndicators(solution, symbols):
    """Check the vorticity, the kinetic energy, the enstrophy and the energy balance"""
    coords, velocity = solution.coords, solution.velocity
    vorticity = curl(velocity, coords)
    success = True
    for i, (component, expected) in enumerate(zip(vorticity, solution.vorticity)):
        success &= report(f"vorticity (component {i})", component - expected)

    energy, enstrophy = solution.indicators
    rho, wavelength = symbols.rho, symbols.wavelength
    kineticEnergy = integrateOverPeriod(rho / 2 * sum(u**2 for u in velocity), coords, wavelength)
    integratedEnstrophy = integrateOverPeriod(
        rho / 2 * sum(w**2 for w in vorticity), coords, wavelength
    )
    success &= report(f"kinetic energy E = {energy}", kineticEnergy - energy)
    success &= report(f"enstrophy eps = {enstrophy}", integratedEnstrophy - enstrophy)
    success &= report(
        "energy balance -dE/dt = 2 nu eps",
        -sp.diff(energy, symbols.t) - 2 * symbols.nu * enstrophy,
    )
    return success


def verify2d():
    """The two-dimensional Taylor-Green vortex (Taylor 1923)"""
    s = makeSymbols()
    x, y = sp.symbols("x y", real=True)
    u0 = sp.symbols("U_0", positive=True)
    k = 2 * sp.pi / s.wavelength
    decay = sp.exp(-2 * s.nu * k**2 * s.t)
    energy = s.rho * u0**2 * s.wavelength**2 / 4 * decay**2
    solution = Solution(
        name="2D: Taylor-Green vortex (Taylor 1923)",
        coords=[x, y],
        velocity=[
            -u0 * sp.cos(k * x) * sp.sin(k * y) * decay,
            u0 * sp.sin(k * x) * sp.cos(k * y) * decay,
        ],
        pressure=-s.rho * u0**2 / 4 * (sp.cos(2 * k * x) + sp.cos(2 * k * y)) * decay**2,
        vorticity=[2 * k * u0 * sp.cos(k * x) * sp.cos(k * y) * decay],
        eigenvalue=2 * k**2,
        indicators=(energy, 2 * k**2 * energy),
    )
    print(solution.name)
    return checkMomentum(solution, s) & checkIndicators(solution, s)


def verify3d():
    """The three-dimensional Beltrami (ABC) flow (Barbato, Berselli & Grisanti 2007)"""
    s = makeSymbols()
    x, y, z, a, b, c = sp.symbols("x y z A B C", real=True)
    k = 2 * sp.pi / s.wavelength
    decay = sp.exp(-s.nu * k**2 * s.t)
    velocity = [
        (a * sp.sin(k * z) + c * sp.cos(k * y)) * decay,
        (b * sp.sin(k * x) + a * sp.cos(k * z)) * decay,
        (c * sp.sin(k * y) + b * sp.cos(k * x)) * decay,
    ]
    energy = s.rho * (a**2 + b**2 + c**2) * s.wavelength**3 / 2 * decay**2
    bernoulliTerms = sum(
        (
            b * c * sp.cos(k * x) * sp.sin(k * y),
            a * b * sp.sin(k * x) * sp.cos(k * z),
            a * c * sp.sin(k * z) * sp.cos(k * y),
        )
    )
    solution = Solution(
        name="3D: Beltrami (ABC) flow (Barbato, Berselli & Grisanti 2007)",
        coords=[x, y, z],
        velocity=velocity,
        pressure=-s.rho * bernoulliTerms * decay**2,
        vorticity=[k * u for u in velocity],
        eigenvalue=k**2,
        indicators=(energy, k**2 * energy),
    )
    print(solution.name)
    return checkMomentum(solution, s) & checkIndicators(solution, s)


def main():
    """Run the checks of the selected dimensions"""
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("--dim", type=int, nargs="+", choices=(2, 3), default=[2, 3])
    args = parser.parse_args()

    success = True
    if 2 in args.dim:
        success &= verify2d()
    if 3 in args.dim:
        success &= verify3d()

    print("All checks passed" if success else "Some checks FAILED")
    sys.exit(0 if success else 1)


if __name__ == "__main__":
    main()
