#!/usr/bin/env python3
# SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

"""
Linearized one-sided convection problem (Riaz et al. 2006, eqs. 2.14-2.17), the reference of the benchmark.

Lengths in units of the height H, depth eta below the saturated top, time in units of phi H / U.
A perturbation c(eta, t) cos(kx), w(eta, t) cos(kx) of the diffusive profile C0 = erfc(eta / (2 sqrt(t / Ra))) satisfies

    (d^2/deta^2 - k^2) w = -k^2 c,                      w(0) = 0,
    dc/dt = (d^2/deta^2 - k^2) c / Ra - dC0/deta w,     c(0) = 0,

where w is the downward velocity. Both decay for eta -> infinity: the layer is semi-infinite, which is exact for the
benchmark as the boundary layer is thin compared with the height (erfc(1 / (2 sqrt(t / Ra))) < 1e-28 for t <= 2, Ra = 500).

The problem is discretized with the Galerkin method and the Laguerre functions of the variable xi = eta / s,

    phi_m(eta) = exp(-xi / 2) (L_{m+1}(xi) - L_m(xi)),   m = 0, ..., n - 1,

for c and w. They vanish at eta = 0 and decay for large eta, so the boundary conditions are built in. Because of the
orthogonality of the Laguerre polynomials, the mass and stiffness matrices are tridiagonal. In time, the system is
integrated with the Radau IIA method.

Usage: python3 linear_stability.py   (prints the checks of the solver and the reference values of the tests)
"""

import numpy as np
import sympy as sp
from scipy.integrate import solve_ivp
from scipy.linalg import solve

# base state, with its gradient dC0/deta = -sqrt(Ra / (pi t)) exp(-Ra eta^2 / (4 t)) derived symbolically
eta_, t_, ra_ = sp.symbols("eta t Ra", positive=True)
C0 = sp.erfc(eta_ / (2 * sp.sqrt(t_ / ra_)))
base_gradient = sp.lambdify((eta_, t_, ra_), sp.diff(C0, eta_), "numpy")


def tridiagonal(n, diagonal, off_diagonal):
    return diagonal * np.eye(n) + off_diagonal * (np.eye(n, k=1) + np.eye(n, k=-1))


class LinearizedProblem:
    def __init__(self, ra, k, t0, t1, n=100):
        """The perturbation at the time t0, evolved until the time t1, with n Laguerre functions"""
        self.ra, self.k, self.n, self.t0, self.t1 = ra, k, n, t0, t1
        # The functions reach about 4 n s into the layer, so s has to be large enough for the longest length of the
        # problem (the perturbation at t1 and the decay of w with 1 / k) and small enough for the thinnest boundary layer.
        longest = max(12.0 * np.sqrt(t1 / ra), 8.0 / k)
        self.s = max(2.0 * np.sqrt(t0 / ra), longest / (4.0 * n))
        # int exp(-xi) L_i L_j dxi = delta_ij gives the integrals of the functions exactly (checked in main)
        self.mass = self.s * tridiagonal(n, 2.0, -1.0)  # int phi_i phi_j deta
        stiffness_k = tridiagonal(n, 2.0, 1.0) / (4.0 * self.s) + k**2 * self.mass  # int phi_i' phi_j' + k^2 phi_i phi_j
        self.velocity = k**2 * solve(stiffness_k, self.mass)  # coefficients of w from those of c
        self.diffusion = solve(self.mass, stiffness_k) / ra

    def basis(self, eta):
        """The functions phi_m at the depths eta, shape (n, len(eta)), by the recurrence of the Laguerre polynomials"""
        xi = np.asarray(eta) / self.s
        laguerre = np.empty((self.n + 1, len(xi)))
        laguerre[0] = 1.0
        laguerre[1] = 1.0 - xi
        for j in range(1, self.n):
            laguerre[j + 1] = ((2 * j + 1 - xi) * laguerre[j] - j * laguerre[j - 1]) / (j + 1)
        return np.exp(-xi / 2) * (laguerre[1:] - laguerre[:-1])

    def solve(self, perturbation, times, rtol=1e-8):
        """Coefficients of the perturbation c(eta, t) at the times t0 <= times <= t1, one row per time, with
        c(eta, t0) = perturbation(eta)"""
        # Gauss-Legendre rule for the integrals with nodes eta = s u^2 clustered at the top, up to the depth where the gradient
        # of the base state (12 sqrt(t1 / Ra)) and the initial perturbation (13 sqrt(t0 / Ra)) are below 1e-15 of their maximum
        u, w = np.polynomial.legendre.leggauss(4 * self.n)
        u_max = np.sqrt(max(12.0 * np.sqrt(self.t1 / self.ra), 13.0 * np.sqrt(self.t0 / self.ra)) / self.s)
        u, w = 0.5 * u_max * (u + 1), 0.5 * u_max * w
        eta, w = self.s * u**2, 2 * self.s * u * w
        phi = self.basis(eta)
        mass_inverse = np.linalg.inv(self.mass)

        def operator(t):
            coupling = (phi * (w * base_gradient(eta, t, self.ra))) @ phi.T  # int phi_i dC0/deta phi_j deta
            return -self.diffusion - mass_inverse @ coupling @ self.velocity

        c0 = solve(self.mass, (phi * w) @ perturbation(eta))
        # absolute tolerance relative to the largest coefficient: small coefficients need not be accurate to rtol
        solution = solve_ivp(lambda t, c: operator(t) @ c, (self.t0, times[-1]), c0, method="Radau", jac=lambda t, c: operator(t),
                             t_eval=times, rtol=rtol, atol=1e-2 * rtol * np.abs(c0).max())
        return solution.y.T

    def l2_norm(self, c):
        """sqrt(int c^2 deta) for the coefficients c (one row per time)"""
        c = np.atleast_2d(c)
        return np.sqrt(np.einsum("ti,ij,tj->t", c, self.mass, c))

    def profile(self, c, eta):
        """The perturbation c(eta) of the coefficients c (one time)"""
        return c @ self.basis(eta)

    def max_norm(self, c, eta_max=0.5):
        """max over the depth of |c| for the coefficients c (one row per time)"""
        return np.abs(np.atleast_2d(c) @ self.basis(np.linspace(0.0, eta_max, 4001))).max(axis=1)


def dominant_mode(ra, t0):
    """Initial perturbation zeta exp(-zeta^2), zeta = eta sqrt(Ra / (4 t0)) (Riaz et al. 2006)"""
    def perturbation(eta):
        zeta = eta * np.sqrt(ra / (4.0 * t0))
        return zeta * np.exp(-zeta**2)
    return perturbation


def main():
    ra, t0 = 500.0, 0.01

    def log_amplification(ra, k, t0, t1, n=100):
        """ln(A(t1) / A(t0)) of the L2 norm A of the dominant-mode perturbation"""
        problem = LinearizedProblem(ra, k, t0, t1, n)
        norm = problem.l2_norm(problem.solve(dominant_mode(ra, t0), np.array([t0, t1])))
        return np.log(norm[1] / norm[0])

    # the matrices are exact: the integrals of the functions with s = 1 are tridiagonal (verified symbolically for 5)
    xi = sp.symbols("xi", positive=True)
    phi = [sp.exp(-xi / 2) * (sp.laguerre(m + 1, xi) - sp.laguerre(m, xi)) for m in range(5)]
    mass = sp.Matrix(5, 5, lambda i, j: sp.integrate(sp.simplify(phi[i] * phi[j]), (xi, 0, sp.oo)))
    stiffness = sp.Matrix(5, 5, lambda i, j: sp.integrate(sp.simplify(sp.diff(phi[i], xi) * sp.diff(phi[j], xi)), (xi, 0, sp.oo)))
    exact = LinearizedProblem(1.0, 1.0, 0.25, 0.25 + 1e-9, 5)  # s = 1
    print("mass matrix, deviation from the symbolic integrals:", np.abs(np.array(mass, dtype=float) - exact.mass).max())
    print("stiffness matrix, deviation:", np.abs(np.array(stiffness, dtype=float) - tridiagonal(5, 2.0, 1.0) / 4.0).max())

    # k = 0: the perturbation eta t^(-3/2) exp(-Ra eta^2 / (4 t)) is an exact solution of the heat equation
    c = eta_ * t_**sp.Rational(-3, 2) * sp.exp(-ra_ * eta_**2 / (4 * t_))
    print("k = 0: residual of the exact solution (symbolic):", sp.simplify(sp.diff(c, t_) - sp.diff(c, eta_, 2) / ra_))

    # reference value of the test: ln of the amplification of the L2 norm, Ra = 500, k = 8 pi, t = 0.01 to 2
    for n in (80, 100, 120):
        print(f"Ra = 500, k = 8 pi, {n} functions: ln(A(2) / A(0.01)) = {log_amplification(ra, 8 * np.pi, t0, 2.0, n):.8f}")

    # growth rate of max|c| for the dominant-mode perturbation, Riaz et al. (2006), figure 4
    times = np.linspace(t0, 4.0, 800)
    for k in (4 * np.pi, 8 * np.pi):
        problem = LinearizedProblem(ra, k, t0, 4.0)
        sigma = np.gradient(np.log(problem.max_norm(problem.solve(dominant_mode(ra, t0), times))), times)
        print(f"Ra = {ra:g}, k = {k / np.pi:g} pi: growth rate of max|c| at t = 1, 2, 4: "
              + ", ".join(f"{np.interp(t, times, sigma):.2f}" for t in (1, 2, 4)))

    # Ra only sets the scales: with k ~ Ra and the times ~ 1 / Ra, the amplification is the same
    for ra in (500.0, 5e3, 5e4, 5e5):
        s = ra / 500.0
        print(f"Ra = {ra:g}, k = {8 * s:g} pi, t = {t0 / s:g} to {2.0 / s:g}: ln(A(t1) / A(t0)) = "
              f"{log_amplification(ra, 8 * np.pi * s, t0 / s, 2.0 / s):.8f}")


if __name__ == "__main__":
    main()
