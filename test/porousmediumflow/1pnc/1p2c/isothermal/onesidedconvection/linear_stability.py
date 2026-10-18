#!/usr/bin/env python3
# SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

"""
Linearized one-sided convection problem (Riaz et al. 2006, eqs. 2.14-2.17), the reference of the benchmark.

Lengths in units of the height H, depth z below the saturated top, time in units of phi H / U.
A perturbation c(z, t) cos(kx), w(z, t) cos(kx) of the diffusive profile C0(z, t) satisfies

    (d^2/dz^2 - k^2) w = -k^2 c,                   w(0) = w(1) = 0,
    dc/dt = (d^2/dz^2 - k^2) c / Ra - dC0/dz w,    c(0) = 0, dc/dz(1) = 0,

where w is the downward velocity. The problem is discretized with Chebyshev collocation in z and
integrated in time with the Radau IIA method; C0 is the diffusive profile of the layer of finite depth.

Usage: python3 linear_stability.py   (prints the checks of the solver and the reference values of the tests)
"""

import numpy as np
from scipy.integrate import solve_ivp
from scipy.interpolate import BarycentricInterpolator
from scipy.special import erfc


def chebyshev(n):
    """Chebyshev points x_j = cos(pi j / n) and the differentiation matrix"""
    x = np.cos(np.pi * np.arange(n + 1) / n)
    c = np.hstack([2.0, np.ones(n - 1), 2.0]) * (-1.0) ** np.arange(n + 1)
    d = np.outer(c, 1.0 / c) / (x[:, None] - x[None, :] + np.eye(n + 1))
    return d - np.diag(d.sum(axis=1)), x


def clenshaw_curtis_weights(n):
    """Quadrature weights on [-1, 1] at the Chebyshev points, n even"""
    theta = np.pi * np.arange(n + 1) / n
    w = np.zeros(n + 1)
    v = np.ones(n - 1)
    w[0] = w[n] = 1.0 / (n**2 - 1)
    for j in range(1, n // 2):
        v -= 2.0 * np.cos(2 * j * theta[1:n]) / (4 * j * j - 1)
    v -= np.cos(n * theta[1:n]) / (n**2 - 1)
    w[1:n] = 2.0 * v / n
    return w


def diffusive_profile_gradient(z, t, ra, terms=40):
    """dC0/dz for C0(0, t) = 1, dC0/dz(1, t) = 0, C0(z, 0) = 0 as a series of images"""
    s = 2.0 * np.sqrt(t / ra)
    gradient = np.zeros_like(z)
    for m in range(terms):
        a, b = 2 * m + z, 2 * m + 2 - z
        gradient += (-1.0) ** m * (np.exp(-((b / s) ** 2)) - np.exp(-((a / s) ** 2)))
    return gradient * 2.0 / (np.sqrt(np.pi) * s)


def dominant_mode(ra, t0):
    """Initial perturbation zeta exp(-zeta^2), zeta = z sqrt(Ra / (4 t0)) (Riaz et al. 2006)"""
    def perturbation(z):
        zeta = z * np.sqrt(ra / (4.0 * t0))
        return zeta * np.exp(-zeta**2)
    return perturbation


class LinearizedProblem:
    def __init__(self, ra, k, n=192):
        self.ra, self.k = ra, k
        dx, x = chebyshev(n)
        self.z = 0.5 * (1.0 - x)
        dz = -2.0 * dx
        dzz = dz @ dz
        # unknowns at the interior points: c(0) = 0, c(1) from the Neumann condition
        self.extension = np.vstack([np.zeros(n - 1), np.eye(n - 1), -dz[n, 1:n] / dz[n, n]])
        if k > 0:
            self.velocity = np.linalg.solve(dzz[1:n, 1:n] - k**2 * np.eye(n - 1), -(k**2) * np.eye(n - 1))
        else:
            self.velocity = np.zeros((n - 1, n - 1))
        self.diffusion = (dzz[1:n, :] @ self.extension - k**2 * np.eye(n - 1)) / ra
        self.weights = 0.5 * clenshaw_curtis_weights(n)

    def operator(self, t):
        return self.diffusion - diffusive_profile_gradient(self.z[1:-1], t, self.ra)[:, None] * self.velocity

    def solve(self, perturbation, t0, times, rtol=1e-10):
        """Perturbation c(z, t) at the given times, one row per time"""
        c0 = perturbation(self.z)[1:-1]
        sol = solve_ivp(lambda t, c: self.operator(t) @ c, (t0, times[-1]), c0, method="Radau",
                        jac=lambda t, c: self.operator(t), t_eval=times, rtol=rtol,
                        atol=1e-14 * np.abs(c0).max())
        return (self.extension @ sol.y).T

    def l2_norm(self, c):
        return np.sqrt((c**2) @ self.weights)

    def max_norm(self, c):
        """Maximum of |c| between the collocation points"""
        result = []
        for ci in np.atleast_2d(c):
            i = np.argmax(np.abs(ci))
            fine = np.linspace(self.z[max(i - 1, 0)], self.z[min(i + 1, len(self.z) - 1)], 201)
            result.append(np.abs(BarycentricInterpolator(self.z, ci)(fine)).max())
        return np.array(result)


def growth_rate(times, norm):
    return np.gradient(np.log(norm), times)


def main():
    # k = 0: the perturbation z t^(-3/2) exp(-Ra z^2 / (4t)) of the semi-infinite layer is exact
    ra, t0, t1 = 500.0, 0.01, 0.5
    problem = LinearizedProblem(ra, 0.0, 128)
    c = problem.solve(dominant_mode(ra, t0), t0, np.array([t0, t1]))
    exact = problem.z * t1**-1.5 * np.exp(-ra * problem.z**2 / (4 * t1))
    exact *= c[0].max() / (problem.z * t0**-1.5 * np.exp(-ra * problem.z**2 / (4 * t0))).max()
    print(f"k = 0: max relative deviation from the exact solution at t = {t1}: "
          f"{np.abs(c[1] - exact).max() / c[1].max():.1e}")

    # growth rate of max|c| for the dominant-mode perturbation, Riaz et al. (2006), figure 4
    times = np.linspace(t0, 4.0, 800)
    for k in (4 * np.pi, 8 * np.pi):
        problem = LinearizedProblem(ra, k)
        sigma = growth_rate(times, problem.max_norm(problem.solve(dominant_mode(ra, t0), t0, times)))
        print(f"Ra = {ra:g}, k = {k / np.pi:g} pi: growth rate of max|c| at t = 1, 2, 4: "
              + ", ".join(f"{np.interp(t, times, sigma):.2f}" for t in (1, 2, 4)))

    # reference value of the test: ln of the amplification of the L2 norm, Ra = 500, k = 8 pi, t = 0.01 to 2
    for n in (128, 192, 256):
        problem = LinearizedProblem(ra, 8 * np.pi, n)
        c = problem.solve(dominant_mode(ra, t0), t0, np.array([t0, 2.0]))
        print(f"Ra = 500, k = 8 pi, {n} Chebyshev points: ln(|c|(2) / |c|(0.01)) = "
              f"{np.log(problem.l2_norm(c[1]) / problem.l2_norm(c[0])):.8f}")


if __name__ == "__main__":
    main()
