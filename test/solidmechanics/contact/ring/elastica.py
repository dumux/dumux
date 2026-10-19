#!/usr/bin/env python3
# SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

"""Reference for an inextensible elastic ring collapsing under uniform external pressure.

usage: elastica.py [h ...]
Units: midline radius R = 1, bending stiffness EI = 1, pressure p in EI/R^3, forces in EI/R^2.
The elliptic mode is symmetric about both axes. On the quarter from (x0, 0) to the minor axis,
with tangent angle theta and arc length s, force and moment balance give
    EI theta' = C + p/2 (x^2 + y^2) - Q/2 x,
where Q is the contact force between the two halves (zero before contact). Contact of a ring with
wall thickness t is modelled by the midline reaching y = h = t/2 on the minor axis.
  - free branch: no contact, parametrised by the semi-minor axis y1, from the bifurcation at p = 3;
  - point contact at (0, h) while the curvature there is negative (the walls curve away from each
    other), parametrised by Q;
  - line contact: a flat segment |x| <= b at y = h with zero curvature, the contact pressure p on
    it and point forces at its ends; the curved part s in [0, pi/2 - b] ends at (b, h) with
    theta = pi and zero curvature, parametrised by b.
"""
import sys

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq, fsolve


def shoot(x0, C, p, Q, length=0.5*np.pi):
    def kappa(x, y):
        return C + 0.5*p*(x*x + y*y) - 0.5*Q*x

    def rhs(s, z):
        x, y, theta, area = z
        return [np.cos(theta), np.sin(theta), kappa(x, y), 0.5*(x*np.sin(theta) - y*np.cos(theta))]
    sol = solve_ivp(rhs, [0.0, length], [x0, 0.0, 0.5*np.pi, 0.0], rtol=1e-12, atol=1e-13)
    x, y, theta, area = sol.y[:, -1]
    # the radial segments along the axes add nothing to the swept area; the flat part adds y*x/2
    return x, y, theta, 4.0*(area + 0.5*x*y), kappa(x, y)


def solve(residual, guess):
    v = fsolve(residual, guess, xtol=1e-11)
    assert np.linalg.norm(residual(v)) < 1e-9, (guess, v, residual(v))
    return v


def pointContact(Q, y1, guess):
    """Semi-major axis, C and pressure of the shape through (0, y1) with the force Q there."""
    def residual(v):
        x, y, theta, _, _ = shoot(v[0], v[1], v[2], Q)
        return [x, theta - np.pi, y - y1]
    return solve(residual, guess)


def lineContact(b, h, guess):
    """Semi-major axis, C, pressure and contact force of the shape flat on |x| <= b at y = h."""
    def residual(v):
        x, y, theta, _, kappa = shoot(v[0], v[1], v[2], v[3], 0.5*np.pi - b)
        return [x - b, y - h, theta - np.pi, kappa]
    return solve(residual, guess)


def reference(h, freeStep=2e-3, forceStep=0.1, widthStep=2e-3, maxPressure=14.0):
    """Branches as an array of rows (p, Q, area/pi, b) and the rows at first contact and at the onset of line contact."""
    rows = []
    v = pointContact(0.0, 1.0 - 1e-4, [1.0 + 1e-4, -0.5, 3.0])
    y1 = 1.0 - 1e-4
    while y1 > h:
        y1 = max(y1 - freeStep, h)
        v = pointContact(0.0, y1, v)
        rows.append((v[2], 0.0, shoot(*v, 0.0)[3]/np.pi, 0.0))
    first = rows[-1]

    Q, previous = 0.0, v
    while True:
        Q += forceStep
        v = pointContact(Q, h, previous)
        if shoot(*v, Q)[4] > 0.0:
            break
        rows.append((v[2], Q, shoot(*v, Q)[3]/np.pi, 0.0))
        previous = v
    QL = brentq(lambda q: shoot(*pointContact(q, h, previous), q)[4], Q - forceStep, Q, xtol=1e-12)
    v = pointContact(QL, h, previous)
    onset = (v[2], QL, shoot(*v, QL)[3]/np.pi, 0.0)
    rows.append(onset)

    w = np.array([v[0], v[1], v[2], QL])
    b = 0.0
    while w[2] < maxPressure:
        b += widthStep
        w = lineContact(b, h, w)
        rows.append((w[2], w[3], shoot(w[0], w[1], w[2], w[3], 0.5*np.pi - b)[3]/np.pi, b))
    return np.array(rows), first, onset


if __name__ == "__main__":
    heights = [float(a) for a in sys.argv[1:]] or [0.0, 0.005, 0.01, 0.02]
    p = pointContact(0.0, 1.0 - 1e-4, [1.0 + 1e-4, -0.5, 3.0])[2]
    print(f"bifurcation check: semi-minor axis 1-1e-4 -> p = {p:.6f} (linear theory 3)")
    for h in heights:
        rows, first, onset = reference(h)
        print(f"\nh = {h}: first contact p = {first[0]:.6f}, area/pi = {first[2]:.6f}")
        print(f"  line contact onset p = {onset[0]:.6f}, Q = {onset[1]:.6f}, area/pi = {onset[2]:.6f}")
        for target in [4, 5, 6, 7, 8, 9, 10, 11, 12, 13]:
            k = np.argmin(np.abs(rows[:, 0] - target))
            print(f"  p {rows[k, 0]:8.4f}  Q {rows[k, 1]:8.4f}  area/pi {rows[k, 2]:.5f}  b {rows[k, 3]:.4f}")
