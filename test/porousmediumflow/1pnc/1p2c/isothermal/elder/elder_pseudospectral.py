#!/usr/bin/env python3
# SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

"""
Pseudo-spectral reference solution of the Elder problem with the Boussinesq approximation.

Python port of the pseudo-spectral method of
M. van Reeuwijk, S. A. Mathias, C. T. Simmons, J. D. Ward (2009), Insights from a pseudospectral approach
to the Elder problem, Water Resour. Res. 45, W04416, doi:10.1029/2008WR007421,
whose MATLAB implementation by M. van Reeuwijk and S. A. Mathias is archived at doi:10.5281/zenodo.3539026
under the Creative Commons Attribution 4.0 International license (CC-BY-4.0).

Dimensionless equations on the symmetric half domain [0, L] x [0, 1] with L = 2 (lengths in units of the
height H, time in units of phi H^2 / D_eff, z upwards):

    dc/dt + Ra (d(q_x c)/dx + d(q_z c)/dz) = lap c,    lap psi = dc/dx,    q_x = dpsi/dz,  q_z = -dpsi/dx,

with psi = 0 on all walls, dc/dx = 0 at x = 0 and x = L, c = 0 at z = 0, c = f(x) at z = 1 with the smoothed
source f(x) = (1 - tanh((x - 1)/delta))/2, and c = 0 at t = 0. The concentration is a cosine series and the
stream function a sine series in x with the wavenumbers k_m = m pi / L. In z, the transport equation is
collocated at the interior Chebyshev-Gauss-Lobatto points and the stream function of each mode is obtained
with the Chebyshev tau method. Time integration uses the implicit-explicit Runge-Kutta scheme ARS(4,4,3)
(Ascher, Ruuth, Spiteri 1997, doi:10.1016/S0168-9274(97)00056-1) with implicit diffusion and explicit advection.

Usage: python3 elder_pseudospectral.py <Ra> <delta> <nx> <nz> <dt> <output.npz> [times ...]
"""

import sys

import numpy as np
from scipy.fft import dct, dst, idct, idst

# ARS(4,4,3): implicit tableau (first column zero) and explicit tableau, both stiffly accurate
ARS443_IMPLICIT = np.array(
    [
        [0.0, 0.0, 0.0, 0.0, 0.0],
        [0.0, 1 / 2, 0.0, 0.0, 0.0],
        [0.0, 1 / 6, 1 / 2, 0.0, 0.0],
        [0.0, -1 / 2, 1 / 2, 1 / 2, 0.0],
        [0.0, 3 / 2, -3 / 2, 1 / 2, 1 / 2],
    ]
)
ARS443_EXPLICIT = np.array(
    [
        [0.0, 0.0, 0.0, 0.0, 0.0],
        [1 / 2, 0.0, 0.0, 0.0, 0.0],
        [11 / 18, 1 / 18, 0.0, 0.0, 0.0],
        [5 / 6, -5 / 6, 1 / 2, 0.0, 0.0],
        [1 / 4, 7 / 4, 3 / 4, -7 / 4, 0.0],
    ]
)
ARS443_GAMMA = 0.5


def smoothed_source(x, delta, source_length=1.0):
    """Source concentration at the top, a smoothed step from 1 to 0 at x = source_length"""
    if delta == 0.0:
        return np.where(x < source_length, 1.0, np.where(x > source_length, 0.0, 0.5))
    return 0.5 * (1.0 - np.tanh((x - source_length) / delta))


def chebyshev_differentiation_matrix(n):
    """Derivative matrix on the points xi_j = cos(j pi / n) (Trefethen, Spectral Methods in MATLAB)"""
    xi = np.cos(np.pi * np.arange(n + 1) / n)
    c = np.hstack([2.0, np.ones(n - 1), 2.0]) * (-1.0) ** np.arange(n + 1)
    dxi = xi[:, None] - xi[None, :]
    d = np.outer(c, 1.0 / c) / (dxi + np.eye(n + 1))
    return d - np.diag(d.sum(axis=1))


def chebyshev_coefficient_derivative(n):
    """Matrix mapping Chebyshev coefficients of p to those of dp/dxi"""
    d = np.zeros((n + 1, n + 1))
    for k in range(n + 1):
        for j in range(k + 1, n + 1, 2):
            d[k, j] = 2.0 * j
    d[0, :] *= 0.5
    return d


class ElderPseudoSpectral:
    """Pseudo-spectral discretization with nx collocation points in x and nz in z"""

    def __init__(self, rayleigh, delta, nx, nz, length=2.0, source_length=1.0):
        self.ra = rayleigh
        self.delta = delta
        self.length = length
        self.mx = nx - 1
        self.nz = nz - 1
        self.x = np.linspace(0.0, length, nx)
        self.k = np.arange(nx) * np.pi / length

        # Chebyshev-Gauss-Lobatto points, index 0 is the top z = 1
        n = self.nz
        xi = np.cos(np.pi * np.arange(n + 1) / n)
        self.z = 0.5 * (xi + 1.0)
        self.dz = 2.0 * chebyshev_differentiation_matrix(n)
        self.dzz = self.dz @ self.dz
        self.cheb = np.cos(np.outer(np.arange(n + 1), np.arange(n + 1)) * np.pi / n)
        self.cheb_inv = np.linalg.inv(self.cheb)

        # integral over z in [0, 1] of the interpolant, from int_{-1}^{1} T_j = 2/(1 - j^2) for even j
        j = np.arange(n + 1)
        moments = np.where(j % 2 == 0, 2.0 / (1.0 - j**2 + (j % 2)), 0.0)
        self.z_weights = 0.5 * moments @ self.cheb_inv

        # stream function and its z-derivative of each mode from the concentration profile of the mode (tau method)
        dc = 2.0 * chebyshev_coefficient_derivative(n)
        dcc = dc @ dc
        bc_rows = np.vstack([(-1.0) ** j, np.ones(n + 1)])  # T_j at z = 0 and z = 1
        keep = np.eye(n + 1)
        keep[n - 1 :, :] = 0.0
        self.psi_op = np.zeros((nx, n + 1, n + 1))
        self.dpsi_op = np.zeros((nx, n + 1, n + 1))
        for m in range(1, self.mx):
            op = dcc - self.k[m] ** 2 * np.eye(n + 1)
            op[n - 1 :, :] = bc_rows
            psi_coeffs = np.linalg.solve(op, -self.k[m] * keep @ self.cheb_inv)
            self.psi_op[m] = self.cheb @ psi_coeffs
            self.dpsi_op[m] = self.cheb @ dc @ psi_coeffs

        self.top = self.cos_coefficients(smoothed_source(self.x, delta, source_length))
        self._implicit_dt = None

    # transforms in x (axis 0): nodal values <-> series coefficients
    def cos_coefficients(self, v):
        a = dct(v, type=1, axis=0) / self.mx
        a[0] *= 0.5
        a[-1] *= 0.5
        return a

    def cos_values(self, a):
        y = a * self.mx
        y[0] *= 2.0
        y[-1] *= 2.0
        return idct(y, type=1, axis=0)

    def sin_coefficients(self, v):
        b = np.zeros_like(v)
        b[1:-1] = dst(v[1:-1], type=1, axis=0) / self.mx
        return b

    def sin_values(self, b):
        v = np.zeros_like(b)
        v[1:-1] = idst(b[1:-1] * self.mx, type=1, axis=0)
        return v

    def full_state(self, interior):
        a = np.zeros((self.mx + 1, self.nz + 1))
        a[:, 0] = self.top
        a[:, 1:-1] = interior
        return a

    def initial_state(self):
        return np.zeros((self.mx + 1, self.nz - 1))

    def stream_function(self, a):
        """psi_m(z_j) and dpsi_m/dz(z_j), the sine coefficients of psi and q_x"""
        return np.einsum("mij,mj->mi", self.psi_op, a), np.einsum("mij,mj->mi", self.dpsi_op, a)

    def advection(self, a):
        if self.ra == 0.0:
            return np.zeros((self.mx + 1, self.nz - 1))
        psi, dpsi = self.stream_function(a)
        qx = self.sin_values(dpsi)
        qz = self.cos_values(-self.k[:, None] * psi)
        c = self.cos_values(a)
        flux_x = self.sin_coefficients(qx * c)
        flux_z = self.cos_coefficients(qz * c)
        return (-self.ra * (self.k[:, None] * flux_x + flux_z @ self.dz.T))[:, 1:-1]

    def diffusion(self, a):
        return (a @ self.dzz.T - self.k[:, None] ** 2 * a)[:, 1:-1]

    def _prepare_implicit(self, dt):
        if self._implicit_dt == dt:
            return
        inner = self.dzz[1:-1, 1:-1]
        identity = np.eye(self.nz - 1)
        ops = identity[None] - ARS443_GAMMA * dt * (
            inner[None] - self.k[:, None, None] ** 2 * identity[None]
        )
        self._implicit_inv = np.linalg.inv(ops)
        self._implicit_dt = dt

    def step(self, u, dt):
        """One ARS(4,4,3) step of the interior values u"""
        self._prepare_implicit(dt)
        boundary = self.dzz[1:-1, 0][None, :] * self.top[:, None]
        stages, explicit, implicit = [u], [self.advection(self.full_state(u))], [None]
        for i in range(1, 5):
            rhs = u.copy()
            for j in range(i):
                rhs += dt * ARS443_EXPLICIT[i, j] * explicit[j]
                if j > 0:
                    rhs += dt * ARS443_IMPLICIT[i, j] * implicit[j]
            rhs += dt * ARS443_GAMMA * boundary
            y = np.einsum("mij,mj->mi", self._implicit_inv, rhs)
            stages.append(y)
            full = self.full_state(y)
            implicit.append(self.diffusion(full))
            if i < 4:
                explicit.append(self.advection(full))
        return stages[-1]

    def diagnostics(self, a):
        """Mean vertical gradient at the top (Sherwood number) and the bottom, and total solute mass"""
        mean = a[0]
        gradient = self.dz @ mean
        return {
            "sherwood": gradient[0],
            "bottom": gradient[-1],
            "mass": self.length * self.z_weights @ mean,
        }

    def run(self, times, dt, initial=None, dt_start=1e-7, growth=1.3, monitor=None):
        """Integrate to the given output times; the step grows geometrically from dt_start to dt"""
        u = self.initial_state() if initial is None else initial
        t, h = 0.0, min(dt_start, dt)
        snapshots, history = [], [(0.0, self.diagnostics(self.full_state(u)))]
        for t_out in times:
            while t < t_out - 1e-14:
                h_step = min(h, t_out - t)
                u = self.step(u, h_step)
                t += h_step
                h = min(h * growth, dt)
                if not np.all(np.isfinite(u)):
                    raise RuntimeError(f"non-finite solution at t = {t}")
                history.append((t, self.diagnostics(self.full_state(u))))
                if monitor:
                    monitor(t, u)
            snapshots.append(self.full_state(u).copy())
        return snapshots, history

    def evaluate(self, a, x, z):
        """Concentration of the series a on the tensor grid x (size p) times z (size q), shape (p, q)"""
        coeffs = a @ self.cheb_inv.T
        basis_x = np.cos(np.outer(x, self.k))
        theta = np.arccos(np.clip(2.0 * np.asarray(z) - 1.0, -1.0, 1.0))
        basis_z = np.cos(np.outer(theta, np.arange(self.nz + 1)))
        return basis_x @ coeffs @ basis_z.T

    def evaluate_velocity(self, a, x, z):
        """Dimensionless Darcy velocity (q_x, q_z) times Ra on the tensor grid x times z"""
        psi, dpsi = self.stream_function(a)
        theta = np.arccos(np.clip(2.0 * np.asarray(z) - 1.0, -1.0, 1.0))
        basis_z = np.cos(np.outer(theta, np.arange(self.nz + 1)))
        psi_z = psi @ self.cheb_inv.T @ basis_z.T
        dpsi_z = dpsi @ self.cheb_inv.T @ basis_z.T
        qx = np.sin(np.outer(x, self.k)) @ dpsi_z
        qz = -np.cos(np.outer(x, self.k)) @ (self.k[:, None] * psi_z)
        return self.ra * qx, self.ra * qz


def source_cosine_coefficients(
    delta, num_modes, length=2.0, source_length=1.0, num_quadrature=200001
):
    """Coefficients of the cosine series f(x) = sum_m f_m cos(m pi x / length) of the smoothed source"""
    x = np.linspace(0.0, length, num_quadrature)
    f = smoothed_source(x, delta, source_length)
    k = np.arange(num_modes) * np.pi / length
    if delta == 0.0:
        coeffs = np.empty(num_modes)
        coeffs[0] = source_length / length
        coeffs[1:] = 2.0 / length * np.sin(k[1:] * source_length) / k[1:]
        return coeffs
    weights = np.full(num_quadrature, x[1] - x[0])
    weights[[0, -1]] *= 0.5
    coeffs = 2.0 / length * (np.cos(np.outer(k, x)) * f) @ weights
    coeffs[0] *= 0.5
    return coeffs


def exact_solution_without_buoyancy(x, z, t, source_coeffs, length=2.0, num_terms=200):
    """Exact solution at Ra = 0 on the tensor grid x times z

    Each cosine mode of the source diffuses into the layer as c_m = f_m (sinh(k_m z)/sinh(k_m)
    - sum_n 2 (-1)^(n+1) n pi/(k_m^2 + n^2 pi^2) sin(n pi z) e^(-(k_m^2 + n^2 pi^2) t)).
    """
    x, z = np.asarray(x, dtype=float), np.asarray(z, dtype=float)
    k = np.arange(len(source_coeffs)) * np.pi / length
    n = np.arange(1, num_terms + 1)
    profiles = np.empty((len(k), len(z)))
    for m, km in enumerate(k):
        steady = (
            z.copy()
            if km == 0.0
            else np.exp(km * (z - 1.0)) * (1.0 - np.exp(-2.0 * km * z)) / (1.0 - np.exp(-2.0 * km))
        )
        rate = km**2 + (n * np.pi) ** 2
        amplitude = 2.0 * (-1.0) ** (n + 1) * n * np.pi / rate * np.exp(-rate * t)
        profiles[m] = steady - np.sin(np.outer(z, n * np.pi)) @ amplitude
    return np.cos(np.outer(x, k)) @ (source_coeffs[:, None] * profiles)


def main():
    if len(sys.argv) < 7:
        sys.exit(__doc__)
    rayleigh, delta = float(sys.argv[1]), float(sys.argv[2])
    nx, nz = int(sys.argv[3]), int(sys.argv[4])
    dt = float(sys.argv[5])
    output = sys.argv[6]
    times = [float(t) for t in sys.argv[7:]] or [0.03, 0.1, 0.3, 1.0]
    solver = ElderPseudoSpectral(rayleigh, delta, nx, nz)
    snapshots, history = solver.run(times, dt)
    keys = history[0][1].keys()
    np.savez(
        output,
        rayleigh=rayleigh,
        delta=delta,
        nx=nx,
        nz=nz,
        dt=dt,
        times=np.array(times),
        snapshots=np.array(snapshots),
        history_time=np.array([h[0] for h in history]),
        **{f"history_{key}": np.array([h[1][key] for h in history]) for key in keys},
    )
    for t, a in zip(times, snapshots):
        d = solver.diagnostics(a)
        print(
            f"t = {t:g}: Sh = {d['sherwood']:.12g}, bottom = {d['bottom']:.12g}, mass = {d['mass']:.12g}"
        )


if __name__ == "__main__":
    main()
