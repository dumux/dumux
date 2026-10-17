#!/usr/bin/env python3
# SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

"""Onset of growth-induced wrinkling of a film on a substrate, plane strain, linearized about the
homogeneous base state.

Layers (bottom to top): substrate 0 < y < H (clamped at y = 0), film H < y < H + h (traction free
at the top). Both compressible neo-Hookean, P_e = mu F_e + (lam ln J_e - mu) F_e^-T, film growth
F_g = diag(theta, 1), stress P = J_g P_e(F F_g^-1) F_g^-T. Base state: u_x = 0, substrate unstrained,
film F = diag(1, ly) with P_yy = 0.

Increment u_x = U(y) sin(kx), u_y = V(y) cos(kx) with the incremental moduli A = dP/dF of each layer:
state Y = (U, V, T1, T2), T1 = A1212 U' - k A1221 V, T2 = k A2211 U + A2222 V', Y' = B Y per layer.
The critical growth theta_c(k) makes the 2x2 traction block of the propagator singular.
"""
import numpy as np
from scipy.linalg import expm
from scipy.optimize import brentq


def stress(F, Fg, mu, lam):
    Fe = F @ np.linalg.inv(Fg)
    FeInvT = np.linalg.inv(Fe).T
    Pe = mu*Fe + (lam*np.log(np.linalg.det(Fe)) - mu)*FeInvT
    return np.linalg.det(Fg)*Pe @ np.linalg.inv(Fg).T


def moduli(F, Fg, mu, lam, eps=1e-6):
    A = np.zeros((2, 2, 2, 2))
    for k in range(2):
        for l in range(2):
            dF = np.zeros((2, 2)); dF[k, l] = eps
            A[:, :, k, l] = (stress(F + dF, Fg, mu, lam) - stress(F - dF, Fg, mu, lam))/(2*eps)
    return A


def filmBaseStretch(theta, mu, lam):
    Fg = np.diag([theta, 1.0])
    return brentq(lambda ly: stress(np.diag([1.0, ly]), Fg, mu, lam)[1, 1], 0.5, 2.0, xtol=1e-15)


def layerMatrix(A, k):
    a1111, a2222 = A[0, 0, 0, 0], A[1, 1, 1, 1]
    a1212, a2121 = A[0, 1, 0, 1], A[1, 0, 1, 0]
    a1122, a2211 = A[0, 0, 1, 1], A[1, 1, 0, 0]
    a1221, a2112 = A[0, 1, 1, 0], A[1, 0, 0, 1]
    # U' = (T1 + k a1221 V)/a1212, V' = (T2 - k a2211 U)/a2222
    # T1' = k^2 a1111 U + k a1122 V', T2' = k^2 a2121 V - k a2112 U'
    dU = np.array([0.0, k*a1221/a1212, 1.0/a1212, 0.0])
    dV = np.array([-k*a2211/a2222, 0.0, 0.0, 1.0/a2222])
    B = np.zeros((4, 4))
    B[0] = dU
    B[1] = dV
    B[2] = np.array([k*k*a1111, 0, 0, 0]) + k*a1122*dV
    B[3] = np.array([0, k*k*a2121, 0, 0]) - k*a2112*dU
    return B


def determinant(theta, k, h, H, film, substrate):
    ly = filmBaseStretch(theta, *film)
    Af = moduli(np.diag([1.0, ly]), np.diag([theta, 1.0]), *film)
    As = moduli(np.eye(2), np.eye(2), *substrate)
    # the film thickness in the reference configuration is h
    P = expm(layerMatrix(Af, k)*h) @ expm(layerMatrix(As, k)*H)
    M = P[2:, 2:]
    return np.linalg.det(M)/np.linalg.norm(M)**2


def criticalGrowth(k, h, H, film, substrate, thetaMax=2.0, dTheta=1e-3):
    t0, d0 = 1.0 + 1e-9, determinant(1.0 + 1e-9, k, h, H, film, substrate)
    t = t0
    while t < thetaMax:
        t1 = t + dTheta
        d1 = determinant(t1, k, h, H, film, substrate)
        if np.sign(d1) != np.sign(d0):
            return brentq(determinant, t, t1, args=(k, h, H, film, substrate), xtol=1e-13)
        t, d0 = t1, d1
    return np.inf


def planeStrainModulus(mu, lam):
    return 4*mu*(lam + mu)/(lam + 2*mu)


def halfSpaceStiffness(k, mu, lam, depth):
    """normal traction per normal displacement amplitude at the top of a clamped layer with u_x = 0 there"""
    P = expm(layerMatrix(moduli(np.eye(2), np.eye(2), mu, lam), k)*depth)
    # top values (U, V) = P[:2, 2:] (a, b); require U = 0, V = 1
    ab = np.linalg.solve(P[:2, 2:], np.array([0.0, 1.0]))
    return (P[2:, 2:] @ ab)[1]


if __name__ == "__main__":
    # the propagator for a deep substrate against the bonded half-space with zero tangential surface
    # displacement, sigma_yy = 4 mu k (1 - nu)/(3 - 4 nu) w
    for mu, lam in [(1.0, 1.5), (1.0, 4.0)]:
        nu = lam/(2*(lam + mu))
        k = 0.7
        print(f"half space nu = {nu:.3f}: propagator {halfSpaceStiffness(k, mu, lam, 60.0):.12f},"
              f" closed form {4*mu*k*(1 - nu)/(3 - 4*nu):.12f}")

    # stiff film on a deep substrate against the small-strain asymptotics with the bonded substrate
    # stiffness, i.e. Es* = Es 4 (1 - nu)^2/(3 - 4 nu) in eps_c = (3 Es*/Ef)^(2/3)/4, k_c h = (3 Es*/Ef)^(1/3)
    for ratio in [1000.0, 3000.0, 10000.0]:
        film, substrate = (ratio, 1.5*ratio), (1.0, 1.5)
        nu = 0.3
        h, H = 1.0, 600.0
        Eratio = planeStrainModulus(*film)/planeStrainModulus(*substrate)/(4*(1 - nu)**2/(3 - 4*nu))
        kAsym = (3.0/Eratio)**(1/3)
        epsAsym = 0.25*(3.0/Eratio)**(2/3)
        from scipy.optimize import minimize_scalar
        res = minimize_scalar(lambda k: criticalGrowth(k, h, H, film, substrate, thetaMax=1.05, dTheta=2e-5),
                              bounds=(0.7*kAsym, 1.4*kAsym), method="bounded", options={"xatol": 1e-6})
        print(f"ratio {ratio:6.0f}: theta_c - 1 = {res.fun - 1:.5e} at k h = {res.x:.5f};"
              f" asymptotic eps_c = {epsAsym:.5e}, k_c h = {kAsym:.5f};"
              f" relative deviation {(res.fun - 1)/epsAsym - 1:+.4f}")

    # the benchmark: film mu = 10, substrate mu = 1, lambda = 1.5 mu, h = 1, H = 30
    film, substrate = (10.0, 15.0), (1.0, 1.5)
    h, H = 1.0, 30.0
    from scipy.optimize import minimize_scalar
    res = minimize_scalar(lambda k: criticalGrowth(k, h, H, film, substrate), bounds=(0.4, 1.2),
                          method="bounded", options={"xatol": 1e-8})
    print(f"benchmark continuum minimum: k_c h = {res.x:.6f}, wavelength {2*np.pi/res.x:.5f}, theta_c = {res.fun:.10f}")
    # sliding sides at x = 0 and x = L admit k = m pi / L
    L = 9.0
    for m in range(1, 5):
        k = m*np.pi/L
        print(f"benchmark L = {L}: m = {m}, k h = {k:.6f}, theta_c = {criticalGrowth(k, h, H, film, substrate):.10f}")
