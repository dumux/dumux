#!/usr/bin/env python3
# SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: MIT
"""
Reference solution of the Henry problem with velocity-dependent dispersion of Fahs et al., Water Resour. Res. 52,
7382-7407 (2016), doi:10.1002/2016WR019288. A Galerkin method solves the steady stream function formulation of
Eqs. (10)-(13),

    a (Psi_XX + Psi_ZZ) - c_X = 0,
    div((bm I + D) grad c) - Q . grad c = 0,  Q = (1 + Psi_Z, -Psi_X),  D = bL ((1 - ra) Q Q^T/|Q| + ra |Q| I),

for Psi = psi/(q d) - Z and C = c - X/xi on 0 < X < xi, 0 < Z < 1, with Psi = C_Z = 0 at Z = 0, 1 and
Psi_X = C = 0 at X = 0, xi. Fahs et al. expand Psi and C in sine and cosine series, Eqs. (14) and (15). These
converge algebraically, because C_XX and Psi_ZZ do not vanish on the boundaries where C and Psi vanish. Here, Psi
and C are expanded in tensor products of Legendre polynomial combinations with the same boundary conditions and the
same truncation orders (Nm, Nn, Nr, Ns). For test case 2, smooth stretchings of X and Z refine the resolution at
the toe of the wedge, where the velocity, and with it the dispersion, vanishes on the bottom. The diffusion terms
are integrated by parts, the integrals are evaluated by Gauss-Legendre quadrature, and Newton's method solves the
nonlinear system. The truncation is raised in steps, each started from the previous solution, and the output
table states how much the isochlors changed in the last step.

usage: python3 henry_legendre.py <test case 1|2> <output.csv>
The Jacobian is dense. The finest truncation of test case 2 needs about 11 GB of memory.
"""
import sys

import numpy as np
import scipy.linalg
from numpy.polynomial.legendre import leggauss
from scipy.optimize import brentq

HEIGHTS = np.round(0.05*np.arange(1, 20), 2)  # heights Z of Tables D1 and D2 of Fahs et al. (2016)
TRUNCATIONS = {1: [(8, 20, 10, 20), (16, 40, 20, 40), (32, 80, 40, 80), (48, 120, 60, 120)],
               2: [(24, 60, 30, 60), (32, 80, 40, 80), (48, 120, 60, 120), (64, 160, 80, 160), (80, 200, 100, 200)]}
# the toe stagnation point of test case 2 lies at X = 1.66 on the bottom
DISCRETIZATION = {1: dict(), 2: dict(xStretching=(0.9, 1.66), zStretching=0.9, quadratureFactor=2.0)}


def parameters(case):
    """Nondimensional parameters from the dimensional ones of Table 2 of Fahs et al. (2016)"""
    q, k, g, densityDifference, viscosity, porosity, depth = 6.6e-5, 1.0204e-9, 9.81, 25.0, 1e-3, 0.35, 1.0
    Dm, aL, aT = {1: (18.86e-6, 0.0, 0.0), 2: (9.43e-8, 0.1, 0.01)}[case]
    return dict(a=q*viscosity/(k*g*densityDifference), bm=porosity*Dm/(q*depth), bL=aL/depth,
                ra=aT/aL if aL > 0 else 0.0)


def gauss(n, length):
    t, w = leggauss(n)
    return 0.5*length*(t + 1), 0.5*length*w


class LegendreBasis:
    """phi_k = (P_k - c_k P_{k+2})/sqrt(4k + 6) of t = 2s/L - 1, k = 0..N-1, on 0 < s < L, with c_k = 1, so that
    phi_k vanishes at both ends, or with c_k = k(k+1)/((k+2)(k+3)), so that its derivative vanishes there"""
    def __init__(self, N, L, zeroDerivative):
        k = np.arange(N)
        self.L = L
        self.c = k*(k + 1)/((k + 2)*(k + 3)) if zeroDerivative else np.ones(N)
        self.scale = 1/np.sqrt(4*k + 6)

    def __call__(self, s):
        """values and s-derivatives at the points s"""
        t = 2*np.asarray(s, dtype=float)/self.L - 1
        N = len(self.c)
        P, dP = np.zeros((len(t), N + 2)), np.zeros((len(t), N + 2))
        P[:, 0], P[:, 1], dP[:, 1] = 1.0, t, 1.0
        for k in range(1, N + 1):
            P[:, k + 1] = ((2*k + 1)*t*P[:, k] - k*P[:, k - 1])/(k + 1)
            dP[:, k + 1] = dP[:, k - 1] + (2*k + 1)*P[:, k]
        k = np.arange(N)
        return (P[:, k] - self.c*P[:, k + 2])*self.scale, (dP[:, k] - self.c*dP[:, k + 2])*self.scale*2/self.L


class Stretching:
    """Monotone map x(s) = s - A/k (sin(k (s - s0)) + sin(k s0)) of [0, L] onto itself, A < 1. It refines the
    resolution by 1/(1 - A) at x = center with k = 2 pi/L, x(s0) = center, or at x = 0 with k = pi/L, s0 = 0."""
    def __init__(self, L, A=0.0, center=None):
        self.L, self.A = L, A
        self.k = np.pi/L if center is None else 2*np.pi/L
        self.s0 = 0.0 if center is None else brentq(lambda s0: s0 - A/self.k*np.sin(self.k*s0) - center,
                                                    center - L/2, center + L/2, xtol=1e-15)

    def __call__(self, s):
        s = np.asarray(s, dtype=float)
        return s - self.A/self.k*(np.sin(self.k*(s - self.s0)) + np.sin(self.k*self.s0))

    def derivative(self, s):
        return 1 - self.A*np.cos(self.k*(np.asarray(s, dtype=float) - self.s0))

    def inverse(self, x):
        x = np.atleast_1d(np.asarray(x, dtype=float))
        grid = np.linspace(0.0, self.L, 4001)
        s = np.interp(x, self(grid), grid)
        for _ in range(8):
            s = s - (self(s) - x)/self.derivative(s)
        return s


class Space:
    """Basis functions b(s(x)) on [0, L], where s(x) inverts a stretching x(s)"""
    def __init__(self, basis, stretching):
        self.basis, self.stretching = basis, stretching

    def atReference(self, s):
        """values and x-derivatives at the points x(s)"""
        V, dV = self.basis(s)
        return V, dV/self.stretching.derivative(s)[:, None]

    def __call__(self, x):
        return self.atReference(self.stretching.inverse(x))


class HenryGalerkin:
    """Galerkin method for Psi = sum_mn A_mn f_m(Z) g_n(X) and C = sum_rs B_rs h_r(Z) k_s(X), with the trial
    functions as test functions. The spaces f, h share one stretching of Z, and g, k one of X."""
    def __init__(self, f, g, h, k, quadraturePoints, a, bm, bL=0.0, ra=0.0, xi=3.0):
        self.f, self.g, self.h, self.k = f, g, h, k
        self.a, self.bm, self.bL, self.ra, self.xi = a, bm, bL, ra, xi
        sz, wz = gauss(quadraturePoints[0], 1.0)
        sx, wx = gauss(quadraturePoints[1], xi)
        self.wz, self.wx = wz*f.stretching.derivative(sz), wx*g.stretching.derivative(sx)
        self.W = self.wz[:, None]*self.wx[None, :]
        F, Fd = f.atReference(sz)
        G, Gd = g.atReference(sx)
        H, Hd = h.atReference(sz)
        K, Kd = k.atReference(sx)
        # (Z factor, X factor) of the stream function (phi) and concentration (chi) basis functions and their
        # X and Z derivatives at the quadrature points
        self.phi, self.phiX, self.phiZ = (F, G), (F, Gd), (Fd, G)
        self.chi, self.chiX, self.chiZ = (H, K), (H, Kd), (Hd, K)
        self.shapeA, self.shapeB = (F.shape[1], G.shape[1]), (H.shape[1], K.shape[1])
        self.nA = F.shape[1]*G.shape[1]
        self.size = self.nA + H.shape[1]*K.shape[1]

    def split(self, u):
        return u[:self.nA].reshape(self.shapeA), u[self.nA:].reshape(self.shapeB)

    def embed(self, coarse, uCoarse):
        """Coefficients of a solution with smaller truncation orders in the same bases"""
        u = np.zeros(self.size)
        for fine, old in zip(self.split(u), coarse.split(uCoarse)):
            m, n = min(fine.shape[0], old.shape[0]), min(fine.shape[1], old.shape[1])
            fine[:m, :n] = old[:m, :n]
        return u

    def fields(self, u):
        """Q and grad c at the quadrature points"""
        A, B = self.split(u)
        Q = np.array([1 + self.phiZ[0] @ A @ self.phiZ[1].T, -(self.phiX[0] @ A @ self.phiX[1].T)])
        gradc = np.array([self.chiX[0] @ B @ self.chiX[1].T + 1/self.xi, self.chiZ[0] @ B @ self.chiZ[1].T])
        return Q, gradc

    def flux(self, Q, gradc):
        """F = (bm I + D(Q)) grad c and its derivatives dF_i/d(grad c)_j and dF_i/dQ_j"""
        q = np.maximum(np.hypot(*Q), 1e-300)
        I = np.eye(2)[:, :, None, None]
        QQ = Q[:, None]*Q[None, :]
        Qg = np.sum(Q*gradc, axis=0)
        dFdgradc = self.bm*I + self.bL*((1 - self.ra)*QQ/q + self.ra*q*I)
        dFdQ = self.bL*((1 - self.ra)*(Q[:, None]*gradc[None, :] + Qg*(I - QQ/q**2))
                        + self.ra*gradc[:, None]*Q[None, :])/q
        return np.einsum("ij...,j...->i...", dFdgradc, gradc), dFdgradc, dFdQ

    def weak(self, test, field):
        return test[0].T @ (self.W*field) @ test[1]

    def residual(self, u):
        Q, gradc = self.fields(u)
        F = self.flux(Q, gradc)[0]
        flow = (-self.a*(self.weak(self.phiX, -Q[1]) + self.weak(self.phiZ, Q[0] - 1))
                - self.weak(self.phi, gradc[0]))
        transport = -(self.weak(self.chiX, F[0]) + self.weak(self.chiZ, F[1])
                      + self.weak(self.chi, np.sum(Q*gradc, axis=0)))
        return np.concatenate([flow.ravel(), transport.ravel()])

    def addForm(self, target, field, test, trial):
        """target[g, h, p, q] += sum over the quadrature points of w field TZ_g TX_h PZ_p PX_q"""
        (TZ, TX), (PZ, PX) = test, trial
        Ng, Np = TZ.shape[1], PZ.shape[1]
        G = (self.wz[:, None, None]*TZ[:, :, None]*PZ[:, None, :]).reshape(len(self.wz), Ng*Np)
        T = (G.T @ field)*self.wx[None, :]
        for i in range(Ng*Np):
            g, p = divmod(i, Np)
            target[g, :, p, :] += (TX.T*T[i]) @ PX

    def jacobian(self, u):
        Q, gradc = self.fields(u)
        _, dFdgradc, dFdQ = self.flux(Q, gradc)
        J = np.zeros((self.size, self.size))
        nA, (nm, nn), (nr, ns) = self.nA, self.shapeA, self.shapeB
        # (test Z, test X, trial Z, trial X) views of the four blocks
        flowA, flowB = J[:nA, :nA].reshape(nm, nn, nm, nn), J[:nA, nA:].reshape(nm, nn, nr, ns)
        transportA, transportB = J[nA:, :nA].reshape(nr, ns, nm, nn), J[nA:, nA:].reshape(nr, ns, nr, ns)
        one = np.ones_like(self.W)
        gradPhi, gradChi = (self.phiX, self.phiZ), (self.chiX, self.chiZ)
        dQ = ((1.0, self.phiZ), (-1.0, self.phiX))  # dQ_x = d Psi_Z, dQ_z = -d Psi_X
        for gradient in gradPhi:
            self.addForm(flowA, -self.a*one, gradient, gradient)
        self.addForm(flowB, -one, self.phi, self.chiX)
        for j, (sign, dQj) in enumerate(dQ):
            self.addForm(transportA, -sign*gradc[j], self.chi, dQj)
            self.addForm(transportB, -Q[j], self.chi, gradChi[j])
            for i in range(2):
                if i == j or self.bL > 0:
                    self.addForm(transportB, -dFdgradc[i, j], gradChi[i], gradChi[j])
                if self.bL > 0:
                    self.addForm(transportA, -sign*dFdQ[i, j], gradChi[i], dQj)
        return J

    def solve(self, u0=None, source=0.0, tol=1e-13, maxIterations=60):
        """Newton's method with backtracking line search for residual(u) = source"""
        u = np.zeros(self.size) if u0 is None else u0.copy()
        R = self.residual(u) - source
        norm = np.linalg.norm(R)
        for iteration in range(1, maxIterations + 1):
            # J is C-ordered, so J.T is Fortran-ordered and is factorized in place without a copy
            lu = scipy.linalg.lu_factor(self.jacobian(u).T, overwrite_a=True, check_finite=False)
            du = scipy.linalg.lu_solve(lu, R, trans=1, check_finite=False)
            del lu  # frees the factorization before the next Jacobian is assembled
            damping = 1.0
            while True:
                uNew = u - damping*du
                RNew = self.residual(uNew) - source
                normNew = np.linalg.norm(RNew)
                if normNew < (1 - 1e-4*damping)*norm or damping < 1/256:
                    break
                damping /= 2
            if not np.isfinite(normNew):
                raise RuntimeError("Newton's method produced a non-finite residual")
            u, R, norm = uNew, RNew, normNew
            if damping*np.max(np.abs(du)) < tol*max(1.0, np.max(np.abs(u))):
                return u, iteration
        raise RuntimeError(f"Newton's method did not converge, residual norm {norm:.3e}")

    def concentration(self, u, X, Z):
        """c and c_X on the grid Z x X (rows Z)"""
        X = np.atleast_1d(np.asarray(X, dtype=float))
        B = self.split(u)[1]
        H = self.h(np.atleast_1d(Z))[0]
        K, Kd = self.k(X)
        return H @ B @ K.T + X[None, :]/self.xi, H @ B @ Kd.T + 1/self.xi

    def velocityX(self, u, X, Z):
        A = self.split(u)[0]
        return 1 + self.f(np.atleast_1d(Z))[1] @ A @ self.g(np.atleast_1d(X))[0].T

    def isochlor(self, u, Z, level, samples=6001):
        """X of the first crossing of c = level from the inland boundary at height Z"""
        X = np.linspace(0.0, self.xi, samples)
        c = self.concentration(u, X, Z)[0][0]
        i = int(np.argmax(c >= level))
        if c[i] < level:
            raise ValueError(f"no isochlor {level} at Z = {Z}")
        return brentq(lambda x: self.concentration(u, x, Z)[0][0, 0] - level, X[i - 1], X[i], xtol=1e-15, rtol=1e-15)

    def isochlorTable(self, u, heights, levels=(0.1, 0.5, 0.9)):
        return np.array([[self.isochlor(u, Z, level) for level in levels] for Z in heights])

    def metrics(self, u):
        """L_toe, L_S, Z_I and Q_S as defined in Section 5 of Fahs et al. (2016), with the advective salt flux
        integrated over the inflow zone Z < Z_I and the diffusive and dispersive ones over the whole sea boundary,
        Eqs. (22)-(24), and the defect of the steady salt balance"""
        X10, X50, X90 = (self.isochlor(u, 0.0, level) for level in (0.1, 0.5, 0.9))
        Zs = np.linspace(0.0, 1.0, 4001)
        i = int(np.argmax(self.velocityX(u, self.xi, Zs)[:, 0] > 0))
        ZI = brentq(lambda z: self.velocityX(u, self.xi, z)[0, 0], Zs[i - 1], Zs[i], xtol=1e-15)
        z, w = gauss(400, ZI)
        advective = -w @ self.velocityX(u, self.xi, z)[:, 0]
        z, w = gauss(800, 1.0)

        def normalFlux(X):
            # Q_z = 0 on the vertical boundaries, so D_xx = bL |Q_x| and D_xz = 0
            Qx = self.velocityX(u, X, z)[:, 0]
            return w @ ((self.bm + self.bL*np.abs(Qx))*self.concentration(u, X, z)[1][:, 0])

        # the net inflow through the sea boundary, with int Q_x dZ = 1, leaves through the inland boundary c = 0
        sea, inland = normalFlux(self.xi), normalFlux(0.0)
        return dict(LToe=self.xi - X50, LS=X90 - X10, ZI=ZI, QS=advective + sea, saltBalance=sea - 1.0 - inland)


def legendreProblem(orders, params, xStretching=(0.0, None), zStretching=0.0, quadratureFactor=1.0, xi=3.0):
    """Psi with f_m = 0 at Z = 0, 1, m = 1..Nm, and g_n' = 0 at X = 0, xi, n = 0..Nn; C with h_r' = 0 at
    Z = 0, 1, r = 0..Nr, and k_s = 0 at X = 0, xi, s = 1..Ns"""
    Nm, Nn, Nr, Ns = orders
    x, z = Stretching(xi, *xStretching), Stretching(1.0, zStretching)
    spaces = (Space(LegendreBasis(Nm, 1.0, False), z), Space(LegendreBasis(Nn + 1, xi, True), x),
              Space(LegendreBasis(Nr + 1, 1.0, True), z), Space(LegendreBasis(Ns, xi, False), x))
    # exact for a stream function times two concentration basis functions; stretchings and the dispersion
    # tensor make the integrands non-polynomial and need the factor
    points = [int(quadratureFactor*((N1 + 2*N2 + 8)/2 + 10)) for N1, N2 in ((Nm, Nr), (Nn, Ns))]
    return HenryGalerkin(*spaces, points, xi=xi, **params)


def solveByContinuation(problem):
    """Newton continuation from bm = 0.1, bL = 0 to the parameters of the problem, with bm interpolated
    geometrically and bL linearly"""
    bm, bL = problem.bm, problem.bL
    problem.bm, problem.bL = 0.1, 0.0
    u, _ = problem.solve()
    progress, step = 0.0, 0.25
    while progress < 1.0:
        trial = min(1.0, progress + step)
        problem.bm, problem.bL = 0.1**(1 - trial)*bm**trial, trial*bL
        try:
            uTrial, iterations = problem.solve(u, maxIterations=25)
        except (RuntimeError, np.linalg.LinAlgError):
            step /= 2
            if step < 1e-3:
                raise
            continue
        progress, u = trial, uTrial
        if iterations <= 6:
            step *= 1.5
    return u


def main(case, output):
    params = parameters(case)
    previous, tables = None, []
    for orders in TRUNCATIONS[case]:
        problem = legendreProblem(orders, params, **DISCRETIZATION[case])
        if previous is not None:
            u, _ = problem.solve(problem.embed(*previous))
        elif params["bL"] > 0:
            u = solveByContinuation(problem)
        else:
            u, _ = problem.solve()
        tables.append(problem.isochlorTable(u, HEIGHTS))
        metrics = problem.metrics(u)
        change = np.max(np.abs(tables[-1] - tables[-2])) if len(tables) > 1 else np.nan
        print(f"truncation {orders}: {problem.size} unknowns, maximum isochlor change {change:.1e}, "
              f"L_toe {metrics['LToe']:.6f}, L_S {metrics['LS']:.6f}, Z_I {metrics['ZI']:.6f}, "
              f"Q_S {metrics['QS']:.6f}, salt balance defect {metrics['saltBalance']:.1e}", flush=True)
        previous = (problem, u)

    with open(output, "w", encoding="utf-8") as f:
        # REUSE-IgnoreStart
        f.write("# SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder\n"
                "# SPDX-License-Identifier: CC-BY-4.0\n"
                f"# Henry problem, test case {case} of Fahs et al. (2016), doi:10.1002/2016WR019288: height Z and\n"
                "# positions X of the 10, 50 and 90 % isochlors, both divided by the aquifer depth, computed by\n"
                f"# henry_legendre.py with the truncation orders {TRUNCATIONS[case][-1]}. They differ by at most\n"
                f"# {change:.1e} from those with the truncation orders {TRUNCATIONS[case][-2]}.\n"
                f"# L_toe = {metrics['LToe']:.6f}, L_S = {metrics['LS']:.6f}, Z_I = {metrics['ZI']:.6f}, "
                f"Q_S = {metrics['QS']:.6f}\n"
                "Z,X10,X50,X90\n")
        # REUSE-IgnoreEnd
        for Z, X in zip(HEIGHTS, tables[-1]):
            f.write(f"{Z:.3f},{X[0]:.6f},{X[1]:.6f},{X[2]:.6f}\n")


if __name__ == "__main__":
    if len(sys.argv) != 3 or sys.argv[1] not in ("1", "2"):
        sys.exit(__doc__)
    main(int(sys.argv[1]), sys.argv[2])
