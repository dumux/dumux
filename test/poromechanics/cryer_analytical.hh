// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup PoromechanicsTests
 * \brief Analytical solution for the Cryer spherical consolidation problem.
 *
 * The Cryer (1963) solution gives the pore pressure at position r and time t
 * in a sphere of radius R₀ subjected to a uniform surface pressure p₀:
 *
 * \f[
 *   \frac{p(r,t)}{p_0} = 2 m B \sum_{i=1}^{\infty}
 *       \frac{\sin z_i - z_i \operatorname{sinc}(z_i r/R_0)}
 *            {m z_i \cos z_i + (2m-1)\sin z_i}\,
 *       \exp\!\left(-\frac{z_i^2 c_v t}{\beta R_0^2}\right),
 *   \qquad \operatorname{sinc}(x) = \frac{\sin x}{x},
 * \f]
 *
 * where the \f$z_i\f$ are the positive real roots of
 * \f[ (1 - m z^2)\sin z - z \cos z = 0 \f]
 *
 * and
 * \f[
 *   m = \frac{K + \tfrac{4}{3}G}{4G}
 *              \!\left(1 + \frac{K S_p}{\alpha_B^2}\right), \quad
 *   c_v = \frac{\kappa}{\mu_f}(K + \tfrac{4}{3}G), \quad
 *   \beta = \alpha_B^2 + (K + \tfrac{4}{3}G) S_p, \quad
 *   B = \frac{\alpha_B}{\alpha_B^2 + K S_p}.
 * \f]
 *
 * \f$B\f$ is Skempton's coefficient: the undrained response to the load is the
 * uniform pore pressure \f$B p_0\f$, and the series sums to \f$B\f$ at \f$t = 0^+\f$.
 * Each term vanishes at the drained surface \f$r = R_0\f$.
 *
 * Reference: Cryer, C.W. (1963) A comparison of the three-dimensional consolidation
 * theories of Biot and Terzaghi. *Q. J. Mech. Appl. Math.* 16(4):401–412.
 */
#ifndef DUMUX_CRYER_ANALYTICAL_HH
#define DUMUX_CRYER_ANALYTICAL_HH

#include <cmath>
#include <cstddef>
#include <vector>

namespace Dumux {

/*!
 * \ingroup PoromechanicsTests
 * \brief Analytical solution for Cryer's spherical consolidation problem.
 */
template<class Scalar = double>
class CryerAnalyticalSolution
{
public:
    /*!
     * \param K drained bulk modulus [Pa]
     * \param G shear modulus [Pa]
     * \param alphaB Biot coefficient [-]
     * \param Sp storage coefficient [1/Pa]
     * \param kappa intrinsic permeability [m²]
     * \param muf fluid dynamic viscosity [Pa·s]
     * \param R0 sphere radius [m]
     * \param nRoots number of terms of the series
     */
    CryerAnalyticalSolution(Scalar K, Scalar G, Scalar alphaB, Scalar Sp,
                             Scalar kappa, Scalar muf, Scalar R0,
                             int nRoots = 120)
    : R0_(R0)
    {
        m_ = (K + 4.0/3.0*G) / (4.0*G) * (1.0 + K*Sp / (alphaB*alphaB));
        cv_ = kappa * (K + 4.0/3.0*G) / muf;
        beta_ = alphaB*alphaB + (K + 4.0/3.0*G)*Sp;
        skemptonB_ = alphaB / (alphaB*alphaB + K*Sp);

        roots_ = findRoots_(nRoots);
    }

    //! Consolidation coefficient c_v [m²/s]
    Scalar cv() const { return cv_; }
    //! Sphere radius R_0 [m]
    Scalar R0() const { return R0_; }

    //! Normalised pore pressure p(r,t)/p_0 at the radial distance 0 <= r <= R_0 and the time t > 0
    Scalar normalizedPressure(Scalar r, Scalar t) const
    {
        using std::sin; using std::cos; using std::exp;

        Scalar sum = 0.0;
        for (Scalar zi : roots_)
        {
            const Scalar denominator = m_*zi*cos(zi) + (2.0*m_ - 1.0)*sin(zi);
            if (std::abs(denominator) < 1e-14) continue;

            const Scalar xi = zi * r / R0_;
            const Scalar sinc = (xi < 1e-12) ? Scalar(1.0) : sin(xi) / xi;
            const Scalar numerator = sin(zi) - zi*sinc;

            const Scalar exponent = -zi*zi * cv_ * t / (beta_ * R0_*R0_);
            if (exponent < -700.0) continue;

            sum += numerator / denominator * exp(exponent);
        }
        return 2.0 * m_ * skemptonB_ * sum;
    }

    //! Normalised pore pressure p(0,t)/p_0 at the sphere centre
    Scalar normalizedCenterPressure(Scalar t) const
    { return normalizedPressure(Scalar(0), t); }

private:
    Scalar charEq_(Scalar z) const
    {
        using std::sin; using std::cos;
        return (1.0 - m_*z*z)*sin(z) - z*cos(z);
    }

    //! The first n positive roots of the characteristic equation, by bisection on the sign changes over sub-intervals of (kπ, (k+1)π)
    std::vector<Scalar> findRoots_(int n) const
    {
        using std::abs;
        constexpr Scalar eps = 1e-10;
        constexpr Scalar rootTol = 1e-6;
        constexpr int maxIter = 200;
        constexpr int nSubIntervals = 128;

        std::vector<Scalar> roots;
        roots.reserve(n);

        for (int k = 0; roots.size() < static_cast<std::size_t>(n); ++k)
        {
            const Scalar intervalLo = k * M_PI + eps;
            const Scalar intervalHi = (k + 1) * M_PI - eps;
            const Scalar step = (intervalHi - intervalLo) / nSubIntervals;

            Scalar lo = intervalLo;
            Scalar flo = charEq_(lo);

            for (int s = 0; s < nSubIntervals
                     && roots.size() < static_cast<std::size_t>(n); ++s)
            {
                Scalar hi = (s == nSubIntervals - 1) ? intervalHi : (lo + step);
                Scalar fhi = charEq_(hi);

                if (flo * fhi > 0.0)
                {
                    lo = hi;
                    flo = fhi;
                    continue;
                }

                Scalar a = lo;
                Scalar b = hi;
                Scalar fa = flo;
                for (int it = 0; it < maxIter; ++it)
                {
                    const Scalar mid = 0.5*(a + b);
                    const Scalar fm = charEq_(mid);
                    if (abs(fm) < 1e-14 || (b - a) < 1e-14)
                    {
                        a = mid;
                        b = mid;
                        break;
                    }

                    if (fa * fm <= 0.0)
                        b = mid;
                    else
                    {
                        a = mid;
                        fa = fm;
                    }
                }

                const Scalar root = 0.5*(a + b);
                if (root > rootTol
                    && (roots.empty() || abs(root - roots.back()) > rootTol))
                    roots.push_back(root);

                lo = hi;
                flo = fhi;
            }
        }

        return roots;
    }

    Scalar R0_, m_, cv_, beta_, skemptonB_;
    std::vector<Scalar> roots_;
};

} // end namespace Dumux
#endif
