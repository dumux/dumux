// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \brief Analytical solution for the Cryer spherical consolidation problem.
 *
 * The Cryer (1963) solution gives the pore pressure at position r and time t
 * in a sphere of radius R₀ subjected to a uniform surface pressure p₀:
 *
 * \f[
 *   \frac{p(r,t)}{p_0} = 2m \sum_{i=1}^{\infty}
 *       \frac{\sin z_i - z_i}
 *            {m z_i \cos z_i + (2m-1)\sin z_i}\,
 *       \frac{\sin(z_i r/R_0)}{z_i r/R_0}\,
 *       \exp\!\left(-\frac{z_i^2 c_v t}{\beta R_0^2}\right)
 * \f]
 *
 * where the \f$z_i\f$ are the positive real roots of
 * \f[ (1 - m z^2)\sin z - z \cos z = 0 \f]
 *
 * and
 * \f[
 *   m      = \frac{K + \tfrac{4}{3}G}{4G}
 *              \!\left(1 + \frac{K S_p}{\alpha_B^2}\right), \quad
 *   c_v    = \frac{\kappa}{\mu_f}(K + \tfrac{4}{3}G), \quad
 *   \beta  = \alpha_B^2 + (K + \tfrac{4}{3}G) S_p.
 * \f]
 *
 * References:
 * - Cryer, C.W. (1963) A comparison of the three-dimensional consolidation
 *   theories of Biot and Terzaghi.  *Q. J. Mech. Appl. Math.* 16(4):401–412.
 * - Verruijt, A. (2009) *An Introduction to Soil Dynamics.* Springer, ch. 3.
 */
#ifndef DUMUX_CRYER_ANALYTICAL_HH
#define DUMUX_CRYER_ANALYTICAL_HH

#include <cmath>
#include <vector>
#include <functional>
#include <stdexcept>

namespace Dumux {

/*!
 * \brief Analytical solution for Cryer's spherical consolidation problem.
 *
 * \tparam Scalar  Floating-point type (double or long double recommended).
 */
template<class Scalar = double>
class CryerAnalyticalSolution
{
public:
    /*!
     * \brief Constructor.
     *
     * \param K       Drained bulk modulus [Pa]
     * \param G       Shear modulus [Pa]
     * \param alphaB  Biot coefficient [-]
     * \param Sp      Storage coefficient [1/Pa]  (0 for incompressible constituents)
     * \param kappa   Intrinsic permeability [m²]
     * \param muf     Fluid dynamic viscosity [Pa·s]
     * \param R0      Sphere radius [m]
     * \param nRoots  Number of eigenvalues to sum (>=120 recommended for Figure-3 range)
     */
    CryerAnalyticalSolution(Scalar K, Scalar G, Scalar alphaB, Scalar Sp,
                             Scalar kappa, Scalar muf, Scalar R0,
                             int nRoots = 120)
    : R0_(R0)
    {
        // dimensionless parameters
        m_    = (K + 4.0/3.0*G) / (4.0*G) * (1.0 + K*Sp / (alphaB*alphaB));
        cv_   = kappa * (K + 4.0/3.0*G) / muf;
        beta_ = alphaB*alphaB + (K + 4.0/3.0*G)*Sp;

        roots_ = findRoots_(nRoots);
    }

    //! Consolidation coefficient c_v [m²/s]
    Scalar cv() const { return cv_; }
    //! Sphere radius R_0 [m]
    Scalar R0() const { return R0_; }

    /*!
     * \brief Normalised pore pressure p(r,t)/p_0 at radial distance r and time t.
     *
     * For r = 0 the limit sinc(0) = 1 is used.
     *
     * \param r  Radial distance from sphere centre [m],  0 <= r <= R0
     * \param t  Time [s],  t > 0
     */
    Scalar normalizedPressure(Scalar r, Scalar t) const
    {
        using std::sin; using std::cos; using std::exp;

        Scalar sum = 0.0;
        for (Scalar zi : roots_)
        {
            const Scalar numerator   = sin(zi) - zi;
            const Scalar denominator = m_*zi*cos(zi) + (2.0*m_ - 1.0)*sin(zi);
            if (std::abs(denominator) < 1e-14) continue;

            // sinc(z_i * r/R0) = sin(z_i*r/R0) / (z_i*r/R0)
            const Scalar xi = zi * r / R0_;
            const Scalar sinc = (xi < 1e-12) ? Scalar(1.0) : sin(xi) / xi;

            const Scalar exponent = -zi*zi * cv_ * t / (beta_ * R0_*R0_);
            if (exponent < -700.0) continue;   // underflow guard

            sum += numerator / denominator * sinc * exp(exponent);
        }
        return 2.0 * m_ * sum;
    }

    //! Convenience: pore pressure at the sphere centre (r = 0)
    Scalar normalizedCenterPressure(Scalar t) const
    { return normalizedPressure(Scalar(0), t); }

private:
    //! Characteristic equation: f(z) = (1 - m*z²)*sin(z) - z*cos(z)
    Scalar charEq_(Scalar z) const
    {
        using std::sin; using std::cos;
        return (1.0 - m_*z*z)*sin(z) - z*cos(z);
    }

    /*!
     * \brief Find the first `n` positive non-trivial roots of charEq_.
     *
     * We scan each interval (kπ + ε, (k+1)π − ε) using sub-intervals and
     * bisection on sign changes. This avoids missing roots in the first
     * interval and skips the trivial z≈0 root.
     */
    std::vector<Scalar> findRoots_(int n) const
    {
        using std::abs;
        constexpr Scalar pi = Scalar(3.14159265358979323846);
        constexpr Scalar eps = 1e-10;
        constexpr Scalar rootTol = 1e-6;
        constexpr int maxIter = 200;
        constexpr int nSubIntervals = 128;

        std::vector<Scalar> roots;
        roots.reserve(n);

        for (int k = 0; roots.size() < static_cast<std::size_t>(n); ++k)
        {
            const Scalar intervalLo = k * pi + eps;
            const Scalar intervalHi = (k + 1) * pi - eps;
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

    Scalar R0_, m_, cv_, beta_;
    std::vector<Scalar> roots_;
};

} // end namespace Dumux
#endif
