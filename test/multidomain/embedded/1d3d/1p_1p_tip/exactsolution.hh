// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup EmbeddedTests
 * \brief Exact solution for a vessel with a closed tip embedded in a three-dimensional domain.
 */
#ifndef DUMUX_TEST_EMBEDDED_1D3D_VESSEL_TIP_EXACT_SOLUTION_HH
#define DUMUX_TEST_EMBEDDED_1D3D_VESSEL_TIP_EXACT_SOLUTION_HH

#include <cmath>
#include <limits>

#include <dune/common/exceptions.hh>

namespace Dumux {

/*!
 * \ingroup EmbeddedTests
 * \brief The distribution of the exchange term in the tissue: line source with circle average
 *        (line), source on the lateral vessel surface (surface), or uniform volume source (kernel)
 */
enum class VesselSourceDistribution { line, surface, kernel };

/*!
 * \ingroup EmbeddedTests
 * \brief Exact solution for a vessel with a closed tip embedded in a three-dimensional domain.
 *
 * The vessel center-line is the segment \f$ s \in [0, L] \f$ on the \f$ z \f$-axis. The solution is
 * mirrored at \f$ z = 0 \f$, i.e. it is the solution for the segment \f$ [-L, L] \f$ restricted to
 * \f$ z \geq 0 \f$, so that \f$ \partial p_t / \partial z = 0 \f$ on the plane \f$ z = 0 \f$.
 * The exchange rate per unit length is the constant \f$ q_0 \f$ and the vessel conductivity is the
 * constant \f$ \kappa \f$. With zero axial flux at the tip, the vessel pressure is
 * \f$ p_v(s) = p_\mathrm{tip} + \frac{q_0}{2\kappa} (L - s)^2 \f$.
 *
 * The tissue pressure solves \f$ -\Delta p_t = q_0 \delta \f$ in free space, where the
 * distribution \f$ \delta \f$ depends on the coupling method: a line source on the center-line,
 * a uniform surface source on the lateral cylinder surface of radius \f$ R \f$ (the end caps carry no
 * source), or a uniform volume source in the cylinder of radius \f$ \varrho \f$. All three are
 * expressed with the potential of a segment parallel to the \f$ z \f$-axis at distance \f$ d \f$,
 * \f$ F(d, z) = \frac{q_0}{4\pi} \left[ \operatorname{asinh}\frac{L - z}{d}
 *   + \operatorname{asinh}\frac{L + z}{d} \right] \f$,
 * averaged over a circle (surface source) or a disc (volume source). The radial integral over the disc
 * is evaluated in closed form, so that each method requires at most one numerical integration over
 * the angle.
 *
 * The exchange term of the coupled problem is \f$ q = \xi \beta(s) (p_v - \bar{p}_t) \f$, where
 * \f$ \bar{p}_t \f$ is the circle average of the tissue pressure at distance \f$ R \f$ (line and surface
 * source) or the tissue pressure on the center-line (volume source). The flux scaling factor of the volume
 * source method is \f$ \xi = \left( 1 + \frac{\beta}{2\pi} \left( \ln\frac{\varrho}{R} - \frac{1}{2} \right) \right)^{-1} \f$
 * (Koch et al. 2020, Eq. 2.12, unit tissue permeability and viscosity), and \f$ \xi = 1 \f$ otherwise.
 * The coefficient \f$ \beta(s) \f$ is chosen such that \f$ q = q_0 \f$, so that \f$ p_v \f$ and
 * \f$ q \f$ are the same for all source distributions.
 */
template<class Scalar>
class VesselTipExactSolution
{
public:
    /*!
     * \param radius the vessel radius \f$ R \f$
     * \param kernelRadius the radius \f$ \varrho \f$ of the volume source
     * \param length the vessel length \f$ L \f$
     * \param exchangeRate the exchange rate per unit length \f$ q_0 \f$
     * \param conductivity the vessel conductivity \f$ \kappa \f$
     * \param tipPressure the vessel pressure at the tip \f$ p_\mathrm{tip} \f$
     */
    VesselTipExactSolution(Scalar radius, Scalar kernelRadius, Scalar length,
                           Scalar exchangeRate, Scalar conductivity, Scalar tipPressure)
    : radius_(radius), kernelRadius_(kernelRadius), length_(length)
    , exchangeRate_(exchangeRate), conductivity_(conductivity), tipPressure_(tipPressure)
    {
        if (radius <= 0.0 || kernelRadius <= 0.0 || length <= 0.0 || conductivity <= 0.0)
            DUNE_THROW(Dune::InvalidStateException, "Radii, length and conductivity have to be positive");
    }

    Scalar length() const { return length_; }
    Scalar exchangeRate() const { return exchangeRate_; }

    //! The vessel pressure \f$ p_v(s) \f$
    Scalar vesselPressure(Scalar s) const
    { return tipPressure_ + 0.5*exchangeRate_/conductivity_*(length_ - s)*(length_ - s); }

    //! The tissue pressure for the line source at radial distance r and height z
    Scalar tissuePressureLineSource(Scalar r, Scalar z) const
    { return segmentPotential_(r, z); }

    //! The tissue pressure for the source on the lateral cylinder surface of radius \f$ R \f$
    Scalar tissuePressureSurfaceSource(Scalar r, Scalar z) const
    {
        using std::sin; using std::sqrt;
        const auto integrand = [&](Scalar phi)
        {
            const Scalar sinHalf = sin(0.5*phi);
            const Scalar d = sqrt((r - radius_)*(r - radius_) + 4.0*r*radius_*sinHalf*sinHalf);
            return segmentPotential_(d, z);
        };
        return integrateFromZeroToPi_(integrand)/M_PI;
    }

    //! The tissue pressure for the uniform volume source in the cylinder of radius \f$ \varrho \f$
    Scalar tissuePressureKernelSource(Scalar r, Scalar z) const
    {
        if (r == 0.0)
            return kernelCenterlinePressure(z);

        const auto integrand = [&](Scalar phi) { return discRayIntegral_(r, phi, z); };
        return integrateFromZeroToPi_(integrand)/M_PI;
    }

    //! The circle average of the line source solution at distance \f$ R \f$
    Scalar wallPressureLineSource(Scalar s) const
    { return tissuePressureLineSource(radius_, s); }

    //! The circle average of the surface source solution on the vessel surface
    Scalar wallPressureSurfaceSource(Scalar s) const
    { return tissuePressureSurfaceSource(radius_, s); }

    //! The volume source solution on the center-line
    Scalar kernelCenterlinePressure(Scalar s) const
    {
        using std::asinh; using std::sqrt; using std::abs;
        const Scalar rho = kernelRadius_;
        Scalar sum = 0.0;
        for (const Scalar a : {length_ - s, length_ + s})
            sum += rho*rho*asinh(a/rho) + a*sqrt(rho*rho + a*a) - a*abs(a);
        return exchangeRate_/(4.0*M_PI*rho*rho)*sum;
    }

    //! The tissue pressure for the given source distribution
    Scalar tissuePressure(VesselSourceDistribution distribution, Scalar r, Scalar z) const
    {
        switch (distribution)
        {
            case VesselSourceDistribution::line: return tissuePressureLineSource(r, z);
            case VesselSourceDistribution::surface: return tissuePressureSurfaceSource(r, z);
            case VesselSourceDistribution::kernel: return tissuePressureKernelSource(r, z);
        }
        DUNE_THROW(Dune::InvalidStateException, "Unknown source distribution");
    }

    //! The tissue pressure \f$ \bar{p}_t \f$ in the exchange term of the given source distribution
    Scalar couplingPressure(VesselSourceDistribution distribution, Scalar s) const
    {
        switch (distribution)
        {
            case VesselSourceDistribution::line: return wallPressureLineSource(s);
            case VesselSourceDistribution::surface: return wallPressureSurfaceSource(s);
            case VesselSourceDistribution::kernel: return kernelCenterlinePressure(s);
        }
        DUNE_THROW(Dune::InvalidStateException, "Unknown source distribution");
    }

    //! The coefficient \f$ \beta(s) \f$ for which the exchange rate is \f$ q_0 \f$
    Scalar couplingCoefficient(VesselSourceDistribution distribution, Scalar s) const
    {
        const Scalar rate = exchangeRate_/(vesselPressure(s) - couplingPressure(distribution, s));
        // for the volume source, beta*xi(beta) = rate is solved for beta
        const Scalar beta = distribution == VesselSourceDistribution::kernel
                            ? rate/(1.0 - rate*kernelShapeFactor_()) : rate;
        if (!(beta > 0.0))
            DUNE_THROW(Dune::InvalidStateException, "The parameters yield a non-positive coupling coefficient at s = " << s);
        return beta;
    }

private:
    Scalar kernelShapeFactor_() const
    {
        using std::log;
        return (log(kernelRadius_/radius_) - 0.5)/(2.0*M_PI);
    }

    /*!
     * \brief Potential of the segment \f$ [-L, L] \f$ at distance d from the axis and height z
     *
     * For points on the axis beyond the tip, the limit \f$ d \to 0 \f$ is finite and the
     * logarithms are combined so that d cancels.
     */
    Scalar segmentPotential_(Scalar d, Scalar z) const
    {
        using std::log; using std::sqrt; using std::asinh;
        const Scalar a = length_ - z;
        const Scalar b = length_ + z;
        const Scalar factor = exchangeRate_/(4.0*M_PI);
        if (a < 0.0 && b > 0.0)
            return factor*log((b + sqrt(b*b + d*d))/(-a + sqrt(a*a + d*d)));
        if (d == 0.0)
            return std::numeric_limits<Scalar>::infinity();
        return factor*(asinh(a/d) + asinh(b/d));
    }

    /*!
     * \brief Integral of the segment potential along a ray through the disc of radius \f$ \varrho \f$
     *
     * Evaluates \f$ \int_0^\varrho \frac{2 r'}{\varrho^2} F(d(r', \varphi), z) \,\mathrm{d}r' \f$ with
     * \f$ d^2 = (r' - c)^2 + e^2 \f$, \f$ c = r\cos\varphi \f$, \f$ e = r\sin\varphi \f$, using the
     * antiderivatives of \f$ (t + c) \operatorname{asinh}(a / \sqrt{t^2 + e^2}) \f$ in \f$ t = r' - c \f$.
     */
    Scalar discRayIntegral_(Scalar r, Scalar phi, Scalar z) const
    {
        using std::cos; using std::sin; using std::abs;
        const Scalar rho = kernelRadius_;
        const Scalar c = r*cos(phi);
        const Scalar e = abs(r*sin(phi));
        const Scalar sinHalf = sin(0.5*phi);
        // rho - r cos(phi) without cancellation for r close to rho and small angles
        const Scalar tUpper = (rho - r) + 2.0*r*sinHalf*sinHalf;
        const Scalar tLower = -c;

        Scalar sum = 0.0;
        for (const Scalar a : {length_ - z, length_ + z})
        {
            sum += antiderivativeT_(tUpper, e, a) + c*antiderivativeOne_(tUpper, e, a);
            sum -= antiderivativeT_(tLower, e, a) + c*antiderivativeOne_(tLower, e, a);
        }
        return exchangeRate_/(4.0*M_PI)*2.0/(rho*rho)*sum;
    }

    //! antiderivative of \f$ t \operatorname{asinh}(a/w) \f$, \f$ w = \sqrt{t^2 + e^2} \f$
    static Scalar antiderivativeT_(Scalar t, Scalar e, Scalar a)
    {
        using std::sqrt; using std::asinh; using std::abs;
        const Scalar w2 = t*t + e*e;
        if (w2 == 0.0)
            return 0.5*a*abs(a);
        return 0.5*(w2*asinh(a/sqrt(w2)) + a*sqrt(w2 + a*a));
    }

    //! antiderivative of \f$ \operatorname{asinh}(a/w) \f$, \f$ w = \sqrt{t^2 + e^2} \f$
    static Scalar antiderivativeOne_(Scalar t, Scalar e, Scalar a)
    {
        using std::sqrt; using std::asinh; using std::atan;
        const Scalar w2 = t*t + e*e;
        Scalar result = 0.0;
        if (w2 > 0.0)
            result += t*asinh(a/sqrt(w2));
        if (a != 0.0)
            result += a*asinh(t/sqrt(a*a + e*e));
        if (e > 0.0)
            result -= e*atan(t*a/(e*sqrt(w2 + a*a)));
        return result;
    }

    /*!
     * \brief Tanh-sinh quadrature on \f$ [0, \pi] \f$
     *
     * The integrands may have an integrable singularity at \f$ \varphi = 0 \f$ (evaluation points
     * on the source surface or inside the source volume). The double-exponential clustering of the
     * nodes at the end points resolves it, and the nodes close to \f$ \varphi = 0 \f$ are computed
     * without cancellation.
     */
    template<class Integrand>
    static Scalar integrateFromZeroToPi_(const Integrand& f)
    {
        using std::sinh; using std::cosh; using std::exp; using std::abs;
        constexpr Scalar tMax = 3.5;
        constexpr Scalar relTol = 1e-13;
        constexpr int maxLevel = 12;

        const auto contribution = [&](Scalar t)
        {
            const Scalar u = 0.5*M_PI*sinh(t);
            const Scalar expTwoU = exp(-2.0*abs(u));
            const Scalar onePlusTanh = u < 0.0 ? 2.0*expTwoU/(1.0 + expTwoU) : 2.0/(1.0 + expTwoU);
            const Scalar coshU = cosh(u);
            const Scalar weight = 0.5*M_PI*0.5*M_PI*cosh(t)/(coshU*coshU);
            return weight*f(0.5*M_PI*onePlusTanh);
        };

        Scalar h = 0.5;
        Scalar sum = contribution(0.0);
        for (int k = 1; k*h <= tMax; ++k)
            sum += contribution(k*h) + contribution(-k*h);
        Scalar integral = h*sum;

        for (int level = 1; level <= maxLevel; ++level)
        {
            h *= 0.5;
            for (int k = 1; k*h <= tMax; k += 2)
                sum += contribution(k*h) + contribution(-k*h);
            const Scalar refined = h*sum;
            if (abs(refined - integral) <= relTol*abs(refined))
                return refined;
            integral = refined;
        }

        DUNE_THROW(Dune::MathError, "Tanh-sinh quadrature did not converge");
    }

    Scalar radius_;
    Scalar kernelRadius_;
    Scalar length_;
    Scalar exchangeRate_;
    Scalar conductivity_;
    Scalar tipPressure_;
};

} // end namespace Dumux

#endif
