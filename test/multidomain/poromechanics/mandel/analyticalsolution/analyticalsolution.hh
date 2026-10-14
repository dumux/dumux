// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup PoromechanicsTests
 * \brief Analytical solution of Mandel's problem
 */
#ifndef DUMUX_TEST_MANDEL_ANALYTICAL_SOLUTION_HH
#define DUMUX_TEST_MANDEL_ANALYTICAL_SOLUTION_HH

#include <cmath>
#include <vector>

#include <dumux/common/parameters.hh>
#include <dumux/nonlinear/findscalarroot.hh>

#include "mandelparameters.hh"

namespace Dumux {

/*!
 * \ingroup PoromechanicsTests
 * \brief Analytical solution of Mandel's problem \cite Mandel1953
 *
 * A poroelastic slab \f$[-a, a] \times [-b, b]\f$ is loaded at \f$t = 0\f$ by two rigid plates
 * with a force \f$2F\f$ per unit length. The lateral sides \f$x = \pm a\f$ are traction-free and drained.
 * Pressure and displacement in the quarter domain \f$[0, a] \times [0, b]\f$ are series over the roots
 * \f$\beta_n\f$ of \f$\tan\beta/\beta = (1-\nu)/(\nu_u-\nu)\f$ \cite Abousleiman1996 \cite Phillips2007.
 */
template<class Scalar>
class MandelAnalyticalSolution
{
public:
    MandelAnalyticalSolution()
    : F_(getParam<Scalar>("Problem.Force"))
    , a_(getParam<std::vector<Scalar>>("Grid.UpperRight")[0])
    {
        const Scalar mobility = getParam<Scalar>("Problem.Permeability")/getParam<Scalar>("Component.LiquidDynamicViscosity");
        const Scalar K = params_.bulkModulus();
        const Scalar G = params_.shearModulus();
        consolidationCoefficient_ = mobility*params_.biotModulus()
                                    *(K + 4.0/3.0*G)/(params_.undrainedBulkModulus() + 4.0/3.0*G);

        const Scalar nu = params_.poissonRatio();
        const Scalar nuu = params_.undrainedPoissonRatio();
        const Scalar ratio = (1.0 - nu)/(nuu - nu);
        beta_.resize(getParam<int>("Problem.NumSeriesTerms"));
        for (std::size_t n = 0; n < beta_.size(); ++n)
            beta_[n] = root_(ratio, n);
    }

    //! The pressure \f$p^+ = FB(1+\nu_u)/(3a)\f$ of the undrained response at \f$t = 0\f$
    template<class GlobalPosition>
    Scalar initialPressure(const GlobalPosition&) const
    { return F_*params_.skemptonCoefficient()*(1.0 + params_.undrainedPoissonRatio())/(3.0*a_); }

    //! The displacement of the undrained response at \f$t = 0\f$
    template<class GlobalPosition>
    GlobalPosition initialDisplacement(const GlobalPosition& globalPos) const
    {
        const Scalar nuu = params_.undrainedPoissonRatio();
        const Scalar scale = F_/(2.0*params_.shearModulus()*a_);
        GlobalPosition u;
        u[0] = scale*nuu*globalPos[0];
        u[1] = -scale*(1.0 - nuu)*globalPos[1];
        return u;
    }

    //! The volumetric strain \f$\nabla \cdot \mathbf{u}\f$ of the undrained response at \f$t = 0\f$
    template<class GlobalPosition>
    Scalar initialVolumetricStrain(const GlobalPosition&) const
    { return F_*(2.0*params_.undrainedPoissonRatio() - 1.0)/(2.0*params_.shearModulus()*a_); }

    //! The pressure at time \f$t > 0\f$
    template<class GlobalPosition>
    Scalar pressure(const GlobalPosition& globalPos, const Scalar t) const
    {
        using std::sin; using std::cos;
        const Scalar x = globalPos[0];
        Scalar p = 0.0;
        for (std::size_t n = 0; n < beta_.size(); ++n)
        {
            const Scalar b = beta_[n];
            p += sin(b)/(b - sin(b)*cos(b))*(cos(b*x/a_) - cos(b))*decay_(b, t);
        }
        return 2.0*initialPressure(globalPos)*p;
    }

    //! The displacement at time \f$t > 0\f$
    template<class GlobalPosition>
    GlobalPosition displacement(const GlobalPosition& globalPos, const Scalar t) const
    {
        using std::sin; using std::cos;
        const Scalar x = globalPos[0];
        const Scalar G = params_.shearModulus();
        const Scalar nu = params_.poissonRatio();
        const Scalar nuu = params_.undrainedPoissonRatio();

        Scalar uniformStrainSum = 0.0;
        Scalar xSum = 0.0;
        for (std::size_t n = 0; n < beta_.size(); ++n)
        {
            const Scalar b = beta_[n];
            const Scalar denominator = b - sin(b)*cos(b);
            uniformStrainSum += sin(b)*cos(b)/denominator*decay_(b, t);
            xSum += cos(b)/denominator*sin(b*x/a_)*decay_(b, t);
        }

        GlobalPosition u;
        u[0] = (nu/2.0 - nuu*uniformStrainSum)*F_/(G*a_)*x + F_/G*xSum;
        u[1] = (-(1.0 - nu)/2.0 + (1.0 - nuu)*uniformStrainSum)*F_/(G*a_)*globalPos[1];
        return u;
    }

    //! The material parameters
    const MandelParameters<Scalar>& params() const
    { return params_; }

private:
    Scalar decay_(const Scalar beta, const Scalar t) const
    {
        using std::exp;
        return exp(-beta*beta*consolidationCoefficient_*t/(a_*a_));
    }

    //! The n-th positive root of \f$\tan\beta/\beta = r\f$ (\f$r > 1\f$) in \f$((n-1/2)\pi, (n+1/2)\pi)\f$
    static Scalar root_(const Scalar ratio, const std::size_t n)
    {
        using std::tan;
        const Scalar delta = 1e-8;
        const Scalar lower = n == 0 ? delta : (n - 0.5)*M_PI + delta;
        const Scalar upper = (n + 0.5)*M_PI - delta;
        return findScalarRootBrent(lower, upper, [&](const Scalar x){ return tan(x)/x - ratio; }, 1e-12);
    }

    MandelParameters<Scalar> params_;
    Scalar F_;
    Scalar a_;
    Scalar consolidationCoefficient_;
    std::vector<Scalar> beta_;
};

} // end namespace Dumux

#endif
