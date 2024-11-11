// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//

#ifndef DUMUX_TEST_MANDEL_ANALYTICAL_SOLUTION_HH
#define DUMUX_TEST_MANDEL_ANALYTICAL_SOLUTION_HH

#include <cmath>
#include <cassert>
#include <vector>
#include <iostream>
#include <algorithm>

#include <dumux/nonlinear/findscalarroot.hh>
#include "mandelelasticparams.hh"

namespace Dumux::Details {

/*!
 * \brief Solves for the nth x for function f(x) = tan(x)/x - a.
 * This function iteratively solves for the nth value using Brent's method.
 */
template <typename Scalar>
Scalar solveNthBeta(const Scalar a, const int n)
{
    // add some tolerance to the lower bound to avoid division by zero
    // the actual bounds are (n-1/2)π and (n+1/2)π
    const Scalar lowerBound = (1.0+1e-9)*(n - 1.0/2) * M_PI;
    const Scalar upperBound = lowerBound + (1.0-1e-5) * M_PI;

    using std::tan;

    return findScalarRootBrent(
        lowerBound, upperBound,
        [&](const Scalar x){ return tan(x)/x - a; },
        1e-8, 200
    );
}

/*!
 * /brief Solves for the beta values for the function f(x) = tan(x)/x - a.
 * This function computes the first n beta values
 */
template <typename Scalar>
std::vector<Scalar> solveBeta(const Scalar a, const int n)
{
    std::vector<Scalar> beta_(n);
    for (int i = 0; i < n; ++i)
        beta_[i] = solveNthBeta(a,i);
    return beta_;
}

} // end namespace Dumux::Details

namespace Dumux {

//! Analytical solution for Mandel's problem.
template<class Scalar>
class MandelAnalyticalSolution
{
public:
    MandelAnalyticalSolution()
    : params_{}
    {
        F_ = getParam<Scalar>("Problem.Force");
        a_ = getParam<std::vector<Scalar>>("Grid.UpperRight")[0];

        const Scalar mu = getParam<Scalar>("Problem.Viscosity");
        const Scalar K = getParam<Scalar>("Problem.Permeability");

        c_ = params_.M() * K / mu
             * (params_.Kb() + 4.0/3* params_.G())
             / (params_.Kbu() + 4.0/3 * params_.G());

        A1_ = 3/ params_.B() /( 1+ params_.Nuu());
        A2_ = params_.alpha() * (1 - 2* params_.Nu()) / (1- params_.Nu());

        const Scalar ratio = (1-params_.Nu())/(params_.Nuu()-params_.Nu());
        beta_ = Details::solveBeta(ratio, getParam<int>("Problem.NBeta"));
    }

    //! div(u) at t = 0
    template<class GlobalPosition>
    const Scalar initialDivU(const GlobalPosition&) const
    { return F_ * (2 * params_.Nuu() -1.0) / 2.0 / params_.G() /a_; }

    //! p at t = 0
    template<class GlobalPosition>
    const Scalar initialPressure(const GlobalPosition&) const
    { return F_ * params_.B() * (1 + params_.Nuu())/3.0/a_; }

    //! p at t = 0
    template<class GlobalPosition>
    const GlobalPosition initialDisplacement(const GlobalPosition& globalPos) const
    {
        static_assert(GlobalPosition::size() == 2, "GlobalPosition is expected to be 2D array");

        GlobalPosition u(0.0);
        const auto x = globalPos[0];
        const auto y = globalPos[1];

        u[0] = F_ * params_.Nuu() * x / 2.0 /params_.G() / a_;
        u[1] = - F_ * (1-params_.Nuu()) * y / 2.0 /params_.G() / a_;

        return u;
    }

    //! p(x,t)
    template<class GlobalPosition>
    Scalar pressure(const GlobalPosition& globalPos, const Scalar t) const
    {
        static_assert(GlobalPosition::size() == 2, "GlobalPosition is expected to be 2D array");

        Scalar p = 0.0;

        const auto x = globalPos[0];
        using std::sin, std::cos, std::exp;
        for (int i = 0; i < beta_.size(); i++)
        {
            p += sin(beta_[i])
                 /(beta_[i] - sin(beta_[i]) * cos(beta_[i]))
                 *(cos(beta_[i]*x/a_) - cos(beta_[i]))
                 * exp(- beta_[i] * beta_[i] * c_ * t / a_ /a_);
        }

        p *= 2.0 * F_* params_.B() * (1.0 + params_.Nuu()) / 3.0 / a_;

        return p;
    }

    //! u(x,t)
    template<class GlobalPosition>
    GlobalPosition displacement(const GlobalPosition& globalPos, const Scalar t) const
    {
        GlobalPosition u(0.0);
        const auto x = globalPos[0];
        const auto y = globalPos[1];

        using std::sin, std::cos, std::exp;
        for (int i = 0; i < beta_.size(); i++)
        {
            u[0] += sin(beta_[i]) * cos(beta_[i])
                    / (beta_[i] - sin(beta_[i]) * cos(beta_[i]))
                    * exp(- beta_[i] * beta_[i] * c_ * t / a_ /a_);


            u[1] += sin(beta_[i]) * cos(beta_[i])
                    / (beta_[i] - sin(beta_[i]) * cos(beta_[i]))
                    * exp(- beta_[i] * beta_[i] * c_ * t / a_ /a_);
        }

        u[0] *= -F_*params_.Nuu() / params_.G() / a_;
        u[1] *= F_*(1-params_.Nuu()) / params_.G() / a_;

        u[0] += F_ * params_.Nu() / 2.0 / params_.G() / a_;
        u[1] -= F_ * (1-params_.Nu()) / 2.0 / params_.G() / a_;

        u[0] *= x;
        u[1] *= y;

        Scalar part2 = 0.0;
        for (int i = 0; i < beta_.size(); i++)
        {
            part2 += cos(beta_[i])
                     / (beta_[i] - sin(beta_[i]) * cos(beta_[i]))
                     * sin(beta_[i] * x / a_)
                     * exp(- beta_[i] * beta_[i] * c_ * t / a_ / a_);
        }

        part2 *= F_/params_.G();

        u[0] += part2;
        return u;
    }

    //! σ_zz(x,t)
    template<class GlobalPosition>
    Scalar normalStressZZ(const GlobalPosition& globalPos, const Scalar t) const
    {
        Scalar sigmaZZ = - F_ /a_;
        Scalar part1 = 0.0;
        Scalar part2 = 0.0;

        const auto x = globalPos[0];
        for (int i = 0; i < beta_.size(); i++)
        {
            part1 += sin(beta_[i])
                     /(beta_[i] - sin(beta_[i]) *cos(beta_[i]))
                     * cos(beta_[i] * x / a_)
                     * exp(-beta_[i] * beta_[i] * c_ * t / a_ / a_);

            part2 += sin(beta_[i]) * cos(beta_[i])
                     /(beta_[i] - sin(beta_[i]) *cos(beta_[i]))
                     * exp(-beta_[i] * beta_[i] * c_ * t / a_ / a_);
        }

        part1 *=  - 2.0 * F_ * params_.B() * (params_.Nuu()-params_.Nu()) * (1 + params_.Nuu())
                  / a_ / (1- params_.Nu()) /  (1 + 2 * params_.Nuu());

        part2 *= 2.0 * F_ / a_;

        sigmaZZ += part1;
        sigmaZZ += part2;

        return sigmaZZ;
    }

    //! poroelatic parameters
    const auto& params() const { return params_; }

private:
    Scalar F_; // force intensity
    Scalar a_; // the length of the domain
    std::vector<Scalar> beta_; // beta values (roots of tan(x)/x - a = 0)

    Scalar c_; // compression coefficient
    Scalar A1_, A2_;

    MandelPoroElasticParameters<Scalar> params_; // poroelastic parameters
};

} // end namespace Dumux

#endif
