// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup PoromechanicsTests
 * \brief The analytical solution for Mandel's problem.
 * Reference:
 *  paper1: "Mandel’s problem revisited." DOI: 10.1680/geot.1996.46.2.187
 *  (Section 6.1) in paper2: "A coupling of mixed and continuous Galerkin finite element methods for poroelasticity I: the continuous in time case" DOI:10.1007/s10596-007-9045-y
 */


#ifndef DUMUX_TEST_MANDEL_ANALYTICAL_SOLUTION_HH
#define DUMUX_TEST_MANDEL_ANALYTICAL_SOLUTION_HH

#include <cmath>
#include <cassert>
#include <vector>
#include <iostream>
#include <algorithm>
#include <fstream>

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
             / (params_.Kbu() + 4.0/3 * params_.G()); // Eq.(21) in paper 1 with M (P-wave Modulus) = K + 4/3 G

        A1_ = 3/ params_.B() /( 1+ params_.Nuu()); // Eq.(27) in paper 1
        A2_ = params_.alpha() * (1 - 2* params_.Nu()) / (1- params_.Nu()); // Eq.(31) in the paper 1

        const Scalar ratio = (1-params_.Nu())/(params_.Nuu()-params_.Nu()); // Eq.(37) in paper 1 as A1/A2, see also the definition of  $\alpha_n$ in paper 2.

        beta_ = Details::solveBeta(ratio, getParam<int>("Problem.NBeta"));
    }

    //! div(u) at t = 0
    template<class GlobalPosition>
    const Scalar initialDivU(const GlobalPosition&) const
    { return F_ * (2 * params_.Nuu() -1.0) / 2.0 / params_.G() /a_; }

    //! div(u) at any time t
    template<class GlobalPosition>
    Scalar divU(const GlobalPosition& globalPos, const Scalar t) const
    {
        const auto x = globalPos[0];

        // From the displacement field u(x,y,t), compute div(u) = du_x/dx + du_y/dy
        //
        // u_x = x * [F_*Nu_/(2*G_*a_) + time_dependent_part] + spatial_part
        // u_y = y * [-F_*(1-Nu_)/(2*G_*a_) + time_dependent_part]
        //
        // du_x/dx = F_*Nu_/(2*G_*a_) + time_dependent_gradient + d(spatial_part)/dx
        // du_y/dy = -F_*(1-Nu_)/(2*G_*a_) + time_dependent_gradient

        using std::sin, std::cos, std::exp;

        // Constant part from initial displacement field
        Scalar dudx_const = F_ * params_.Nu() / 2.0 / params_.G() / a_;
        Scalar dudy_const = -F_ * (1-params_.Nu()) / 2.0 / params_.G() / a_;

        // Time-dependent part - derivative of the series terms
        Scalar time_dependent = 0.0;
        for (int i = 0; i < beta_.size(); i++)
        {
            time_dependent += sin(beta_[i]) * cos(beta_[i])
                            / (beta_[i] - sin(beta_[i]) * cos(beta_[i]))
                            * exp(- beta_[i] * beta_[i] * c_ * t / a_ /a_);
        }

        // For u_x: coefficient is -F_*Nuu_/G_/a_, for u_y: coefficient is F_*(1-Nuu_)/G_/a_
        Scalar dudx_time = -F_*params_.Nuu() / params_.G() / a_ * time_dependent;
        Scalar dudy_time = F_*(1-params_.Nuu()) / params_.G() / a_ * time_dependent;

        // Spatial derivative of the additional x-dependent term in u_x
        Scalar dudx_spatial = 0.0;
        for (int i = 0; i < beta_.size(); i++)
        {
            dudx_spatial += cos(beta_[i])
                          / (beta_[i] - sin(beta_[i]) * cos(beta_[i]))
                          * (beta_[i]/a_) * cos(beta_[i] * x / a_)
                          * exp(- beta_[i] * beta_[i] * c_ * t / a_ / a_);
        }
        dudx_spatial *= F_/params_.G();

        return dudx_const + dudy_const + dudx_time + dudy_time + dudx_spatial;
    }

    //! p at t = 0, see $p^{+}$ in section 6.1 in paper 2
    template<class GlobalPosition>
    const Scalar initialPressure(const GlobalPosition&) const
    { return F_ * params_.B() * (1 + params_.Nuu())/3.0/a_; }

    //! u at t = 0, see $S_{xo}, S_{yo}$ in section 6.1 in paper 2
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

    //! p(x,t), see p in paper 2
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

    //! u(x,t), see u_x, u_y in paper 2
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

    //! σ_yy(x,t) in paper 2
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

    //! Print all beta values to beta.txt file
    void printBetaValues() const
    {
        std::ofstream file("beta.txt");
        if (file.is_open())
        {
            file << "Beta values for Mandel's problem:\n";
            file << "Number of beta values: " << beta_.size() << "\n\n";
            for (size_t i = 0; i < beta_.size(); ++i)
            {
                file << "beta[" << i << "] = " << beta_[i] << "\n";
            }
            file.close();
            std::cout << "Beta values written to beta.txt" << std::endl;
        }
        else
        {
            std::cerr << "Unable to open beta.txt for writing" << std::endl;
        }
    }

private:
    Scalar F_; // force intensity
    Scalar a_; // the length of the domain
    std::vector<Scalar> beta_; // beta values (roots of tan(x)/x - a = 0)

    Scalar c_; // compression coefficient, see Eq.(21) in paper 1
    Scalar A1_, A2_; // some coefficients, see Eq.(27) and Eq.(31) in paper 1

    MandelPoroElasticParameters<Scalar> params_; // poroelastic parameters
};

} // end namespace Dumux

#endif
