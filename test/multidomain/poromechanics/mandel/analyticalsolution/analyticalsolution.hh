// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup PoromechanicsTests
 * \brief Solver for the equation tan(x)/x = a
 */


#ifndef TEST_MANDEL_ANALYTICAL_SOLUTION_HH
#define TEST_MANDEL_ANALYTICAL_SOLUTION_HH

#include <cmath>
#include <cassert>
#include <vector>
#include <iostream>
#include <algorithm>

#include "mandelelasticparams.hh"
namespace Dumux{
namespace Details{

using namespace std;


/**
 * @brief Computes the value of the function f(x) = tan(x) / x - a.
 */
template<class Scalar>
Scalar fx(const Scalar a, const Scalar x)
{
    return tan(x) / x - a;
}


/**
 * @brief Computes the derivative of the function f(x) = tan(x) / x - a.
 */
template<class Scalar>
Scalar fxx(const Scalar a, const Scalar x)
{
    return 1.0 / cos(x) / cos(x) / x - tan(x) / x / x;
}

/**
 * @brief Solves for the nth x for function f(x) = tan(x) / x - a.
 *
 * This function iteratively solves for the nth value using a combination
 * of Newton's method and a clamping strategy to ensure convergence within
 * specified bounds [(n-1/2)*Pi, (n+1/2)*Pi].
 *
 */
template <typename Scalar>
Scalar solveNthBeta(const Scalar a, const int n)
{
    // Initial guess
    Scalar xOld = (n+0.1) * M_PI;
    Scalar xNew = xOld;

    Scalar lowerBound = (n - 1.0/2) * M_PI;
    Scalar upperBound = lowerBound + M_PI;

    do
    {
        xOld = xNew;
        Scalar factor = 1.0;
        do
        {
            xNew = - fx(a,xOld) / fxx(a,xOld)/factor + xOld;
            factor *= 2;
        } while (xNew != std::clamp(xNew,lowerBound,upperBound));

    } while (abs(fx(a,xNew)) > 1e-14 && abs(xNew-xOld)/xOld > 1e-30);

    return xNew;
}

/**
 * @brief Solves for the beta values for the function f(x) = tan(x) / x - a.
 *
 * This function computes the first n beta values for the function f(x) = tan(x) / x - a.
 */
template <typename Scalar>
std::vector<Scalar> solveBeta(const Scalar a, const int n)
{
    std::vector<Scalar> beta_(n);
    for (int i = 0; i < n; i++)
    {
        beta_[i] = solveNthBeta(a,i);
    }
    return beta_;
}
}//end namespace Details

/*!
 * \brief Analytical solution for Mandel's problem.
 *
 * \tparam Scalar
 */
template<class Scalar>
class MandelAnalyticalSolution
{
public:
    MandelAnalyticalSolution()
    :param_()
    {
        F_ = getParam<Scalar>("Problem.Force");
        a_ = getParam<std::vector<Scalar>>("Grid.UpperRight")[0];
        nBeta_ = getParam<int>("Analytical.NBeta");
        beta_.resize(nBeta_);

        const Scalar mu = getParam<Scalar>("MaterialParameters.Viscosity");
        const Scalar K = getParam<Scalar>("SpatialParams.Permeability");

        c_ = param_.M() * K / mu
             * (param_.Kb() + 4.0/3* param_.G())
             / (param_.Kbu() + 4.0/3 * param_.G());

        A1_ = 3/ param_.B() /( 1+ param_.Nuu());
        A2_ = param_.alpha() * (1 - 2* param_.Nu()) / (1- param_.Nu());

        const Scalar ratio = (1-param_.Nu())/(param_.Nuu()-param_.Nu());
        beta_ = Details::solveBeta(ratio,nBeta_);
    }

    template<class GlobalPosition>
    Scalar pressureAtPos(const GlobalPosition& globalPos,
                         const Scalar time) const
    {
        Scalar p = 0.0;
        using namespace std;
        for (int i = 0; i < nBeta_; i++)
        {
            p += sin(beta_[i])
                 /(beta_[i] - sin(beta_[i]) * cos(beta_[i]))
                 *(cos(beta_[i]*globalPos[0]/a_) - cos(beta_[i]))
                 * exp(- beta_[i] * beta_[i] * c_ * time / a_ /a_);
        }

        p *= 2.0 * F_* param_.B() * (1.0 + param_.Nuu()) / 3.0 / a_;

        return p;
    }

    const Scalar initialPressure() const
    { return F_ * param_.B() * (1+param_.Nuu())/3.0/a_;}

    template<class GlobalPosition>
    const GlobalPosition initialDisplacementAtPos(const GlobalPosition& globalPos) const
    {
        GlobalPosition u(0.0);

        u[0] = F_ * param_.Nuu() * globalPos[0] /2.0 /param_.G() / a_;
        u[1] = - F_ * (1-param_.Nuu()) * globalPos[1] /2.0 /param_.G() / a_;
        return u;
    }

    template<class GlobalPosition>
    GlobalPosition displacementAtPos(const GlobalPosition& globalPos,
                                     const Scalar time) const
    {
        GlobalPosition u(0.0);

        using namespace std;

        for (int i = 0; i < nBeta_; i++)
        {
            u[0] += sin(beta_[i]) * cos(beta_[i])
                    / (beta_[i] - sin(beta_[i]) * cos(beta_[i]))
                    * exp(- beta_[i] * beta_[i] * c_ * time / a_ /a_);


            u[1] += sin(beta_[i]) * cos(beta_[i])
                    / (beta_[i] - sin(beta_[i]) * cos(beta_[i]))
                    * exp(- beta_[i] * beta_[i] * c_ * time / a_ /a_);
        }

        u[0] *= -F_*param_.Nuu()/param_.G()/a_;
        u[1] *= F_*(1-param_.Nuu())/param_.G()/a_;

        u[0] += F_ * param_.Nu() / 2.0 / param_.G() /a_;
        u[1] -= F_ * (1-param_.Nu()) / 2.0 / param_.G() / a_;

        u[0] *= globalPos[0];
        u[1] *= globalPos[1];

        Scalar part2 = 0.0;
        for (int i = 0; i < nBeta_; i++)
        {
            part2 += cos(beta_[i])
                     / (beta_[i] - sin(beta_[i]) * cos(beta_[i]))
                     * sin(beta_[i] * globalPos[0] / a_)
                     * exp(- beta_[i] * beta_[i] * c_ * time / a_ /a_);
        }

        part2 *= F_/param_.G();

        u[0] += part2;
        return u;
    }

    template<class GlobalPosition>
    Scalar sigmaZZAtPos(const GlobalPosition& globalPos,
                        const Scalar time) const
    {
        Scalar sigmaZZ = - F_ /a_;
        Scalar part1 = 0.0;
        Scalar part2 = 0.0;
        for (int i = 0; i < nBeta_; i++)
        {
            part1 += sin(beta_[i])
                     /(beta_[i] - sin(beta_[i]) *cos(beta_[i]))
                     * cos(beta_[i] * globalPos[0] / a_)
                     * exp(-beta_[i] * beta_[i] * c_ * time / a_ / a_);

            part2 += sin(beta_[i]) * cos(beta_[i])
                     /(beta_[i] - sin(beta_[i]) *cos(beta_[i]))
                     * exp(-beta_[i] * beta_[i] * c_ * time / a_ / a_);
        }

        part1 *=  - 2.0 * F_ * param_.B() * (param_.Nuu()-param_.Nu()) * (1 + param_.Nuu())
                  / a_ / (1- param_.Nu()) /  (1 + 2 * param_.Nuu());

        part2 *= 2.0 * F_ / a_;

        sigmaZZ += part1;
        sigmaZZ += part2;

        return sigmaZZ;
    }

    const Scalar& alpha() const
    { return param_.alpha();}

    const auto& param() const
    { return param_;}

    const Scalar initialDivU() const
    {
        return F_ * (2 * param_.Nuu() -1) / 2 / param_.G() /a_;
    }

    const Scalar Kdr() const
    { return 2*param_.lambda();}
private:
    Scalar F_; // force intensity
    Scalar a_; // the length of the domain
    int nBeta_;
    std::vector<Scalar> beta_;

    Scalar c_; //compression coefficient
    Scalar A1_;
    Scalar A2_;
    MandelPoroElasticParameters<Scalar> param_;
}; // end class MandelAnalyticalSolution

}//end namespace Dumux

#endif //TEST_MANDEL_ANALYTICAL_SOLUTION_HH
