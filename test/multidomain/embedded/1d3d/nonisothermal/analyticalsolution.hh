// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \brief Reference solution of the wellbore heat transport benchmark.
 *
 * Exact solution of the benchmark problem without axial heat conduction in the fluid and in the
 * rock, with an infinite rock in the radial direction, a constant mean velocity
 * v = q/(pi r_i^2) and a constant heat transfer coefficient U_o referred to the outer radius r_o,
 * as in the numerical model. It keeps the heat storage in the fluid and the transient thermal
 * resistance of the rock.
 *
 * The rock temperature is eliminated in the Laplace domain, where the fluid temperature reads
 *
 *   T_f(l, s) = T_s0/s + dT/(s + v kappa) + dT (1/s - 1/(s + v kappa)) exp(-kappa l) exp(-s l/v),
 *   kappa(s)  = 1/(q rho_f c_p_f (R_U + R_s(s))),   R_U = 1/(2 pi r_o U_o),
 *
 * with dT = T_in - T_s0 and the distance l from the inlet. The transient thermal resistance of
 * the rock R_s(s) with rho = r_o sqrt(s/alpha_s) depends on the borehole interior r < r_o:
 *
 *   benchmark: K0(rho)/(2 pi lambda_s rho K1(rho))   rock only at r >= r_o, flux at r = r_o
 *   embedded:  I0(rho) K0(rho)/(2 pi lambda_s)       rock everywhere, source on r = r_o
 *
 * The factor exp(-s l/v) is a shift in time by the travel time l/v. The two remaining transforms
 * are inverted numerically with the fixed Talbot method (Abate & Valko 2004).
 */
#ifndef WELLBORE_ANALYTICAL_SOLUTION_HH
#define WELLBORE_ANALYTICAL_SOLUTION_HH

#include <cmath>
#include <complex>

#include <dumux/common/parameters.hh>

namespace Dumux {

template<class Scalar>
class WellboreReferenceSolution
{
    using Complex = std::complex<Scalar>;

public:
    //! the borehole interior r < r_o: pipe and grout without heat capacity, or filled with rock as in the 1D-3D model
    enum class Variant { benchmark, embedded };

    WellboreReferenceSolution(const Scalar innerRadius, const Scalar outerRadius)
    : outerRadius_(outerRadius)
    {
        const Scalar injectionRate = getParam<Scalar>("Well.Problem.InjectionRate");
        const Scalar fluidDensity = getParam<Scalar>("Component.LiquidDensity");
        const Scalar fluidHeatCapacity = getParam<Scalar>("Component.LiquidHeatCapacity");
        const Scalar fluidConductivity = getParam<Scalar>("Component.LiquidThermalConductivity");
        const Scalar wallConductivity = getParam<Scalar>("MixedDimension.PipeWallThermalConductivity");
        injectionTemperature_ = getParam<Scalar>("Well.Problem.InjectionTemperature");
        rockTemperature_ = getParam<Scalar>("Rock.Problem.InitialTemperature");
        rockConductivity_ = getParam<Scalar>("Component.SolidThermalConductivity");
        rockDiffusivity_ = rockConductivity_/(getParam<Scalar>("Component.SolidDensity")*getParam<Scalar>("Component.SolidHeatCapacity"));

        velocity_ = injectionRate/(M_PI*innerRadius*innerRadius);
        capacityRate_ = injectionRate*fluidDensity*fluidHeatCapacity;

        // the heat transfer coefficient of ExchangeFluxCalculator, referred to the outer radius
        const Scalar nusselt = 4.364;
        const Scalar film = fluidConductivity*nusselt/(2*innerRadius);
        const Scalar transferCoefficient = 1.0/(outerRadius/(innerRadius*film) + outerRadius*std::log(outerRadius/innerRadius)/wallConductivity);
        pipeResistance_ = 1.0/(2*M_PI*outerRadius*transferCoefficient);
    }

    //! the time [s] the fluid needs from the inlet to the distance l [m]
    Scalar travelTime(const Scalar l) const
    { return l/velocity_; }

    //! the fluid temperature [K] at the distance l [m] from the inlet at the time t > 0 [s]
    Scalar temperature(const Scalar l, const Scalar t, const Variant variant) const
    {
        const auto kappa = [&](const Complex& s)
        { return 1.0/(capacityRate_*(pipeResistance_ + rockResistance_(s, variant))); };

        // the fluid that was in the pipe at t = 0, heated in place
        const Scalar initial = talbotInverse_([&](const Complex& s)
        { return 1.0/(s + velocity_*kappa(s)); }, t);

        // the fluid that entered at the inlet. Its contribution starts from zero at the arrival,
        // so a difference in the order of round-off errors must not reach the inverse transform
        const Scalar delay = t - travelTime(l);
        const Scalar inflow = delay < 1e-3 ? 0.0 : talbotInverse_([&](const Complex& s)
        {
            const Complex k = kappa(s);
            return (1.0/s - 1.0/(s + velocity_*k))*std::exp(-k*l);
        }, delay);

        return rockTemperature_ + (injectionTemperature_ - rockTemperature_)*(initial + inflow);
    }

private:
    //! the thermal resistance of the rock per unit length in the Laplace domain, I0 K1 + I1 K0 = 1/rho gives the embedded variant
    Complex rockResistance_(const Complex& s, const Variant variant) const
    {
        const Complex rho = outerRadius_*std::sqrt(s/rockDiffusivity_);
        const Complex ratio = variant == Variant::benchmark ? ratioK1K0_(rho) : ratioK1K0_(rho) + ratioI1I0_(rho);
        return 1.0/(2*M_PI*rockConductivity_*rho*ratio);
    }

    //! K1(z)/K0(z) for Re(z) > 0, with Steed's method for the continued fraction of Temme (Numerical Recipes, Bessel functions of fractional order)
    static Complex ratioK1K0_(const Complex& z)
    {
        Complex b = 2.0*(1.0 + z);
        Complex d = 1.0/b;
        Complex delta = d;
        Complex h = d;
        Scalar a = -0.25;
        for (int i = 2; std::abs(delta) > 1e-15*std::abs(h); ++i)
        {
            a -= 2.0*(i - 1);
            b += 2.0;
            d = 1.0/(b + a*d);
            delta *= b*d - 1.0;
            h += delta;
        }
        return (0.5 + z - 0.25*h)/z;
    }

    //! I1(z)/I0(z), backward evaluation of the continued fraction from I_{n-1} = I_{n+1} + 2n/z I_n
    static Complex ratioI1I0_(const Complex& z)
    {
        Complex ratio(0.0);
        for (int n = static_cast<int>(2*std::abs(z)) + 50; n > 0; --n)
            ratio = 1.0/(2.0*n/z + ratio);
        return ratio;
    }

    //! the inverse Laplace transform at the time t > 0 with the fixed Talbot method
    template<class Transform>
    static Scalar talbotInverse_(const Transform& transform, const Scalar t)
    {
        constexpr int numNodes = 24;
        const Scalar r = 2.0*numNodes/(5.0*t);
        Complex sum = 0.5*std::exp(r*t)*transform(Complex(r));
        for (int k = 1; k < numNodes; ++k)
        {
            const Scalar theta = M_PI*k/numNodes;
            const Scalar cot = 1.0/std::tan(theta);
            const Complex s = r*theta*Complex(cot, 1.0);
            sum += std::exp(t*s)*transform(s)*Complex(1.0, theta + (theta*cot - 1.0)*cot);
        }
        return r/numNodes*sum.real();
    }

    Scalar outerRadius_;
    Scalar injectionTemperature_, rockTemperature_;
    Scalar rockConductivity_, rockDiffusivity_;
    Scalar velocity_, capacityRate_, pipeResistance_;
};

} // end namespace Dumux

#endif
