// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup TwoPOneCTests
 * \brief Steady solution of the horizontal heat pipe of Udell (1985) inside the two-phase zone.
 */
#ifndef DUMUX_TEST_2P1C_HEATPIPE_SOLUTION_HH
#define DUMUX_TEST_2P1C_HEATPIPE_SOLUTION_HH

#include <cmath>

#include <dune/common/exceptions.hh>

#include <dumux/common/parameters.hh>
#include <dumux/common/integrate.hh>
#include <dumux/material/constants.hh>

namespace Dumux {

/*!
 * \ingroup TwoPOneCTests
 * \brief Steady solution of the horizontal heat pipe of Udell (1985) \cite Udell1985HeatPipe inside the two-phase zone.
 *
 * Without heat conduction, the vapor mass flux \f$ q/h_{fg} \f$ flows towards the cold end at \f$ x = 0 \f$
 * and the same liquid mass flux flows back. Darcy's law for both phases gives the capillary pressure
 * gradient (eq. (9) of Udell)
 * \f[ \frac{\mathrm{d} p_c}{\mathrm{d} x} = \frac{q}{K h_{fg}} \left( \frac{\nu_g}{k_{rg}} + \frac{\nu_l}{k_{rl}} \right), \f]
 * with \f$ p_c = \sigma \sqrt{\phi/K}\, f(s) \f$. The position is a quadrature in the saturation,
 * the gas pressure follows from \f$ \mathrm{d} p_g / \mathrm{d} s = \sigma \sqrt{\phi/K}\, f'(s) / (1 + \beta k_{rg}/k_{rl}) \f$
 * with \f$ \beta = \nu_l/\nu_g \f$, and the temperature from the Clapeyron equation with Kelvin's
 * vapor-pressure lowering (eq. (19) of Udell).
 */
template<class Scalar>
class TwoPOneCHeatPipeSolution
{
public:
    TwoPOneCHeatPipeSolution()
    {
        const Scalar permeability = getParam<Scalar>("SpatialParams.Permeability");
        const Scalar porosity = getParam<Scalar>("SpatialParams.Porosity");
        pcScale_ = getParam<Scalar>("SpatialParams.SurfaceTension")*std::sqrt(porosity/permeability);

        latentHeat_ = getParam<Scalar>("Component.EnthalpyOfVaporization");
        liquidDensity_ = getParam<Scalar>("Component.LiquidDensity");
        nuL_ = getParam<Scalar>("Component.LiquidDynamicViscosity")/liquidDensity_;
        nuG_ = getParam<Scalar>("Component.GasDynamicViscosity")/getParam<Scalar>("Component.GasDensity");
        specificGasConstant_ = Constants<Scalar>::R/getParam<Scalar>("Component.MolarMass");

        fluxScale_ = getParam<Scalar>("Problem.HeatFlux")/(permeability*latentHeat_);
        coldEndGasPressure_ = getParam<Scalar>("Problem.ColdEndGasPressure");
        coldEndSaturation_ = getParam<Scalar>("Problem.ColdEndSaturation");

        referenceTemperature_ = getParam<Scalar>("FluidSystem.ReferenceTemperature");
        referencePressure_ = getParam<Scalar>("FluidSystem.ReferencePressure");
    }

    //! The liquid saturation at the distance x from the cold end
    Scalar saturation(const Scalar x) const
    {
        if (x <= 0.0)
            return coldEndSaturation_;

        Scalar lower = minSaturation_;
        Scalar upper = coldEndSaturation_;
        if (position_(lower) < x)
            DUNE_THROW(Dune::InvalidStateException, "Position " << x << " is beyond the end of the two-phase zone");

        while (upper - lower > 1e-14)
        {
            const Scalar middle = 0.5*(lower + upper);
            if (position_(middle) > x)
                lower = middle;
            else
                upper = middle;
        }
        return 0.5*(lower + upper);
    }

    struct Values
    {
        Scalar saturation;
        Scalar gasPressure;
        Scalar liquidPressure;
        Scalar temperature;
    };

    //! Liquid saturation, phase pressures and temperature at the distance x from the cold end
    Values operator() (const Scalar x) const
    {
        const Scalar s = saturation(x);
        const Scalar pc = pcScale_*f_(s);
        const auto dpgds = [&](const Scalar u){ return pcScale_*dfds_(u)/(1.0 + nuL_/nuG_*krg_(u)/krl_(u)); };
        const Scalar pg = coldEndGasPressure_ - integrateScalarFunction(dpgds, s, coldEndSaturation_, 1e-12*pcScale_);

        using std::log;
        const Scalar temperature = referenceTemperature_*(1.0 + pc/(latentHeat_*liquidDensity_))
                                   /(1.0 - referenceTemperature_*specificGasConstant_/latentHeat_*log(pg/referencePressure_));
        return {s, pg, pg - pc, temperature};
    }

private:
    Scalar position_(const Scalar s) const
    {
        const auto dxds = [&](const Scalar u){ return -pcScale_*dfds_(u)/(fluxScale_*(nuG_/krg_(u) + nuL_/krl_(u))); };
        return integrateScalarFunction(dxds, s, coldEndSaturation_, 1e-14);
    }

    // Leverett function of Udell's eq. (20) and its derivative
    static Scalar f_(const Scalar s)
    {
        const Scalar sn = 1.0 - s;
        return ((1.263*sn - 2.120)*sn + 1.417)*sn;
    }

    static Scalar dfds_(const Scalar s)
    {
        const Scalar sn = 1.0 - s;
        return -((3.0*1.263*sn - 2.0*2.120)*sn + 1.417);
    }

    static Scalar krl_(const Scalar s) { return s*s*s; }
    static Scalar krg_(const Scalar s) { return (1.0 - s)*(1.0 - s)*(1.0 - s); }

    static constexpr Scalar minSaturation_ = 1e-6;
    Scalar pcScale_, latentHeat_, liquidDensity_, nuL_, nuG_, specificGasConstant_;
    Scalar fluxScale_, coldEndGasPressure_, coldEndSaturation_;
    Scalar referenceTemperature_, referencePressure_;
};

} // end namespace Dumux

#endif
