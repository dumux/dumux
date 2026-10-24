// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup TwoPOneCTests
 * \brief Single-component two-phase fluid system with the assumptions of the heat pipe solution of Udell (1985).
 */
#ifndef DUMUX_TEST_2P1C_HEATPIPE_FLUIDSYSTEM_HH
#define DUMUX_TEST_2P1C_HEATPIPE_FLUIDSYSTEM_HH

#include <cmath>

#include <dumux/common/parameters.hh>
#include <dumux/material/constants.hh>
#include <dumux/material/components/constant.hh>
#include <dumux/material/fluidsystems/2p1c.hh>

namespace Dumux::FluidSystems {

/*!
 * \ingroup TwoPOneCTests
 * \brief Constant densities and viscosities of both phases, a vapor enthalpy that exceeds the liquid
 *        enthalpy by the constant latent heat, and the saturation temperature from the Clapeyron
 *        equation with Kelvin's vapor-pressure lowering, eq. (19) in Udell (1985) \cite Udell1985HeatPipe.
 */
template<class Scalar>
class HeatPipe
: public TwoPOneC<Scalar, Components::Constant<1, Scalar>>
{
    using ParentType = TwoPOneC<Scalar, Components::Constant<1, Scalar>>;
    using Component = Components::Constant<1, Scalar>;

public:
    using typename ParentType::ParameterCache;
    using ParentType::liquidPhaseIdx;
    using ParentType::gasPhaseIdx;

    /*!
     * \brief The temperature of the two-phase state from the vapor pressure \f$p_g\f$ and the capillary pressure \f$p_c\f$,
     *        \f$ T = T_0 (1 + p_c/(h_{fg} \varrho_l)) / (1 - T_0 R_w/h_{fg} \ln(p_g/p_0)) \f$.
     */
    template <class FluidState>
    static Scalar vaporTemperature(const FluidState& fluidState, const unsigned int phaseIdx)
    {
        static const Scalar referenceTemperature = getParam<Scalar>("FluidSystem.ReferenceTemperature");
        static const Scalar referencePressure = getParam<Scalar>("FluidSystem.ReferencePressure");
        const Scalar latentHeat = Component::vaporizationEnthalpy();
        const Scalar liquidDensity = Component::liquidDensity(referenceTemperature, referencePressure);
        const Scalar specificGasConstant = Constants<Scalar>::R/Component::molarMass();

        const Scalar pg = fluidState.pressure(gasPhaseIdx);
        const Scalar pc = pg - fluidState.pressure(liquidPhaseIdx);

        using std::log;
        return referenceTemperature*(1.0 + pc/(latentHeat*liquidDensity))
               /(1.0 - referenceTemperature*specificGasConstant/latentHeat*log(pg/referencePressure));
    }

    template <class FluidState>
    static Scalar enthalpy(const FluidState& fluidState, const int phaseIdx)
    {
        static const Scalar referenceTemperature = getParam<Scalar>("Component.ReferenceTemperature");
        const Scalar temperature = fluidState.temperature(phaseIdx);
        const Scalar sensible = Component::liquidHeatCapacity(temperature, fluidState.pressure(phaseIdx))
                                *(temperature - referenceTemperature);
        return phaseIdx == gasPhaseIdx ? sensible + Component::vaporizationEnthalpy() : sensible;
    }

    template <class FluidState>
    static Scalar enthalpy(const FluidState& fluidState, const ParameterCache&, const int phaseIdx)
    { return enthalpy(fluidState, phaseIdx); }
};

} // end namespace Dumux::FluidSystems

#endif
