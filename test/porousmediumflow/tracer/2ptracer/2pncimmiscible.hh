// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   See the file COPYING for full copying permissions.                      *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 2 of the License, or       *
 *   (at your option) any later version.                                     *
 *                                                                           *
 *   This program is distributed in the hope that it will be useful,         *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the            *
 *   GNU General Public License for more details.                            *
 *                                                                           *
 *   You should have received a copy of the GNU General Public License       *
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.   *
 *****************************************************************************/
/*!
 * \file
 * \ingroup Fluidsystems
 * \brief @copybrief Dumux::FluidSystems::TwoPNCImmiscible
 */
#ifndef DUMUX_2P_NC_IMMISCIBLE_FLUID_SYSTEM_HH
#define DUMUX_2P_NC_IMMISCIBLE_FLUID_SYSTEM_HH

#include <limits>
#include <cassert>

#include <dune/common/exceptions.hh>

#include <dumux/material/fluidsystems/1pliquid.hh>
#include <dumux/material/fluidsystems/1pgas.hh>
#include <dumux/material/fluidstates/immiscible.hh>
#include <dumux/material/components/base.hh>
#include <dumux/material/components/nacl.hh>

#include <dumux/material/fluidsystems/base.hh>

namespace Dumux {
namespace FluidSystems {

/*!
 * \ingroup Fluidsystems
 * \brief A fluid system for two-phase models assuming immiscibility and
 *        thermodynamic equilibrium
 *
 * The fluid phases are completely specified by means of their
 * constituting components.
 * The fluids can be defined individually via FluidSystem::OnePLiquid<Scalar, Component> and
 * FluidSystem::OnePGas<Scalar, Component>. These fluids consist of one component.
 * \tparam Scalar the scalar type
 * \tparam Fluid0 a one-phase fluid system (use FluidSystem::OnePLiquid<Scalar, Component> / FluidSystem::OnePGas<Scalar, Component>)
 * \tparam Fluid1 a one-phase fluid system (use FluidSystem::OnePLiquid<Scalar, Component> / FluidSystem::OnePGas<Scalar, Component>)
 */
template <class Scalar, class Fluid0, class Fluid1>
class TwoPNCImmiscible
: public Base<Scalar, TwoPNCImmiscible<Scalar, Fluid0, Fluid1> >
{
    static_assert((Fluid0::numPhases == 1), "Fluid0 has more than one phase");
    static_assert((Fluid1::numPhases == 1), "Fluid1 has more than one phase");
    static_assert((Fluid0::numComponents == 1), "Fluid0 has more than one component");
    static_assert((Fluid1::numComponents == 1), "Fluid1 has more than one component");
    // two gaseous phases at once do not make sense physically! (but two liquids are fine)
    static_assert(!Fluid0::isGas() || !Fluid1::isGas(), "One phase has to be a liquid!");

    using ThisType = TwoPNCImmiscible<Scalar, Fluid0, Fluid1>;
    using Base = Dumux::FluidSystems::Base<Scalar, ThisType>;
    using NaCl = Components::NaCl<Scalar>;

public:
    static constexpr int numPhases = 2; //!< Number of phases in the fluid system
    // TODO: here, set the number of tracers, twhich will be added to the total number of components automatically:
    static constexpr int numTracers = 1;
    static constexpr int numComponents = 2 + numTracers; //!< Number of components in the fluid system

    static constexpr int phase0Idx = 0; //!< index of the first phase
    static constexpr int phase1Idx = 1; //!< index of the second phase
    static constexpr int comp0Idx = 0; //!< index of the frist component
    static constexpr int comp1Idx = 1; //!< index of the second component

    //TODO: for each tracercomponent, define a unique component index for all added tracers
    static constexpr int comp2Idx = 2; //!< index of the first tracer
//     static constexpr int comp3Idx = 3; //!< index of the second tracer
//     static constexpr int comp4Idx = 4; //!< index of the ... tracer
//     static constexpr int comp5Idx = 5;
//     static constexpr int comp6Idx = 6;
//     static constexpr int comp7Idx = 7;
//     static constexpr int comp8Idx = 8;
//     static constexpr int comp9Idx = 9;
//     static constexpr int comp10Idx = 10;
//     static constexpr int comp11Idx = 11; // index of the 10th tracer

    /****************************************
     * Fluid phase related static parameters
     ****************************************/
    /*!
     * \brief Return the human readable name of a fluid phase
     * \param phaseIdx The index of the fluid phase to consider
     */
    static std::string phaseName(int phaseIdx)
    {
        assert(0 <= phaseIdx && phaseIdx < numPhases);

        if (!Fluid0::isGas() && !Fluid1::isGas())
        {
            if (phaseIdx == phase0Idx)
                return Components::IsAqueous<typename Fluid0::Component>::value ? "aq" : "napl";
            else
                return Components::IsAqueous<typename Fluid1::Component>::value ? "aq" : "napl";
        }
        else
        {
            if (phaseIdx == phase0Idx)
                return Fluid0::isGas() ? "gas" : "liq";
            else
                return Fluid1::isGas() ? "gas" : "liq";
        }
    }

    /*!
     * \brief Returns whether the fluids are miscible
     */
    static constexpr bool isMiscible()
    { return false; }

    /*!
     * \brief Return whether a phase is gaseous
     * \param phaseIdx The index of the fluid phase to consider
     */
    static constexpr bool isGas(int phaseIdx)
    {
        assert(0 <= phaseIdx && phaseIdx < numPhases);

        if (phaseIdx == phase0Idx)
            return Fluid0::isGas();
            return Fluid1::isGas();
    }

    /*!
     * \brief Returns true if and only if a fluid phase is assumed to
     *        be an ideal mixture.
     * \param phaseIdx The index of the fluid phase to consider
     *
     * We define an ideal mixture as a fluid phase where the fugacity
     * coefficients of all components times the pressure of the phase
     * are independent on the fluid composition. This assumption is true
     * if immiscibility is assumed. If you are unsure what
     * this function should return, it is safe to return false. The
     * only damage done will be (slightly) increased computation times
     * in some cases.
     */
    static bool isIdealMixture(int phaseIdx)
    {
        assert(0 <= phaseIdx && phaseIdx < numPhases);

        // we assume immisibility
        return true;
    }

    /*!
     * \brief Returns true if and only if a fluid phase is assumed to
     *        be an ideal gas.
     *
     * \param phaseIdx The index of the fluid phase to consider
     */
    static constexpr bool isIdealGas(int phaseIdx)
    {
        assert(0 <= phaseIdx && phaseIdx < numPhases);

        // let the fluids decide
        if (phaseIdx == phase0Idx)
            return Fluid0::isIdealGas();
        return Fluid1::isIdealGas();
    }

    /*!
     * \brief Returns true if and only if a fluid phase is assumed to
     *        be compressible.
     *
     * Compressible means. that the partial derivative of the density
     * to the fluid pressure is always larger than zero.
     *
     * \param phaseIdx The index of the fluid phase to consider
     */
    static constexpr bool isCompressible(int phaseIdx)
    {
        assert(0 <= phaseIdx && phaseIdx < numPhases);

        // let the fluids decide
        if (phaseIdx == phase0Idx)
            return Fluid0::isCompressible();
            return Fluid1::isCompressible();
    }

    /*!
     * \brief Returns true if the liquid phase viscostiy is constant
     *
     * \param phaseIdx The index of the fluid phase to consider
     */
    static constexpr bool viscosityIsConstant(int phaseIdx)
    {
        assert(0 <= phaseIdx && phaseIdx < numPhases);

        // let the fluids decide
        if (phaseIdx == phase0Idx)
            return Fluid0::viscosityIsConstant();
            return Fluid1::viscosityIsConstant();
    }

    /*!
     * \brief Returns true if and only if a fluid phase is assumed to
     *        be an ideal gas.
     *
     * \param phaseIdx The index of the fluid phase to consider
     */
    static bool isIdealFluid1(int phaseIdx)
    {
        assert(0 <= phaseIdx && phaseIdx < numPhases);

        // let the fluids decide
        if (phaseIdx == phase0Idx)
            return Fluid0::isIdealFluid1();
            return Fluid1::isIdealFluid1();
    }

    /****************************************
     * Component related static parameters
     ****************************************/
    /*!
     * \brief Return the human readable name of a component
     *
     * \param compIdx index of the component
     */
    static std::string componentName(int compIdx)
    {
        assert(0 <= compIdx && compIdx < numComponents);

        if (compIdx == comp0Idx)
            return Fluid0::name();
        else if (compIdx == comp1Idx)
            return Fluid1::name();
        //TODO: for each tracercomponent, define a unique name! Otherwise Paraview will crash.
        else if (compIdx == comp2Idx)
            return "NaCl1";
//         else if (compIdx == comp3Idx)
//             return "NaCl2";
//         else if (compIdx == comp4Idx)
//             return "NaCl3";
//         else if (compIdx == comp5Idx)
//             return "NaCl4";
//         else if (compIdx == comp5Idx)
//             return "NaCl5";
//         else if (compIdx == comp5Idx)
//             return "NaCl6";
//         else if (compIdx == comp5Idx)
//             return "NaCl7";
//         else if (compIdx == comp5Idx)
//             return "NaCl8";
//         else if (compIdx == comp5Idx)
//             return "NaCl9";
//         else if (compIdx == comp5Idx)
//             return "NaCl10";
    }


    /*!
     * \brief Return the molar mass of a component in \f$\mathrm{[kg/mol]}\f$.
     * \param compIdx index of the component
     */
    static Scalar molarMass(int compIdx)
    {
        assert(0 <= compIdx && compIdx < numComponents);

        if (compIdx == comp0Idx)
            return Fluid0::molarMass();
        else if (compIdx == comp1Idx)
            return Fluid1::molarMass();
        //else, tracer, so returning NaCl value
        else // (compIdx > comp1Idx)
            return NaCl::molarMass();
    }

    /*!
     * \brief Critical temperature of a component \f$\mathrm{[K]}\f$.
     * \param compIdx index of the component
     */
    static Scalar criticalTemperature(int compIdx)
    {
        assert(0 <= compIdx && compIdx < numComponents);

        if (compIdx == comp0Idx)
            return Fluid0::criticalTemperature();
        return Fluid1::criticalTemperature();
    }

    /*!
     * \brief Critical pressure of a component \f$\mathrm{[Pa]}\f$.
     * \param compIdx index of the component
     */
    static Scalar criticalPressure(int compIdx)
    {
        assert(0 <= compIdx && compIdx < numComponents);

        if (compIdx == comp0Idx)
            return Fluid0::criticalPressure();
        return Fluid1::criticalPressure();
    }

    /*!
     * \brief The acentric factor of a component \f$\mathrm{[-]}\f$.
     * \param compIdx index of the component
     */
    static Scalar acentricFactor(int compIdx)
    {
        assert(0 <= compIdx && compIdx < numComponents);

        if (compIdx == comp0Idx)
            return Fluid0::acentricFactor();
        return Fluid1::acentricFactor();
    }

    /****************************************
     * thermodynamic relations
     ****************************************/

    /*!
     * \brief Initialize the fluid system's static parameters
     */
    static void init()
    {
        // initialize with some default values
        init(/*tempMin=*/273.15, /*tempMax=*/623.15, /*numTemp=*/100,
             /*pMin=*/-10.0, /*pMax=*/20e6, /*numP=*/200);
    }

    /*!
     * \brief Initialize the fluid system's static parameters using
     *        problem specific temperature and pressure ranges
     *
     * \param tempMin The minimum temperature used for tabulation of water \f$\mathrm{[K]}\f$
     * \param tempMax The maximum temperature used for tabulation of water\f$\mathrm{[K]}\f$
     * \param nTemp The number of ticks on the temperature axis of the  table of water
     * \param pressMin The minimum pressure used for tabulation of water \f$\mathrm{[Pa]}\f$
     * \param pressMax The maximum pressure used for tabulation of water \f$\mathrm{[Pa]}\f$
     * \param nPress The number of ticks on the pressure axis of the  table of water
     */
    static void init(Scalar tempMin, Scalar tempMax, std::size_t nTemp,
                     Scalar pressMin, Scalar pressMax, std::size_t nPress)
    {
        if (Fluid0::Component::isTabulated)
            Fluid0::Component::init(tempMin, tempMax, nTemp, pressMin, pressMax, nPress);

        if (Fluid1::Component::isTabulated)
            Fluid1::Component::init(tempMin, tempMax, nTemp, pressMin, pressMax, nPress);
    }

    using Base::density;
    /*!
     * \brief Calculate the density \f$\mathrm{[kg/m^3]}\f$ of a fluid phase
     *
     */
    template <class FluidState>
    static Scalar density(const FluidState &fluidState,
                          int phaseIdx)
    {
        assert(0 <= phaseIdx && phaseIdx < numPhases);

        Scalar temperature = fluidState.temperature(phaseIdx);
        Scalar pressure = fluidState.pressure(phaseIdx);
        if (phaseIdx == phase0Idx)
            return Fluid0::density(temperature, pressure);
        return Fluid1::density(temperature, pressure);
    }

    using Base::molarDensity;
    /*!
     * \brief The molar density \f$\rho_{mol,\alpha}\f$
     *   of a fluid phase \f$\alpha\f$ in \f$\mathrm{[mol/m^3]}\f$
     *
     * The molar density is defined by the
     * mass density \f$\rho_\alpha\f$ and the component molar mass \f$M_\alpha\f$:
     *
     * \f[\rho_{mol,\alpha} = \frac{\rho_\alpha}{M_\alpha} \;.\f]
     */
    template <class FluidState>
    static Scalar molarDensity(const FluidState &fluidState, int phaseIdx)
    {
       Scalar temperature = fluidState.temperature(phaseIdx);
       Scalar pressure = fluidState.pressure(phaseIdx);
       if (phaseIdx == phase0Idx)
           return Fluid0::molarDensity(temperature, pressure);
       return Fluid1::molarDensity(temperature, pressure);
    }

    using Base::viscosity;
    /*!
     * \brief Return the viscosity of a phase \f$\mathrm{[Pa*s]}\f$.
     * \param fluidState The fluid state of the two-phase model
     * \param phaseIdx Index of the fluid phase
     */
    template <class FluidState>
    static Scalar viscosity(const FluidState &fluidState,
                            int phaseIdx)
    {
        assert(0 <= phaseIdx  && phaseIdx < numPhases);

        Scalar temperature = fluidState.temperature(phaseIdx);
        Scalar pressure = fluidState.pressure(phaseIdx);
        if (phaseIdx == phase0Idx)
            return Fluid0::viscosity(temperature, pressure);
        return Fluid1::viscosity(temperature, pressure);
    }

    using Base::fugacityCoefficient;
    /*!
     * \brief Calculate the fugacity coefficient \f$\mathrm{[-]}\f$ of an individual
     *        component in a fluid phase
     *
     * The fugacity coefficient \f$\mathrm{\phi^\kappa_\alpha}\f$ is connected to the
     * fugacity \f$\mathrm{f^\kappa_\alpha}\f$ and the component's mole
     * fraction \f$\mathrm{x^\kappa_\alpha}\f$ by means of the relation
     *
     * \f[
     f^\kappa_\alpha = \phi^\kappa_\alpha\;x^\kappa_\alpha\;p_\alpha
     * \f]
     *
     * \param fluidState The fluid state of the two-phase model
     * \param phaseIdx Index of the fluid phase
     * \param compIdx index of the component
     */
    template <class FluidState>
    static Scalar fugacityCoefficient(const FluidState &fluidState,
                                      int phaseIdx,
                                      int compIdx)
    {
        assert(0 <= phaseIdx  && phaseIdx < numPhases);
        assert(0 <= compIdx  && compIdx < numComponents);

        if (phaseIdx == phase1Idx && compIdx == comp1Idx)
            return 1e-2;
            // We could calculate the real fugacity coefficient of
            // the component in the fluid. Probably that's not worth
            // the effort, since the fugacity coefficient of the other
            // component is infinite anyway...
        else if (phaseIdx == phase1Idx && compIdx != comp1Idx)
            return 1.0;
        else if (phaseIdx == phase0Idx && compIdx == comp0Idx)
            //water likes to stay in the water phase
            return 1e-2;
        else if (phaseIdx == phase0Idx && compIdx == comp1Idx)
            //NAPL doesn't want to be in the water phase
            return fluidState.pressure(phase0Idx);
        else if (compIdx > numPhases && phaseIdx == phase0Idx)
            //NaCl/tracer always stays in the water phase
            return 0.0;

        return 1.0;
    }

    using Base::diffusionCoefficient;
    /*!
     * \brief Calculate the binary molecular diffusion coefficient for
     *        a component in a fluid phase \f$\mathrm{[mol^2 * s / (kg*m^3)]}\f$
     * \param fluidState The fluid state of the two-phase model
     * \param phaseIdx Index of the fluid phase
     * \param compIdx index of the component
     *
     * Molecular diffusion of a compoent \f$\mathrm{\kappa}\f$ is caused by a
     * gradient of the chemical potential and follows the law
     *
     * \f[ J = - D \mathbf{grad} \mu_\kappa \f]
     *
     * where \f$\mathrm{\mu_\kappa]}\f$ is the component's chemical potential,
     * \f$\mathrm{D}\f$ is the diffusion coefficient and \f$\mathrm{J}\f$ is the
     * diffusive flux. \f$\mathrm{\mu_\kappa}\f$ is connected to the component's
     * fugacity \f$\mathrm{f_\kappa}\f$ by the relation
     *
     * \f[ \mu_\kappa = R T_\alpha \mathrm{ln} \frac{f_\kappa}{p_\alpha} \f]
     *
     * where \f$\mathrm{p_\alpha}\f$ and \f$\mathrm{T_\alpha}\f$ are the fluid phase'
     * pressure and temperature.
     */
    template <class FluidState>
    static Scalar diffusionCoefficient(const FluidState &fluidState,
                                       int phaseIdx,
                                       int compIdx)
    {
        DUNE_THROW(Dune::InvalidStateException,
                   "Diffusion coefficients of components are meaningless if"
                   " immiscibility is assumed");
    }

    using Base::binaryDiffusionCoefficient;
    /*!
     * \brief Given a phase's composition, temperature and pressure,
     *        return the binary diffusion coefficient \f$\mathrm{[m^2/s]}\f$ for components
     *        \f$\mathrm{i}\f$ and \f$\mathrm{j}\f$ in this phase.
     * \param fluidState The fluid state of the two-phase model
     * \param phaseIdx Index of the fluid phase
     * \param compIIdx index of the component i
     * \param compJIdx index of the component j
     */
    template <class FluidState>
    static Scalar binaryDiffusionCoefficient(const FluidState &fluidState,
                                             int phaseIdx,
                                             int compIIdx,
                                             int compJIdx)

    {
        //TODO: below is the function from brineair, however, in the 2ptracerprblem, diffusion is set to zero, so doing the same here:
        return 0.0;
//         assert(0 <= phaseIdx && phaseIdx < numPhases);
//         assert(0 <= compIIdx && compIIdx < numComponents);
//         assert(0 <= compJIdx && compJIdx < numComponents);
//
//         const auto T = fluidState.temperature(phaseIdx);
//         const auto p = fluidState.pressure(phaseIdx);
//
//         if (compIIdx > compJIdx)
//             std::swap(compIIdx, compJIdx);
//
//         if (phaseIdx == phase0Idx)
//         {
// //             if(compIIdx == H2OIdx && compJIdx == AirIdx)
// //                 return H2O_Air::liquidDiffCoeff(T, p);
//             if (compIIdx == comp0Idx && compJIdx > numPhases)
//                 return Brine::binaryDiffusionCoefficient(BrineAdapter<FluidState>(fluidState), Brine::liquidPhaseIdx, Brine::H2OIdx, Brine::NaClIdx);
//             else
//                 DUNE_THROW(Dune::NotImplemented, "Binary diffusion coefficient of components "
//                                                  << compIIdx << " and " << compJIdx
//                                                  << " in phase " << phaseIdx);
//         }
//         else if (phaseIdx == gasPhaseIdx)
//         {
//             if (compIIdx == H2OIdx && compJIdx == AirIdx)
//                 return H2O_Air::gasDiffCoeff(T, p);
//
//             // NaCl is expected to never be present in the gas phase. we need to
//             // return a diffusion coefficient that does not case numerical problems.
//             // We choose a very small value here.
//             else if (compIIdx == AirIdx && compJIdx == NaClIdx)
//                 return 1e-12;
//
//             else
//                 DUNE_THROW(Dune::NotImplemented, "Binary diffusion coefficient of components "
//                                                  << compIIdx << " and " << compJIdx
//                                                  << " in phase " << phaseIdx);
//         }
//
//         DUNE_THROW(Dune::InvalidStateException, "Invalid phase index " << phaseIdx);
    }

    using Base::enthalpy;
    /*!
     * \brief Return the specific enthalpy of a fluid phase \f$\mathrm{[J/kg]}\f$.
     * \param fluidState The fluid state of the two-phase model
     * \param phaseIdx Index of the fluid phase
     */
    template <class FluidState>
    static Scalar enthalpy(const FluidState &fluidState,
                                 int phaseIdx)
    {
        assert(0 <= phaseIdx  && phaseIdx < numPhases);

        Scalar temperature = fluidState.temperature(phaseIdx);
        Scalar pressure = fluidState.pressure(phaseIdx);
        if (phaseIdx == phase0Idx)
            return Fluid0::enthalpy(temperature, pressure);
        return Fluid1::enthalpy(temperature, pressure);
    }

    using Base::thermalConductivity;
    /*!
     * \brief Thermal conductivity of a fluid phase \f$\mathrm{[W/(m K)]}\f$.
     * \param fluidState The fluid state of the two-phase model
     * \param phaseIdx Index of the fluid phase
     */
    template <class FluidState>
    static Scalar thermalConductivity(const FluidState &fluidState,
                                      int phaseIdx)
    {
        assert(0 <= phaseIdx  && phaseIdx < numPhases);

        Scalar temperature = fluidState.temperature(phaseIdx);
        Scalar pressure = fluidState.pressure(phaseIdx);
        if (phaseIdx == phase0Idx)
            return Fluid0::thermalConductivity(temperature, pressure);
        return Fluid1::thermalConductivity(temperature, pressure);
    }

    using Base::heatCapacity;
    /*!
     * @copybrief Base::thermalConductivity
     *
     * Additional comments:
     *
     * Specific isobaric heat capacity of a fluid phase.
     *        \f$\mathrm{[J/(kg*K)]}\f$.
     *
     * \param fluidState The fluid state of the two-phase model
     * \param phaseIdx for which phase to give back the heat capacity
     */
    template <class FluidState>
    static Scalar heatCapacity(const FluidState &fluidState,
                               int phaseIdx)
    {
        assert(0 <= phaseIdx  && phaseIdx < numPhases);

        Scalar temperature = fluidState.temperature(phaseIdx);
        Scalar pressure = fluidState.pressure(phaseIdx);
        if (phaseIdx == phase0Idx)
            return Fluid0::heatCapacity(temperature, pressure);
        return Fluid1::heatCapacity(temperature, pressure);
    }
};

} // end namespace FluidSystems
} // end namespace Dumux

#endif
