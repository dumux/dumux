// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup FluidSystems
 * \copybrief Dumux::FluidSystems::TwoPImmiscible
 */
#ifndef DUMUX_2P_IMMISCIBLE_FLUID_SYSTEM_HH
#define DUMUX_2P_IMMISCIBLE_FLUID_SYSTEM_HH

#include <limits>
#include <cassert>

#include <dune/common/exceptions.hh>

#include <dumux/material/fluidsystems/1pliquid.hh>
#include <dumux/material/fluidsystems/1pgas.hh>
#include <dumux/material/fluidstates/immiscible.hh>
#include <dumux/material/components/base.hh>
#include <dumux/io/name.hh>

#include "base.hh"

namespace Dumux {
namespace FluidSystems {

/*!
 * \ingroup FluidSystems
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
class TwoPImmiscible
: public Base<Scalar, TwoPImmiscible<Scalar, Fluid0, Fluid1> >
{
    static_assert((Fluid0::numPhases == 1), "Fluid0 has more than one phase");
    static_assert((Fluid1::numPhases == 1), "Fluid1 has more than one phase");
    static_assert((Fluid0::numComponents == 1), "Fluid0 has more than one component");
    static_assert((Fluid1::numComponents == 1), "Fluid1 has more than one component");
    // two gaseous phases at once do not make sense physically! (but two liquids are fine)
    static_assert(!Fluid0::isGas() || !Fluid1::isGas(), "One phase has to be a liquid!");

    using ThisType = TwoPImmiscible<Scalar, Fluid0, Fluid1>;

public:
    static constexpr int numPhases = 2; //!< Number of phases in the fluid system
    static constexpr int numComponents = 2; //!< Number of components in the fluid system

    static constexpr int phase0Idx = 0; //!< index of the first phase
    static constexpr int phase1Idx = 1; //!< index of the second phase
    static constexpr int comp0Idx = 0; //!< index of the first component
    static constexpr int comp1Idx = 1; //!< index of the second component

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
            static const auto name0 = Components::IsAqueous<typename Fluid0::Component>::value ? IOName::aqueousPhase() : IOName::naplPhase();
            static const auto name1 = Components::IsAqueous<typename Fluid1::Component>::value ? IOName::aqueousPhase() : IOName::naplPhase();

            if (name0 != name1)
                return (phaseIdx == phase0Idx) ? name0 : name1;
            else
                return (phaseIdx == phase0Idx) ? name0 + "_0" : name1 + "_1";
        }
        else
        {
            if (phaseIdx == phase0Idx)
                return Fluid0::isGas() ? IOName::gaseousPhase() : IOName::liquidPhase();
            else
                return Fluid1::isGas() ? IOName::gaseousPhase() : IOName::liquidPhase();
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
        return Fluid1::name();
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
        return Fluid1::molarMass();
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

    using Base<Scalar, ThisType>::density;
    //! \copybrief Base<Scalar,ThisType>::density(const FluidState&,int)
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

    using Base<Scalar, ThisType>::molarDensity;
    //! \copybrief Base<Scalar,ThisType>::molarDensity(const FluidState&,int)
    template <class FluidState>
    static Scalar molarDensity(const FluidState &fluidState, int phaseIdx)
    {
       Scalar temperature = fluidState.temperature(phaseIdx);
       Scalar pressure = fluidState.pressure(phaseIdx);
       if (phaseIdx == phase0Idx)
           return Fluid0::molarDensity(temperature, pressure);
       return Fluid1::molarDensity(temperature, pressure);
    }

    using Base<Scalar, ThisType>::viscosity;
    //! \copybrief Base<Scalar,ThisType>::viscosity(const FluidState&,int)
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

    using Base<Scalar, ThisType>::fugacityCoefficient;
    //! \copybrief Base<Scalar,ThisType>::fugacityCoefficient(const FluidState&,int,int)
    template <class FluidState>
    static Scalar fugacityCoefficient(const FluidState &fluidState,
                                      int phaseIdx,
                                      int compIdx)
    {
        assert(0 <= phaseIdx  && phaseIdx < numPhases);
        assert(0 <= compIdx  && compIdx < numComponents);

        if (phaseIdx == compIdx)
            // We could calculate the real fugacity coefficient of
            // the component in the fluid. Probably that's not worth
            // the effort, since the fugacity coefficient of the other
            // component is infinite anyway...
            return 1.0;
        return std::numeric_limits<Scalar>::infinity();
    }

    using Base<Scalar, ThisType>::diffusionCoefficient;
    //! \copybrief Base<Scalar,ThisType>::diffusionCoefficient(const FluidState&,int,int)
    template <class FluidState>
    static Scalar diffusionCoefficient(const FluidState &fluidState,
                                       int phaseIdx,
                                       int compIdx)
    {
        DUNE_THROW(Dune::InvalidStateException,
                   "Diffusion coefficients of components are meaningless if"
                   " immiscibility is assumed");
    }

    using Base<Scalar, ThisType>::binaryDiffusionCoefficient;
    //! \copybrief Base<Scalar,ThisType>::binaryDiffusionCoefficient(const FluidState&,int,int,int)
    template <class FluidState>
    static Scalar binaryDiffusionCoefficient(const FluidState &fluidState,
                                             int phaseIdx,
                                             int compIIdx,
                                             int compJIdx)

    {
        DUNE_THROW(Dune::InvalidStateException,
                   "Binary diffusion coefficients of components are meaningless if"
                   " immiscibility is assumed");
    }

    using Base<Scalar, ThisType>::enthalpy;
    //! \copybrief Base<Scalar,ThisType>::enthalpy(const FluidState&,int)
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

    using Base<Scalar, ThisType>::thermalConductivity;
    //! \copybrief Base<Scalar,ThisType>::thermalConductivity(const FluidState&,int)
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

    using Base<Scalar, ThisType>::heatCapacity;
    //! \copybrief Base<Scalar,ThisType>::heatCapacity(const FluidState&,int)
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
