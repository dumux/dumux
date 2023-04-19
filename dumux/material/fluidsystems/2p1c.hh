// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup FluidSystems
 * \copybrief Dumux::FluidSystems::TwoPOneC
 */
#ifndef DUMUX_2P_1C_FLUID_SYSTEM_HH
#define DUMUX_2P_1C_FLUID_SYSTEM_HH

#include <limits>
#include <cassert>
#include <iostream>

#include <dune/common/exceptions.hh>

#include <dumux/io/name.hh>

#include "base.hh"

namespace Dumux {
namespace FluidSystems {

/*!
 * \ingroup FluidSystems
 * \brief A two-phase fluid system with only one component.
 */
template <class Scalar, class ComponentType>
class TwoPOneC
    : public Base<Scalar, TwoPOneC<Scalar, ComponentType> >
{
    using ThisType = TwoPOneC<Scalar, ComponentType>;
    using Component = ComponentType;

public:
    static constexpr int numPhases = 2; //!< Number of phases in the fluid system
    static constexpr int numComponents = 1; //!< Number of components in the fluid system

    static constexpr int liquidPhaseIdx = 0; //!< index of the liquid phase
    static constexpr int gasPhaseIdx = 1; //!< index of the gas phase
    static constexpr int phase0Idx = liquidPhaseIdx; //!< index of the first phase
    static constexpr int phase1Idx = gasPhaseIdx; //!< index of the second phase

    static constexpr int comp0Idx = 0; //!< index of the only component

    /****************************************
    * Fluid phase related static parameters
    ****************************************/
    /*!
     * \brief Return the human readable name of a fluid phase
     *
     * \param phaseIdx The index of the fluid phase to consider
     */
    static std::string phaseName(int phaseIdx)
    {
        static std::string name[] = {
            std::string(IOName::liquidPhase()),
            std::string(IOName::gaseousPhase()),
        };

        assert(0 <= phaseIdx && phaseIdx < numPhases);
        return name[phaseIdx];
    }

    /*!
     * \brief Returns whether the fluids are miscible
     */
    static constexpr bool isMiscible()
    { return false; }

    /*!
     * \brief Return whether a phase is gaseous
     *
     * \param phaseIdx The index of the fluid phase to consider
     */
    static constexpr bool isGas(int phaseIdx)
    {
        assert(0 <= phaseIdx && phaseIdx < numPhases);
        return phaseIdx == gasPhaseIdx;
    }

    /*!
     * \brief Returns true if and only if a fluid phase is assumed to
     *        be an ideal mixture.
     *
     * We define an ideal mixture as a fluid phase where the fugacity
     * coefficients of all components times the pressure of the phase
     * are independent on the fluid composition. This assumption is true
     * if Henry's law and Raoult's law apply. If you are unsure what
     * this function should return, it is safe to return false. The
     * only damage done will be (slightly) increased computation times
     * in some cases.
     *
     * \param phaseIdx The index of the fluid phase to consider
     */
    static constexpr bool isIdealMixture(int phaseIdx)
    {
        assert(0 <= phaseIdx && phaseIdx < numPhases);
        // we assume Henry's and Raoult's laws for the water phase and
        // and no interaction between gas molecules of different
        // components, so all phases are ideal mixtures!
        return true;
    }

    /*!
     * \brief Returns true if and only if a fluid phase is assumed to
     *        be compressible.
     *
     * Compressible means that the partial derivative of the density
     * to the fluid pressure is always larger than zero.
     *
     * \param phaseIdx The index of the fluid phase to consider
     */
    static constexpr bool isCompressible(int phaseIdx)
    {
        assert(0 <= phaseIdx && phaseIdx < numPhases);
        // gases are always compressible
        if (phaseIdx == gasPhaseIdx)
            return true;
        // the component decides for the liquid phase...
        return Component::liquidIsCompressible();
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

        if (phaseIdx == gasPhaseIdx)
            // let the components decide
            return Component::gasIsIdeal();
        return false; // not a gas
    }

    /****************************************
     * Component related static parameters
     ****************************************/

    /*!
     * \brief Return the human readable name of a component
     *
     * \param compIdx The index of the component to consider
     */
    static std::string componentName(int compIdx)
    {
        assert(0 <= compIdx && compIdx < numComponents);
        return Component::name();
    }

    /*!
     * \brief Return the molar mass of a component in [kg/mol].
     *
     * \param compIdx The index of the component to consider
     */
    static Scalar molarMass(int compIdx)
    {
        assert(0 <= compIdx && compIdx < numComponents);
        return Component::molarMass();
    }

    /*!
     * \brief Critical temperature of a component [K].
     *
     * \param compIdx The index of the component to consider
     */
    static Scalar criticalTemperature(int compIdx)
    {
        assert(0 <= compIdx && compIdx < numComponents);
        return Component::criticalTemperature();
    }

    /*!
     * \brief Critical pressure of a component [Pa].
     *
     * \param compIdx The index of the component to consider
     */
    static Scalar criticalPressure(int compIdx)
    {
        assert(0 <= compIdx && compIdx < numComponents);
        return Component::criticalPressure();
    }

    /*!
     * \brief Molar volume of a component at the critical point [m^3/mol].
     *
     * \param compIdx The index of the component to consider
     */
    static Scalar criticalMolarVolume(int compIdx)
    {
        DUNE_THROW(Dune::NotImplemented,
                   "TwoPOneC::criticalMolarVolume()");
        return 0;
    }

    /*!
     * \brief The acentric factor of a component [].
     *
     * \param compIdx The index of the component to consider
     */
    static Scalar acentricFactor(int compIdx)
    {
        assert(0 <= compIdx && compIdx < numComponents);
        return Component::acentricFactor();
    }

    /****************************************
     * thermodynamic relations
     ****************************************/

    /*!
     * \brief Initialize the fluid system's static parameters generically
     *
     * If a tabulated H2O component is used, we do our best to create
     * tables that always work.
     */
    static void init()
    {
        init(/*tempMin=*/273.15,
             /*tempMax=*/623.15,
             /*numTemp=*/100,
             /*pMin=*/0.0,
             /*pMax=*/20e6,
             /*numP=*/200);
    }

    /*!
     * \brief Initialize the fluid system's static parameters using
     *        problem specific temperature and pressure ranges
     *
     * \param tempMin The minimum temperature used for tabulation of water [K]
     * \param tempMax The maximum temperature used for tabulation of water [K]
     * \param nTemp The number of ticks on the temperature axis of the  table of water
     * \param pressMin The minimum pressure used for tabulation of water [Pa]
     * \param pressMax The maximum pressure used for tabulation of water [Pa]
     * \param nPress The number of ticks on the pressure axis of the  table of water
     */
    static void init(Scalar tempMin, Scalar tempMax, unsigned nTemp,
                     Scalar pressMin, Scalar pressMax, unsigned nPress)
    {
        if (Component::isTabulated)
        {
            Component::init(tempMin, tempMax, nTemp,
                            pressMin, pressMax, nPress);
        }
    }

    using Base<Scalar, ThisType>::density;
    //! \copydoc Base<Scalar,ThisType>::density(const FluidState&,int)
    template <class FluidState>
    static Scalar density(const FluidState &fluidState,
                          int phaseIdx)
    {
        assert(0 <= phaseIdx  && phaseIdx < numPhases);

        Scalar temperature = fluidState.temperature(phaseIdx);
        Scalar pressure = fluidState.pressure(phaseIdx);

        // liquid phase
        if (phaseIdx == liquidPhaseIdx) {
            return Component::liquidDensity(temperature, pressure);
        }
        else if (phaseIdx == gasPhaseIdx)// gas phase
        {
            return Component::gasDensity(temperature, pressure);
        }
        else
            DUNE_THROW(Dune::InvalidStateException, "Invalid phase index.");
    }

    using Base<Scalar, ThisType>::molarDensity;
    //! \copydoc Base<Scalar,ThisType>::molarDensity(const FluidState&,int)
    template <class FluidState>
    static Scalar molarDensity(const FluidState &fluidState, int phaseIdx)
    {
        Scalar temperature = fluidState.temperature(phaseIdx);
        Scalar pressure = fluidState.temperature(phaseIdx);
        if (phaseIdx == liquidPhaseIdx)
            return Component::liquidMolarDensity(temperature, pressure);
        else if (phaseIdx == gasPhaseIdx)
            return Component::gasMolarDensity(temperature, pressure);
        else
            DUNE_THROW(Dune::InvalidStateException, "Invalid phase index.");
    }

    using Base<Scalar, ThisType>::viscosity;
    //! \copydoc Base<Scalar,ThisType>::viscosity(const FluidState&,int)
    template <class FluidState>
    static Scalar viscosity(const FluidState &fluidState,
                            int phaseIdx)
    {
        assert(0 <= phaseIdx  && phaseIdx < numPhases);

        Scalar temperature = fluidState.temperature(phaseIdx);
        Scalar pressure = fluidState.pressure(phaseIdx);

        // liquid phase
        if (phaseIdx == liquidPhaseIdx)
            return Component::liquidViscosity(temperature, pressure);
        else if (phaseIdx == gasPhaseIdx) // gas phase
            return Component::gasViscosity(temperature, pressure);
        else
            DUNE_THROW(Dune::InvalidStateException, "Invalid phase index.");
    }

    /*!
     * \brief calculate the temperature of vapor at a given pressure on the vapor pressure curve.
     *
     * \param fluidState An arbitrary fluid state
     * \param phaseIdx The index of the fluid phase to consider
     */
    template <class FluidState>
    static Scalar vaporTemperature(const FluidState &fluidState,
                                   const unsigned int phaseIdx)
    {
        assert(0 <= phaseIdx  && phaseIdx < numPhases);
        Scalar pressure = fluidState.pressure(gasPhaseIdx);

        return Component::vaporTemperature(pressure);
    }

    using Base<Scalar, ThisType>::fugacityCoefficient;
    /*!
     * \brief Calculate the fugacity coefficient [-] of an individual
     *        component in a fluid phase
     *
     * The fugacity coefficient \f$\phi^\kappa_\alpha\f$ of
     * component \f$\kappa\f$ in phase \f$\alpha\f$ is connected to
     * the fugacity \f$f^\kappa_\alpha\f$ and the component's mole
     * fraction \f$x^\kappa_\alpha\f$ by means of the relation
     *
     * \f[
     f^\kappa_\alpha = \phi^\kappa_\alpha\;x^\kappa_\alpha\;p_\alpha
     \f]
     * where \f$p_\alpha\f$ is the pressure of the fluid phase.
     *
     * The quantity "fugacity" itself is just an other way to express
     * the chemical potential \f$\zeta^\kappa_\alpha\f$ of the
     * component. It is defined via
     *
     * \f[
     f^\kappa_\alpha := \exp\left\{\frac{\zeta^\kappa_\alpha}{k_B T_\alpha} \right\}
     \f]
     * where \f$k_B = 1.380\cdot10^{-23}\;J/K\f$ is the Boltzmann constant.
     *
     * \param fluidState An arbitrary fluid state
     * \param phaseIdx The index of the fluid phase to consider
     * \param compIdx The index of the component to consider
     */
    template <class FluidState>
    static Scalar fugacityCoefficient(const FluidState &fluidState,
                                      int phaseIdx,
                                      int compIdx)
    {
        assert(0 <= phaseIdx  && phaseIdx < numPhases);
        assert(0 <= compIdx  && compIdx < numComponents);

        Scalar temperature = fluidState.temperature(phaseIdx);
        Scalar pressure = fluidState.pressure(phaseIdx);

        // liquid phase
        if (phaseIdx == liquidPhaseIdx)
           return Component::vaporPressure(temperature)/pressure;

        // for the gas phase, assume an ideal gas when it comes to
        // fugacity (-> fugacity == partial pressure)
        else if (phaseIdx == gasPhaseIdx)
            return 1.0;

        else
            DUNE_THROW(Dune::InvalidStateException, "Invalid phase index.");
    }


    using Base<Scalar, ThisType>::diffusionCoefficient;
    //! \copydoc Base<Scalar,ThisType>::diffusionCoefficient(const FluidState&,int,int)
    template <class FluidState>
    static Scalar diffusionCoefficient(const FluidState &fluidState,
                                       int phaseIdx,
                                       int compIdx)
    {
        // TODO!
        DUNE_THROW(Dune::NotImplemented, "Diffusion coefficients");
    }

    using Base<Scalar, ThisType>::binaryDiffusionCoefficient;
    //! \copydoc Base<Scalar,ThisType>::binaryDiffusionCoefficient(const FluidState&,int,int,int)
    template <class FluidState>
    static Scalar binaryDiffusionCoefficient(const FluidState &fluidState,
                                             int phaseIdx,
                                             int compIIdx,
                                             int compJIdx)

    { DUNE_THROW(Dune::NotImplemented, "Binary Diffusion coefficients"); }

    using Base<Scalar, ThisType>::enthalpy;
    //! \copydoc Base<Scalar,ThisType>::enthalpy(const FluidState&,int)
    template <class FluidState>
    static Scalar enthalpy(const FluidState &fluidState,
                           int phaseIdx)
    {
        assert(0 <= phaseIdx  && phaseIdx < numPhases);

        // liquid phase
        if (phaseIdx == liquidPhaseIdx)
            return Component::liquidEnthalpy(fluidState.temperature(phaseIdx),
                                             fluidState.pressure(phaseIdx));

        else if (phaseIdx == gasPhaseIdx) // gas phase
            return Component::gasEnthalpy(fluidState.temperature(phaseIdx),
                                          fluidState.pressure(phaseIdx));

        else
            DUNE_THROW(Dune::InvalidStateException, "Invalid phase index.");
    }

    using Base<Scalar, ThisType>::thermalConductivity;
    //! \copydoc Base<Scalar,ThisType>::thermalConductivity(const FluidState&,int)
    template <class FluidState>
    static Scalar thermalConductivity(const FluidState &fluidState,
                                      const int phaseIdx)
    {
        assert(0 <= phaseIdx  && phaseIdx < numPhases);
        // liquid phase
        if (phaseIdx == liquidPhaseIdx)
            return Component::liquidThermalConductivity(fluidState.temperature(phaseIdx),
                                                        fluidState.pressure(phaseIdx)); //0.68 ;

        else if (phaseIdx == gasPhaseIdx) // gas phase
            return Component::gasThermalConductivity(fluidState.temperature(phaseIdx),
                                                     fluidState.pressure(phaseIdx)); //0.0248;

        else
            DUNE_THROW(Dune::InvalidStateException, "Invalid phase index.");
    }

    using Base<Scalar, ThisType>::heatCapacity;
    //! \copydoc Base<Scalar,ThisType>::heatCapacity(const FluidState&,int)
    template <class FluidState>
    static Scalar heatCapacity(const FluidState &fluidState,
                               int phaseIdx)
    {
        assert(0 <= phaseIdx  && phaseIdx < numPhases);
        // liquid phase
        if (phaseIdx == liquidPhaseIdx)
            return Component::liquidHeatCapacity(fluidState.temperature(phaseIdx),
                                                 fluidState.pressure(phaseIdx));//4.217e3 ;

        else if (phaseIdx == gasPhaseIdx) // gas phase
            return Component::gasHeatCapacity(fluidState.temperature(phaseIdx),
                                              fluidState.pressure(phaseIdx));//2.029e3;

        else
            DUNE_THROW(Dune::InvalidStateException, "Invalid phase index.");
    }
};

} // end namespace FluidSystems
} // end namespace Dumux

#endif
