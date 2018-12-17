// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   See the file COPYING for full copying permissions.                      *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 3 of the License, or       *
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
 * \brief @copybrief Dumux::FluidSystems::OnePLiquid
 */
#ifndef DUMUX_LIQUID_PHASE_HH
#define DUMUX_LIQUID_PHASE_HH

#include <cassert>
#include <limits>

#include <dune/common/exceptions.hh>

#include <dumux/material/fluidsystems/base.hh>
#include <dumux/material/components/componenttraits.hh>
#include <dumux/io/name.hh>

namespace Dumux {
namespace FluidSystems {

/*!
 * \ingroup Fluidsystems
 * \brief A liquid phase consisting of a single component
 */
template <class Scalar, class ComponentT>
class OnePLiquid
: public Base<Scalar, OnePLiquid<Scalar, ComponentT> >
{
    using ThisType = OnePLiquid<Scalar, ComponentT>;
    using Base = Dumux::FluidSystems::Base<Scalar, ThisType>;

    static_assert(ComponentTraits<ComponentT>::hasLiquidState, "The component does not implement a liquid state!");

public:
    using Component = ComponentT;
    using ParameterCache = NullParameterCache;

    static constexpr int numPhases = 1;  //!< Number of phases in the fluid system
    static constexpr int numComponents = 1; //!< Number of components in the fluid system

    static constexpr int phase0Idx = 0; //!< index of the only phase
    static constexpr int comp0Idx = 0; //!< index of the only component

    /*!
     * \brief Initialize the fluid system's static parameters generically
     */
    static void init()
    { }

    /****************************************
     * Fluid phase related static parameters
     ****************************************/
    /*!
     * \brief Return the human readable name of a fluid phase
     *
     * \param phaseIdx The index of the fluid phase to consider
     */
    static std::string phaseName(int phaseIdx = 0)
    { return IOName::liquidPhase(); }

    /*!
     * \brief A human readable name for the component.
     *
     * \param compIdx The index of the component to consider
     */
    static std::string componentName(int compIdx = 0)
    { return Component::name(); }

    /*!
     * \brief A human readable name for the component.
     */
    static std::string name()
    { return Component::name(); }

    /*!
     * \brief There is only one phase, so not mass transfer between phases can occur
     */
    static constexpr bool isMiscible()
    { return false; }

    /*!
     * \brief Returns whether the fluid is a liquid
     */
    static constexpr bool isGas(int phaseIdx = 0)
    { return false; }

    /*!
     * \brief Returns true if and only if a fluid phase is assumed to
     *        be an ideal mixture.
     *
     * We define an ideal mixture as a fluid phase where the fugacity
     * coefficients of all components times the pressure of the phase
     * are independent on the fluid composition. This assumption is true
     * if only a single component is involved. If you are unsure what
     * this function should return, it is safe to return false. The
     * only damage done will be (slightly) increased computation times
     * in some cases.
     *
     * \param phaseIdx The index of the fluid phase to consider
     */
    static constexpr bool isIdealMixture(int phaseIdx = 0)
    { return true; }

    /*!
     * \brief Returns true if the fluid is assumed to be compressible
     */
    static constexpr bool isCompressible(int phaseIdx = 0)
    { return Component::liquidIsCompressible(); }

    /*!
     * \brief Returns true if the fluid viscosity is constant
     */
    static constexpr bool viscosityIsConstant(int phaseIdx = 0)
    { return Component::liquidViscosityIsConstant(); }

    /*!
     * \brief Returns true if the fluid is assumed to be an ideal gas
     */
    static constexpr bool isIdealGas(int phaseIdx = 0)
    { return false; /* we're a liquid! */ }

    /*!
     * \brief The mass in \f$\mathrm{[kg]}\f$ of one mole of the component.
     */
    static Scalar molarMass(int compIdx = 0)
    {  return Component::molarMass(); }

    /*!
     * \brief Returns the critical temperature \f$\mathrm{[K]}\f$ of the component
     */
    static Scalar criticalTemperature(int compIdx = 0)
    {  return Component::criticalTemperature(); }

    /*!
     * \brief Returns the critical pressure \f$\mathrm{[Pa]}\f$ of the component
     */
    static Scalar criticalPressure(int compIdx = 0)
    {  return Component::criticalPressure(); }

    /*!
     * \brief Returns the temperature \f$\mathrm{[K]}\f$ at the component's triple point.
     */
    static Scalar tripleTemperature(int compIdx = 0)
    {  return Component::tripleTemperature(); }

    /*!
     * \brief Returns the pressure \f$\mathrm{[Pa]}\f$ at the component's triple point.
     */
    static Scalar triplePressure(int compIdx = 0)
    { return Component::triplePressure(); }

    /*!
     * \brief The vapor pressure in \f$\mathrm{[Pa]}\f$ of the component at a given
     *        temperature.
     */
    static Scalar vaporPressure(Scalar T)
    {  return Component::vaporPressure(T); }

    /*!
     * \brief The density \f$\mathrm{[kg/m^3]}\f$ of the component at a given pressure and temperature.
     */
    static Scalar density(Scalar temperature, Scalar pressure)
    {  return Component::liquidDensity(temperature, pressure); }

    using Base::density;
    /*!
     * \brief The density \f$\mathrm{[kg/m^3]}\f$ of the component at a given pressure and temperature.
     */
    template <class FluidState>
    static Scalar density(const FluidState &fluidState,
                          const int phaseIdx)
    {
        return density(fluidState.temperature(phaseIdx),
                       fluidState.pressure(phaseIdx));
    }

    using Base::molarDensity;
    /*!
     * \brief The molar density \f$\rho_{mol,\alpha}\f$
     *   of a fluid phase \f$\alpha\f$ in \f$\mathrm{[mol/m^3]}\f$
     *
     * The molar density is defined by the
     * mass density \f$\rho_\alpha\f$ and the main component molar mass \f$M_\alpha\f$:
     *
     * \f[\rho_{mol,\alpha} = \frac{\rho_\alpha}{M_\alpha} \;.\f]
     */
    template <class FluidState>
    static Scalar molarDensity(const FluidState &fluidState, const int phaseIdx)
    {
        return molarDensity(fluidState.temperature(phaseIdx),
                            fluidState.pressure(phaseIdx));
    }

    /*!
     * \brief The density \f$\mathrm{[kg/m^3]}\f$ of the component at a given pressure and temperature.
     */
    static Scalar molarDensity(Scalar temperature, Scalar pressure)
    {  return Component::liquidMolarDensity(temperature, pressure); }

    /*!
     * \brief The pressure \f$\mathrm{[Pa]}\f$ of the component at a given density and temperature.
     */
    static Scalar pressure(Scalar temperature, Scalar density)
    {  return Component::liquidPressure(temperature, density); }

    /*!
     * \brief Specific enthalpy \f$\mathrm{[J/kg]}\f$ the pure component as a liquid.
     */
    static const Scalar enthalpy(Scalar temperature, Scalar pressure)
    {  return Component::liquidEnthalpy(temperature, pressure); }

    using Base::enthalpy;
    /*!
     * \brief Specific enthalpy \f$\mathrm{[J/kg]}\f$ the pure component as a liquid.
     */
    template <class FluidState>
    static Scalar enthalpy(const FluidState &fluidState,
                           const int phaseIdx)
    {
        return enthalpy(fluidState.temperature(phaseIdx),
                        fluidState.pressure(phaseIdx));
    }

    /*!
     * \brief Specific internal energy \f$\mathrm{[J/kg]}\f$ the pure component as a liquid.
     */
    static const Scalar internalEnergy(Scalar temperature, Scalar pressure)
    { return Component::liquidInternalEnergy(temperature, pressure); }

    /*!
     * \brief The dynamic liquid viscosity \f$\mathrm{[N/m^3*s]}\f$ of the pure component.
     */
    static Scalar viscosity(Scalar temperature, Scalar pressure)
    {  return Component::liquidViscosity(temperature, pressure); }

    using Base::viscosity;
    /*!
     * \brief The dynamic liquid viscosity \f$\mathrm{[N/m^3*s]}\f$ of the pure component.
     */
    template <class FluidState>
    static Scalar viscosity(const FluidState &fluidState,
                            const int phaseIdx)
    {
        return viscosity(fluidState.temperature(phaseIdx),
                         fluidState.pressure(phaseIdx));
    }

    using Base::fugacityCoefficient;
    /*!
     * \copybrief Base::fugacityCoefficient
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

        if (phaseIdx == compIdx)
            // We could calculate the real fugacity coefficient of
            // the component in the fluid. Probably that's not worth
            // the effort, since the fugacity coefficient of the other
            // component is infinite anyway...
            return 1.0;
        return std::numeric_limits<Scalar>::infinity();
    }

    using Base::diffusionCoefficient;
    /*!
     * \copybrief Base::diffusionCoefficient
     *
     * \param fluidState An arbitrary fluid state
     * \param phaseIdx The index of the fluid phase to consider
     * \param compIdx The index of the component to consider
     */
    template <class FluidState>
    static Scalar diffusionCoefficient(const FluidState &fluidState,
                                       int phaseIdx,
                                       int compIdx)
    {
        DUNE_THROW(Dune::InvalidStateException, "Not applicable: Diffusion coefficients");
    }

    using Base::binaryDiffusionCoefficient;
    /*!
     * \copybrief Base::binaryDiffusionCoefficient
     *
     * \param fluidState An arbitrary fluid state
     * \param phaseIdx The index of the fluid phase to consider
     * \param compIIdx The index of the component to consider
     * \param compJIdx The index of the component to consider
     */
    template <class FluidState>
    static Scalar binaryDiffusionCoefficient(const FluidState &fluidState,
                                             int phaseIdx,
                                             int compIIdx,
                                             int compJIdx)
    {
        DUNE_THROW(Dune::InvalidStateException, "Not applicable: Binary diffusion coefficients");
    }

    /*!
     * \brief Thermal conductivity of the fluid \f$\mathrm{[W/(m K)]}\f$.
     */
    static Scalar thermalConductivity(Scalar temperature, Scalar pressure)
    { return Component::liquidThermalConductivity(temperature, pressure); }

    using Base::thermalConductivity;
    /*!
     * \brief Thermal conductivity of the fluid \f$\mathrm{[W/(m K)]}\f$.
     */
    template <class FluidState>
    static Scalar thermalConductivity(const FluidState &fluidState,
                                      const int phaseIdx)
    {
        return thermalConductivity(fluidState.temperature(phaseIdx),
                                   fluidState.pressure(phaseIdx));
    }

    /*!
     * \brief Specific isobaric heat capacity of the fluid \f$\mathrm{[J/(kg K)]}\f$.
     */
    static Scalar heatCapacity(Scalar temperature, Scalar pressure)
    { return Component::liquidHeatCapacity(temperature, pressure); }

    using Base::heatCapacity;
    /*!
     * \brief Specific isobaric heat capacity of the fluid \f$\mathrm{[J/(kg K)]}\f$.
     */
    template <class FluidState>
    static Scalar heatCapacity(const FluidState &fluidState,
                               const int phaseIdx)
    {
        return heatCapacity(fluidState.temperature(phaseIdx),
                            fluidState.pressure(phaseIdx));
    }
};

} // namespace FluidSystems
} // namespace

#endif
