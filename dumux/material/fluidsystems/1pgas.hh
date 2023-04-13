// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup FluidSystems
 * \brief @copybrief Dumux::FluidSystems::OnePGas
 */
#ifndef DUMUX_GAS_PHASE_HH
#define DUMUX_GAS_PHASE_HH

#include <cassert>
#include <limits>

#include <dune/common/exceptions.hh>

#include <dumux/material/fluidsystems/base.hh>
#include <dumux/material/components/componenttraits.hh>
#include <dumux/io/name.hh>

namespace Dumux {
namespace FluidSystems {

/*!
 * \ingroup FluidSystems
 * \brief A gaseous phase consisting of a single component
 */
template <class Scalar, class ComponentT>
class OnePGas
: public Base<Scalar, OnePGas<Scalar, ComponentT> >
{
    using ThisType = OnePGas<Scalar, ComponentT>;

    static_assert(ComponentTraits<ComponentT>::hasGasState, "The component does not implement a gas state!");

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
    { return IOName::gaseousPhase(); }

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
     * \brief Returns whether the fluid is gaseous
     */
    static constexpr bool isGas(int phaseIdx = 0)
    { return true; }

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
    { return Component::gasIsCompressible(); }

    /*!
     * \brief Returns true if the fluid is assumed to be an ideal gas
     */
    static constexpr bool isIdealGas(int phaseIdx = 0)
    { return Component::gasIsIdeal(); }

    /*!
     * \brief Returns true if and only if a fluid phase is assumed to
     *        have a constant viscosity.
     *
     * \param phaseIdx The index of the fluid phase to consider
     */
    static constexpr bool viscosityIsConstant(int phaseIdx)
    { return Component::gasViscosityIsConstant(); }

    /*!
     * \brief The mass in \f$\mathrm{[kg]}\f$ of one mole of the component.
     */
    static Scalar molarMass(int compIdx = 0)
    {  return Component::molarMass(); }

    /*!
     * \brief Returns the critical temperature in \f$\mathrm{[K]}\f$ of the component
     */
    static Scalar criticalTemperature(int compIdx = 0)
    {  return Component::criticalTemperature(); }

    /*!
     * \brief Returns the critical pressure in \f$\mathrm{[Pa]}\f$ of the component
     */
    static Scalar criticalPressure(int compIdx = 0)
    {  return Component::criticalPressure(); }

    /*!
     * \brief Returns the temperature in \f$\mathrm{[K]}\f$ at the component's triple point.
     */
    static Scalar tripleTemperature(int compIdx = 0)
    {  return Component::tripleTemperature(); }

    /*!
     * \brief Returns the pressure in \f$\mathrm{[Pa]}\f$ at the component's triple point.
     */
    static Scalar triplePressure(int compIdx = 0)
    { return Component::triplePressure(); }

    /*!
     * \brief The vapor pressure in \f$\mathrm{[Pa]}\f$ of the component at a given
     *        temperature.
     * \param T temperature \f$\mathrm{[K]}\f$
     */
    static Scalar vaporPressure(Scalar T)
    { return Component::vaporPressure(T); }

    /*!
     * \brief The density \f$\mathrm{[kg/m^3]}\f$ of the component at a given pressure and temperature.
     * \param temperature The given temperature \f$\mathrm{[K]}\f$
     * \param pressure The given pressure \f$\mathrm{[Pa]}\f$
     */
    static Scalar density(Scalar temperature, Scalar pressure)
    {  return Component::gasDensity(temperature, pressure); }

    using Base<Scalar, ThisType>::density;
    //! \copydoc Base<Scalar,ThisType>::density(const FluidState&,int)
    template <class FluidState>
    static Scalar density(const FluidState &fluidState,
                          const int phaseIdx)
    {
        return density(fluidState.temperature(phaseIdx),
                       fluidState.pressure(phaseIdx));
    }

    /*!
     * \brief The molar density \f$\rho_{mol,\alpha}\f$
     *   of a fluid phase \f$\alpha\f$ in \f$\mathrm{[mol/m^3]}\f$
     *
     * The molar density is defined by the
     * mass density \f$\rho_\alpha\f$ and the component molar mass \f$M_\alpha\f$:
     *
     * \f[\rho_{mol,\alpha} = \frac{\rho_\alpha}{M_\alpha} \;.\f]
     * \param temperature The temperature at which to evaluate the molar density
     * \param pressure The pressure at which to evaluate the molar density
     */
    static Scalar molarDensity(Scalar temperature, Scalar pressure)
    { return Component::gasMolarDensity(temperature, pressure); }

    using Base<Scalar, ThisType>::molarDensity;
    //! \copydoc Base<Scalar,ThisType>::molarDensity(const FluidState&,int)
    template <class FluidState>
    static Scalar molarDensity(const FluidState &fluidState,
                               const int phaseIdx)
    {
        return molarDensity(fluidState.temperature(phaseIdx),
                            fluidState.pressure(phaseIdx));
    }

    /*!
     * \brief The pressure \f$\mathrm{[Pa]}\f$ of the component at a given density and temperature.
     * \param temperature The given temperature \f$\mathrm{[K]}\f$
     * \param density The given density \f$\mathrm{[kg/m^3]}\f$
     */
    static Scalar pressure(Scalar temperature, Scalar density)
    {  return Component::gasPressure(temperature, density); }

    /*!
     * \brief Specific enthalpy \f$\mathrm{[J/kg]}\f$ of the pure component as a gas.
     * \param temperature The given temperature \f$\mathrm{[K]}\f$
     * \param pressure The given pressure \f$\mathrm{[Pa]}\f$
     */
    static const Scalar enthalpy(Scalar temperature, Scalar pressure)
    {  return Component::gasEnthalpy(temperature, pressure); }

    using Base<Scalar, ThisType>::enthalpy;
    //! \copydoc Base<Scalar,ThisType>::enthalpy(const FluidState&,int)
    template <class FluidState>
    static Scalar enthalpy(const FluidState &fluidState,
                           const int phaseIdx)
    {
        return enthalpy(fluidState.temperature(phaseIdx),
                        fluidState.pressure(phaseIdx));
    }

    /*!
     * \brief Specific internal energy \f$\mathrm{[J/kg]}\f$ of the pure component as a gas.
     * \param temperature The given temperature \f$\mathrm{[K]}\f$
     * \param pressure The given pressure \f$\mathrm{[Pa]}\f$
     */
    static const Scalar internalEnergy(Scalar temperature, Scalar pressure)
    { return Component::gasInternalEnergy(temperature, pressure); }

    /*!
     * \brief The dynamic viscosity \f$\mathrm{[Pa s]}\f$ of the pure component at a given pressure and temperature.
     * \param temperature The given temperature \f$\mathrm{[K]}\f$
     * \param pressure The given pressure \f$\mathrm{[Pa]}\f$
     */
    static Scalar viscosity(Scalar temperature, Scalar pressure)
    {  return Component::gasViscosity(temperature, pressure); }

    using Base<Scalar, ThisType>::viscosity;
    //! \copydoc Base<Scalar,ThisType>::viscosity(const FluidState&,int)
    template <class FluidState>
    static Scalar viscosity(const FluidState &fluidState,
                            const int phaseIdx)
    {
        return viscosity(fluidState.temperature(phaseIdx),
                         fluidState.pressure(phaseIdx));
    }

    using Base<Scalar, ThisType>::fugacityCoefficient;
    //! \copydoc Base<Scalar,ThisType>::fugacityCoefficient(const FluidState&,int,int)
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
    //! \copydoc Base<Scalar,ThisType>::diffusionCoefficient(const FluidState&,int,int)
    template <class FluidState>
    static Scalar diffusionCoefficient(const FluidState &fluidState,
                                       int phaseIdx,
                                       int compIdx)
    {
        DUNE_THROW(Dune::InvalidStateException, "Not applicable: Diffusion coefficients");
    }

    using Base<Scalar, ThisType>::binaryDiffusionCoefficient;
    //! \copydoc Base<Scalar,ThisType>::binaryDiffusionCoefficient(const FluidState&,int,int,int)
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
     * \param temperature The given temperature \f$\mathrm{[K]}\f$
     * \param pressure The given pressure \f$\mathrm{[Pa]}\f$
     */
    static Scalar thermalConductivity(Scalar temperature, Scalar pressure)
    { return Component::gasThermalConductivity(temperature, pressure); }

    using Base<Scalar, ThisType>::thermalConductivity;
    //! \copydoc Base<Scalar,ThisType>::thermalConductivity(const FluidState&,int)
    template <class FluidState>
    static Scalar thermalConductivity(const FluidState &fluidState,
                                      const int phaseIdx)
    {
        return thermalConductivity(fluidState.temperature(phaseIdx),
                                   fluidState.pressure(phaseIdx));
    }

    /*!
     * \brief Specific isobaric heat capacity of the fluid \f$\mathrm{[J/(kg K)]}\f$.
     * \param temperature The given temperature \f$\mathrm{[K]}\f$
     * \param pressure The given pressure \f$\mathrm{[Pa]}\f$
     */
    static Scalar heatCapacity(Scalar temperature, Scalar pressure)
    { return Component::gasHeatCapacity(temperature, pressure); }

    using Base<Scalar, ThisType>::heatCapacity;
    //! \copydoc Base<Scalar,ThisType>::heatCapacity(const FluidState&,int)
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
