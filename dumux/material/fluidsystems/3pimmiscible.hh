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
 * \brief @copybrief Dumux::FluidSystems::ThreePImmiscible
 */
#ifndef DUMUX_3P_IMMISCIBLE_FLUID_SYSTEM_HH
#define DUMUX_3P_IMMISCIBLE_FLUID_SYSTEM_HH

#include <cassert>
#include <limits>
#include <iostream>

#include <dune/common/exceptions.hh>

#include <dumux/material/fluidsystems/1pliquid.hh>
#include <dumux/material/fluidsystems/1pgas.hh>
#include <dumux/material/fluidstates/immiscible.hh>
#include <dumux/material/components/base.hh>
#include <dumux/io/name.hh>

namespace Dumux {
namespace FluidSystems {

/*!
 * \ingroup Fluidsystems
 * \brief A fluid system for three-phase models assuming immiscibility and
 *        thermodynamic equilibrium
 *
 * The fluid phases are completely specified by means of their
 * constituting components.
 * The wetting and the nonwetting phase can be defined individually
 * via FluidSystem::OnePLiquid<Scalar, Component>. The gas phase can be defined via
 * FluidSystems::OnePGas<Scalar, Component>
 * These phases consist of one pure component.
 *
 * \tparam Scalar the scalar type
 * \tparam WettingFluid the wetting phase fluid system (use FluidSystem::OnePLiquid<Scalar, Component>)
 * \tparam NonwettingFluid the wetting phase fluid system (use FluidSystem::OnePLiquid<Scalar, Component>)
 * \tparam Gas the gas phase fluid system (use FluidSystem::OnePGas<Scalar, Component>)
 */
template <class Scalar, class WettingFluid, class NonwettingFluid, class Gas>
class ThreePImmiscible
: public Base<Scalar, ThreePImmiscible<Scalar, WettingFluid, NonwettingFluid, Gas> >
{
    static_assert((WettingFluid::numPhases == 1), "WettingFluid has more than one phase");
    static_assert((NonwettingFluid::numPhases == 1), "NonwettingFluid has more than one phase");
    static_assert((Gas::numPhases == 1), "Gas has more than one phase");
    static_assert((WettingFluid::numComponents == 1), "WettingFluid has more than one component");
    static_assert((NonwettingFluid::numComponents == 1), "NonwettingFluid has more than one component");
    static_assert((Gas::numComponents == 1), "Gas has more than one component");

    using ThisType = ThreePImmiscible<Scalar, WettingFluid, NonwettingFluid, Gas>;
    using Base = Dumux::FluidSystems::Base<Scalar, ThisType>;
public:
    /****************************************
     * Fluid phase related static parameters
     ****************************************/

    //! Number of phases in the fluid system
    static constexpr int numPhases = 3;

    //! Index of the wetting phase
    static constexpr int wPhaseIdx = 0;
    //! Index of the nonwetting phase
    static constexpr int nPhaseIdx = 1;
    //! Index of the gas phase
    static constexpr int gPhaseIdx = 2;

    /*!
     * \brief Return the human readable name of a fluid phase
     * \param phaseIdx The index of the fluid phase to consider
     */
    static std::string phaseName(int phaseIdx)
    {
        assert(0 <= phaseIdx && phaseIdx < numPhases);
        switch (phaseIdx)
        {
            case wPhaseIdx: return Components::IsAqueous<typename WettingFluid::Component>::value
                            ? IOName::aqueousPhase() : IOName::naplPhase();
            case nPhaseIdx: return Components::IsAqueous<typename NonwettingFluid::Component>::value
                            ? IOName::aqueousPhase() : IOName::naplPhase();
            case gPhaseIdx: return IOName::gaseousPhase();
        }
        DUNE_THROW(Dune::InvalidStateException, "Invalid phase index " << phaseIdx);
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

        switch (phaseIdx)
        {
            case wPhaseIdx: return WettingFluid::isGas(); break;
            case nPhaseIdx: return NonwettingFluid::isGas(); break;
            case gPhaseIdx: return Gas::isGas(); break;
            default: return false; // TODO: constexpr-compatible throw
        }
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
    static constexpr bool isIdealMixture(int phaseIdx)
    { return true; }

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
        switch(phaseIdx)
        {
            case wPhaseIdx: return WettingFluid::isCompressible(); break;
            case nPhaseIdx: return NonwettingFluid::isCompressible(); break;
            case gPhaseIdx: return Gas::isCompressible(); break;
            default: return false; // TODO: constexpr-compatible throw
        }
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
        switch(phaseIdx)
        {
            case wPhaseIdx: return WettingFluid::isIdealGas(); break;
            case nPhaseIdx: return NonwettingFluid::isIdealGas(); break;
            case gPhaseIdx: return Gas::isIdealGas(); break;
            default: return false; // TODO: constexpr-compatible throw
        }
    }

    /****************************************
     * Component related static parameters
     ****************************************/

    //! Number of components in the fluid system
    static constexpr int numComponents = 3;//WettingFluid::numComponents + NonwettingFluid::numComponents + Gas::numComponents; // TODO: 3??

    //! Index of the wetting phase's component
    static constexpr int wCompIdx = 0;
    //! Index of the nonwetting phase's component
    static constexpr int nCompIdx = 1;
    //! Index of the gas phase's component
    static constexpr int gCompIdx = 2; // TODO: correct??

    /*!
     * \brief Return the human readable name of a component
     *
     * \param compIdx index of the component
     */
    static std::string componentName(int compIdx)
    {
        assert(0 <= compIdx && compIdx < numComponents);

        switch(compIdx)
        {
            case wCompIdx: return WettingFluid::name(); break;
            case nCompIdx: return NonwettingFluid::name(); break;
            case gCompIdx: return Gas::name(); break;
            default: DUNE_THROW(Dune::InvalidStateException, "Invalid component index");
        }
    }

    /*!
     * \brief Return the molar mass of a component in \f$\mathrm{[kg/mol]}\f$.
     * \param compIdx index of the component
     */
    static Scalar molarMass(int compIdx)
    {
        assert(0 <= compIdx && compIdx < numComponents);

        switch(compIdx)
        {
            case wCompIdx: return WettingFluid::molarMass(); break;
            case nCompIdx: return NonwettingFluid::molarMass(); break;
            case gCompIdx: return Gas::molarMass(); break;
            default: DUNE_THROW(Dune::InvalidStateException, "Invalid component index");
        }
    }

    /*!
     * \brief Critical temperature of a component \f$\mathrm{[K]}\f$.
     * \param compIdx index of the component
     */
    static Scalar criticalTemperature(int compIdx)
    {
        assert(0 <= compIdx && compIdx < numComponents);

        switch(compIdx)
        {
            case wCompIdx: return WettingFluid::criticalTemperature(); break;
            case nCompIdx: return NonwettingFluid::criticalTemperature(); break;
            case gCompIdx: return Gas::criticalTemperature(); break;
            default: DUNE_THROW(Dune::InvalidStateException, "Invalid component index");
        }
    }

    /*!
     * \brief Critical pressure of a component \f$\mathrm{[Pa]}\f$.
     * \param compIdx index of the component
     */
    static Scalar criticalPressure(int compIdx)
    {
        assert(0 <= compIdx && compIdx < numComponents);

        switch(compIdx)
        {
            case wCompIdx: return WettingFluid::criticalPressure(); break;
            case nCompIdx: return NonwettingFluid::criticalPressure(); break;
            case gCompIdx: return Gas::criticalPressure(); break;
            default: DUNE_THROW(Dune::InvalidStateException, "Invalid component index");
        }
    }

    /*!
     * \brief The acentric factor of a component \f$\mathrm{[-]}\f$.
     * \param compIdx index of the component
     */
    static Scalar acentricFactor(int compIdx)
    {
        assert(0 <= compIdx && compIdx < numComponents);

        switch(compIdx)
        {
            case wCompIdx: return WettingFluid::Component::acentricFactor(); break;
            case nCompIdx: return NonwettingFluid::Component::acentricFactor(); break;
            case gCompIdx: return Gas::Component::acentricFactor(); break;
            default: DUNE_THROW(Dune::InvalidStateException, "Invalid component index");
        }
    }

    /****************************************
     * thermodynamic relations
     ****************************************/

    /*!
     * \brief Initialize the fluid system's static parameters
     */
    static constexpr void init()
    {
        // two gaseous phases at once do not make sense physically!
        // (But two liquids are fine)
        static_assert(!WettingFluid::isGas() && !NonwettingFluid::isGas() && Gas::isGas(), "There can only be one gaseous phase!");
    }

    /*!
     * \brief Initialize the fluid system's static parameters using
     *        problem specific temperature and pressure ranges
     *
     * \param tempMin The minimum temperature used for tabulation of water \f$\mathrm{[K]}\f$
     * \param tempMax The maximum temperature used for tabulation of water \f$\mathrm{[K]}\f$
     * \param nTemp The number of ticks on the temperature axis of the  table of water
     * \param pressMin The minimum pressure used for tabulation of water \f$\mathrm{[Pa]}\f$
     * \param pressMax The maximum pressure used for tabulation of water \f$\mathrm{[Pa]}\f$
     * \param nPress The number of ticks on the pressure axis of the  table of water
     */
    static void init(Scalar tempMin, Scalar tempMax, unsigned nTemp,
                     Scalar pressMin, Scalar pressMax, unsigned nPress)
    {
        // two gaseous phases at once do not make sense physically!
        static_assert(!WettingFluid::isGas() && !NonwettingFluid::isGas() && Gas::isGas(), "There can only be one gaseous phase!");

        if (WettingFluid::Component::isTabulated)
        {
            std::cout << "Initializing tables for the wetting fluid properties ("
                      << nTemp*nPress
                      << " entries).\n";

            WettingFluid::Component::init(tempMin, tempMax, nTemp,
                                          pressMin, pressMax, nPress);

        }

        if (NonwettingFluid::Component::isTabulated)
        {
            std::cout << "Initializing tables for the nonwetting fluid properties ("
                      << nTemp*nPress
                      << " entries).\n";

            NonwettingFluid::Component::init(tempMin, tempMax, nTemp,
                                             pressMin, pressMax, nPress);

        }

        if (Gas::Component::isTabulated)
        {
            std::cout << "Initializing tables for the gas fluid properties ("
                      << nTemp*nPress
                      << " entries).\n";

            Gas::Component::init(tempMin, tempMax, nTemp,
                                 pressMin, pressMax, nPress);

        }
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

        switch(phaseIdx)
        {
            case wPhaseIdx: return WettingFluid::density(temperature, pressure); break;
            case nPhaseIdx: return NonwettingFluid::density(temperature, pressure); break;
            case gPhaseIdx: return Gas::density(temperature, pressure); break;
            default: DUNE_THROW(Dune::InvalidStateException, "Invalid phase index");
        }
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
    static Scalar molarDensity(const FluidState &fluidState,
                               int phaseIdx)
    {
        assert(0 <= phaseIdx && phaseIdx < numPhases);

        Scalar temperature = fluidState.temperature(phaseIdx);
        Scalar pressure = fluidState.pressure(phaseIdx);

        switch(phaseIdx)
        {
            case wPhaseIdx: return WettingFluid::molarDensity(temperature, pressure);
            case nPhaseIdx: return NonwettingFluid::molarDensity(temperature, pressure);
            case gPhaseIdx: return Gas::molarDensity(temperature, pressure);
            default: DUNE_THROW(Dune::InvalidStateException, "Invalid phase index");
        }
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

        switch(phaseIdx)
        {
            case wPhaseIdx: return WettingFluid::viscosity(temperature, pressure); break;
            case nPhaseIdx: return NonwettingFluid::viscosity(temperature, pressure); break;
            case gPhaseIdx: return Gas::viscosity(temperature, pressure); break;
            default: DUNE_THROW(Dune::InvalidStateException, "Invalid phase index");
        }
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
        DUNE_THROW(Dune::InvalidStateException,
                   "Binary diffusion coefficients of components are meaningless if"
                   " immiscibility is assumed");
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

        switch(phaseIdx)
        {
            case wPhaseIdx: return WettingFluid::enthalpy(temperature, pressure); break;
            case nPhaseIdx: return NonwettingFluid::enthalpy(temperature, pressure); break;
            case gPhaseIdx: return Gas::enthalpy(temperature, pressure); break;
            default: DUNE_THROW(Dune::InvalidStateException, "Invalid phase index");
        }
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

        switch(phaseIdx)
        {
            case wPhaseIdx: return WettingFluid::thermalConductivity(temperature, pressure); break;
            case nPhaseIdx: return NonwettingFluid::thermalConductivity(temperature, pressure); break;
            case gPhaseIdx: return Gas::thermalConductivity(temperature, pressure); break;
            default: DUNE_THROW(Dune::InvalidStateException, "Invalid phase index");
        }
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

        switch(phaseIdx)
        {
            case wPhaseIdx: return WettingFluid::heatCapacity(temperature, pressure); break;
            case nPhaseIdx: return NonwettingFluid::heatCapacity(temperature, pressure); break;
            case gPhaseIdx: return Gas::heatCapacity(temperature, pressure); break;
            default: DUNE_THROW(Dune::InvalidStateException, "Invalid phase index");
        }
    }
};

} // end namespace FluidSystems
} // end namespace Dumux

#endif
