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
 * \brief @copybrief Dumux::FluidSystems::TwoPImmiscible
 */
#ifndef DUMUX_2P_IMMISCIBLE_FLUID_SYSTEM_HH
#define DUMUX_2P_IMMISCIBLE_FLUID_SYSTEM_HH

#include <limits>
#include <cassert>

#include <dune/common/exceptions.hh>

#include <dumux/material/fluidsystems/liquidphase.hh>
#include <dumux/material/fluidsystems/gasphase.hh>
#include <dumux/material/fluidstates/immiscible.hh>

#include "base.hh"

namespace Dumux {
namespace FluidSystems {

/*!
 * \ingroup Fluidsystems
 * \brief A fluid system for two-phase models assuming immiscibility and
 *        thermodynamic equilibrium
 *
 * The fluid phases are completely specified by means of their
 * constituting components.
 * The wetting and the non-wetting phase can be defined individually
 * via FluidSystem::LiquidPhase<Component> and
 * FluidSystem::GasPhase<Component>. These phases consist of one pure
 * component.
 * \tparam Scalar the scalar type
 * \tparam WettingPhase the wetting phase fluid system (use FluidSystem::LiquidPhase<Component> / FluidSystem::GasPhase<Component>)
 * \tparam NonwettingPhase the wetting phase fluid system (use FluidSystem::LiquidPhase<Component> / FluidSystem::GasPhase<Component>)
 */
template <class Scalar, class WettingPhase, class NonwettingPhase>
class TwoPImmiscible
: public BaseFluidSystem<Scalar, TwoPImmiscible<Scalar, WettingPhase, NonwettingPhase> >
{
    using ThisType = TwoPImmiscible<Scalar, WettingPhase, NonwettingPhase>;
    using Base = BaseFluidSystem<Scalar, ThisType>;
public:
    /****************************************
     * Fluid phase related static parameters
     ****************************************/

    //! Number of phases in the fluid system
    static constexpr int numPhases = 2;

    //! Index of the wetting phase
    static constexpr int wPhaseIdx = 0;
    //! Index of the non-wetting phase
    static constexpr int nPhaseIdx = 1;

    /*!
     * \brief Return the human readable name of a fluid phase
     * \param phaseIdx The index of the fluid phase to consider
     */
    static std::string phaseName(int phaseIdx)
    {
        assert(0 <= phaseIdx && phaseIdx < numPhases);

        static std::string name[] = {
            std::string("w"),
            std::string("n")
        };
        return name[phaseIdx];
    }

    /*!
     * \brief Return whether a phase is liquid
     * \param phaseIdx The index of the fluid phase to consider
     */
    static bool isLiquid(int phaseIdx)
    {
        assert(0 <= phaseIdx && phaseIdx < numPhases);

        if (phaseIdx == wPhaseIdx)
            return WettingPhase::isLiquid();
        return NonwettingPhase::isLiquid();
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
        if (phaseIdx == wPhaseIdx)
            return WettingPhase::isCompressible();
        return NonwettingPhase::isCompressible();
    }

    /*!
     * \brief Returns true if and only if a fluid phase is assumed to
     *        be an ideal gas.
     *
     * \param phaseIdx The index of the fluid phase to consider
     */
    static bool isIdealGas(int phaseIdx)
    {
        assert(0 <= phaseIdx && phaseIdx < numPhases);

        // let the fluids decide
        if (phaseIdx == wPhaseIdx)
            return WettingPhase::isIdealGas();
        return NonwettingPhase::isIdealGas();
    }

    /****************************************
     * Component related static parameters
     ****************************************/

    //! Number of components in the fluid system
    static constexpr int numComponents = 2;

    //! Index of the wetting phase's component
    static constexpr int wCompIdx = 0;
    //! Index of the non-wetting phase's component
    static constexpr int nCompIdx = 1;

    /*!
     * \brief Return the human readable name of a component
     *
     * \param compIdx index of the component
     */
    static std::string componentName(int compIdx)
    {
        assert(0 <= compIdx && compIdx < numComponents);

        if (compIdx == wCompIdx)
            return WettingPhase::name();
        return NonwettingPhase::name();
    }

    /*!
     * \brief Return the molar mass of a component in \f$\mathrm{[kg/mol]}\f$.
     * \param compIdx index of the component
     */
    static Scalar molarMass(int compIdx)
    {
        assert(0 <= compIdx && compIdx < numComponents);

        if (compIdx == wCompIdx)
            return WettingPhase::molarMass();
        return NonwettingPhase::molarMass();
    }

    /*!
     * \brief Critical temperature of a component \f$\mathrm{[K]}\f$.
     * \param compIdx index of the component
     */
    static Scalar criticalTemperature(int compIdx)
    {
        assert(0 <= compIdx && compIdx < numComponents);

        if (compIdx == wCompIdx)
            return WettingPhase::criticalTemperature();
        return NonwettingPhase::criticalTemperature();
    }

    /*!
     * \brief Critical pressure of a component \f$\mathrm{[Pa]}\f$.
     * \param compIdx index of the component
     */
    static Scalar criticalPressure(int compIdx)
    {
        assert(0 <= compIdx && compIdx < numComponents);

        if (compIdx == wCompIdx)
            return WettingPhase::criticalPressure();
        return NonwettingPhase::criticalPressure();
    }

    /*!
     * \brief The acentric factor of a component \f$\mathrm{[-]}\f$.
     * \param compIdx index of the component
     */
    static Scalar acentricFactor(int compIdx)
    {
        assert(0 <= compIdx && compIdx < numComponents);

        if (compIdx == wCompIdx)
            return WettingPhase::acentricFactor();
        return NonwettingPhase::acentricFactor();
    }

    /****************************************
     * thermodynamic relations
     ****************************************/

    /*!
     * \brief Initialize the fluid system's static parameters
     */
    static void init()
    {
        // two gaseous phases at once do not make sense physically!
        // (But two liquids are fine)
        assert(WettingPhase::isLiquid() || NonwettingPhase::isLiquid());
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
        if (phaseIdx == wPhaseIdx)
            return WettingPhase::density(temperature, pressure);
        return NonwettingPhase::density(temperature, pressure);
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
        if (phaseIdx == wPhaseIdx)
            return WettingPhase::viscosity(temperature, pressure);
        return NonwettingPhase::viscosity(temperature, pressure);
    }

    using Base::fugacityCoefficient;
    /*!
     * \brief Calculate the fugacity coefficient \f$\mathrm{[-]}\f$ of an individual
     *        component in a fluid phase
     *
     * The fugacity coefficient \f$\mathrm{\phi_\kappa_\alpha}\f$ is connected to the
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
        if (phaseIdx == wPhaseIdx)
            return WettingPhase::enthalpy(temperature, pressure);
        return NonwettingPhase::enthalpy(temperature, pressure);
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
        if (phaseIdx == wPhaseIdx)
            return WettingPhase::thermalConductivity(temperature, pressure);
        return NonwettingPhase::thermalConductivity(temperature, pressure);
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
        if (phaseIdx == wPhaseIdx)
            return WettingPhase::heatCapacity(temperature, pressure);
        return NonwettingPhase::heatCapacity(temperature, pressure);
    }
};

} // end namespace FluidSystems
} // end namespace Dumux

#endif
