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
 *
 * \brief @copybrief Dumux::FluidSystems::OneP
 */
#ifndef DUMUX_1P_FLUIDSYSTEM_HH
#define DUMUX_1P_FLUIDSYSTEM_HH

#include <limits>
#include <cassert>

#include <dumux/material/fluidsystems/liquidphase.hh>
#include <dumux/material/fluidsystems/gasphase.hh>

#include <dune/common/exceptions.hh>

#include "base.hh"
#include <dumux/material/components/simpleh2o.hh>
#include <dumux/material/components/h2o.hh>
#include <dumux/material/components/n2.hh>
#include <dumux/material/components/tabulatedcomponent.hh>

#ifdef DUMUX_PROPERTIES_HH
#include <dumux/common/basicproperties.hh>
#include <dumux/material/fluidsystems/defaultcomponents.hh>
#endif

namespace Dumux {
namespace FluidSystems {

/*!
 * \ingroup Fluidsystems
 *
 * \brief A fluid system for single phase models.
 *
 * The fluid is defined as a template parameter. For existing
 * components the FluidSystems::LiquidPhase<Component> and
 * FluidSystems::GasPhase<Component> may be used.
 */
template <class Scalar, class Fluid>
class OneP
    : public BaseFluidSystem<Scalar, OneP<Scalar, Fluid> >
{
    typedef OneP<Scalar, Fluid> ThisType;
    typedef BaseFluidSystem<Scalar, ThisType> Base;

public:
    /****************************************
     * Fluid phase related static parameters
     ****************************************/

    //! Number of phases in the fluid system
    static constexpr int numPhases = 1;

    /*!
     * \brief Return the human readable name of a fluid phase
     *
     * \param phaseIdx The index of the fluid phase to consider
     */
    static std::string phaseName(int phaseIdx)
    {
        assert(0 <= phaseIdx && phaseIdx < numPhases);

        return Fluid::name();
    }

    /*!
     * \brief Return whether a phase is liquid
     *
     * \param phaseIdx The index of the fluid phase to consider
     */
    static bool isLiquid(int phaseIdx)
    {
        assert(0 <= phaseIdx && phaseIdx < numPhases);

        return Fluid::isLiquid();
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
    static bool isCompressible(int phaseIdx)
    {
        assert(0 <= phaseIdx && phaseIdx < numPhases);

        // let the fluid decide
        return Fluid::isCompressible();
    }

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
    static bool isIdealGas(int phaseIdx)
    {
        assert(0 <= phaseIdx && phaseIdx < numPhases);

        // let the fluid decide
        return Fluid::isIdealGas();
    }

    /****************************************
     * Component related static parameters
     ****************************************/

    //! Number of components in the fluid system
    static constexpr int numComponents = 1;

    /*!
     * \brief Return the human readable name of a component
     *
     * \param compIdx The index of the component to consider
     */
    static std::string componentName(int compIdx)
    {
        assert(0 <= compIdx && compIdx < numComponents);

        return Fluid::name();
    }

    /*!
     * \brief Return the molar mass of a component in \f$\mathrm{[kg/mol]}\f$.
     *
     * \param compIdx index of the component
     */
    static Scalar molarMass(int compIdx)
    {
        assert(0 <= compIdx && compIdx < numComponents);

        return Fluid::molarMass();
    }

    /*!
     * \brief Critical temperature of a component \f$\mathrm{[K]}\f$.
     *
     * \param compIdx The index of the component to consider
     */
    static Scalar criticalTemperature(int compIdx)
    {
        assert(0 <= compIdx && compIdx < numComponents);

        return Fluid::criticalTemperature();
    }

    /*!
     * \brief Critical pressure of a component \f$\mathrm{[Pa]}\f$.
     *
     * \param compIdx The index of the component to consider
     */
    static Scalar criticalPressure(int compIdx)
    {
        assert(0 <= compIdx && compIdx < numComponents);

        return Fluid::criticalPressure();
    }

    /*!
     * \brief The acentric factor of a component \f$\mathrm{[-]}\f$.
     *
     * \param compIdx The index of the component to consider
     */
    static Scalar acentricFactor(int compIdx)
    {
        assert(0 <= compIdx && compIdx < numComponents);

        return Fluid::acentricFactor();
    }

    /****************************************
     * thermodynamic relations
     ****************************************/

    /*!
     * \brief Initialize the fluid system's static parameters
     */
    static void init()
    { }

    using Base::density;
    /*!
     * \brief Calculate the density \f$\mathrm{[kg/m^3]}\f$ of a fluid phase
     *
     * \param fluidState An arbitrary fluid state
     * \param phaseIdx The index of the fluid phase to consider
     */
    template <class FluidState>
    static Scalar density(const FluidState &fluidState,
                          int phaseIdx)
    {
        assert(0 <= phaseIdx && phaseIdx < numPhases);

        Scalar temperature = fluidState.temperature(phaseIdx);
        Scalar pressure = fluidState.pressure(phaseIdx);
        return Fluid::density(temperature, pressure);
    }

    using Base::viscosity;
    /*!
     * \brief Calculate the dynamic viscosity of a fluid phase \f$\mathrm{[Pa*s]}\f$
     *
     * \param fluidState An arbitrary fluid state
     * \param phaseIdx The index of the fluid phase to consider
     */
    template <class FluidState>
    static Scalar viscosity(const FluidState &fluidState,
                            int phaseIdx)
    {
        assert(0 <= phaseIdx && phaseIdx < numPhases);

        Scalar temperature = fluidState.temperature(phaseIdx);
        Scalar pressure = fluidState.pressure(phaseIdx);
        return Fluid::viscosity(temperature, pressure);
    }

    using Base::fugacityCoefficient;
    /*!
     * \brief Calculate the fugacity coefficient \f$\mathrm{[Pa]}\f$ of an individual
     *        component in a fluid phase
     *
     * The fugacity coefficient \f$\mathrm{\phi_\kappa_\alpha}\f$ is connected to the
     * fugacity \f$\mathrm{f^\kappa_\alpha}\f$ and the component's mole
     * fraction \f$\mathrm{x^\kappa_\alpha}\f$ by means of the relation
     *
     * \f[
     f^\kappa_\alpha = \phi^\kappa_\alpha\;x^\kappa_\alpha\;p_\alpha
     * \f]
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
        using std::numeric_limits;
        return numeric_limits<Scalar>::infinity();
    }

    using Base::diffusionCoefficient;
    /*!
     * \brief Calculate the molecular diffusion coefficient for a
     *        component in a fluid phase \f$\mathrm{[mol^2 * s / (kg*m^3)]}\f$
     *
     * Molecular diffusion of a compoent \f$\mathrm{\kappa}\f$ is caused by a
     * gradient of the chemical potential and follows the law
     *
     * \f[ J = - D \mathbf{grad} \mu_\kappa \f]
     *
     * where \f$\mathrm{\mu_\kappa}\f$ is the component's chemical potential,
     * \f$\mathrm{D}\f$ is the diffusion coefficient and \f$\mathrm{J}\f$ is the
     * diffusive flux. \f$\mathrm{\mu_\kappa}\f$ is connected to the component's
     * fugacity \f$\mathrm{f_\kappa}\f$ by the relation
     *
     * \f[ \mu_\kappa = R T_\alpha \mathrm{ln} \frac{f_\kappa}{p_\alpha} \f]
     *
     * where \f$\mathrm{p_\alpha}\f$ and \f$\mathrm{T_\alpha}\f$ are the fluid phase'
     * pressure and temperature.
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
     * \brief Given a phase's composition, temperature and pressure,
     *        return the binary diffusion coefficient \f$\mathrm{[m^2/s]}\f$ for components
     *        \f$\mathrm{i}\f$ and \f$\mathrm{j}\f$ in this phase.
     *
     * \param fluidState An arbitrary fluid state
     * \param phaseIdx The index of the fluid phase to consider
     * \param compIIdx The index of the first component to consider
     * \param compJIdx The index of the second component to consider
     */
    template <class FluidState>
    static Scalar binaryDiffusionCoefficient(const FluidState &fluidState,
                                             int phaseIdx,
                                             int compIIdx,
                                             int compJIdx)

    {
        DUNE_THROW(Dune::InvalidStateException, "Not applicable: Binary diffusion coefficients");
    }

    using Base::enthalpy;
    /*!
     * \brief Given a phase's composition, temperature, pressure and
     *        density, calculate its specific enthalpy \f$\mathrm{[J/kg]}\f$.
     *
     * \param fluidState An arbitrary fluid state
     * \param phaseIdx The index of the fluid phase to consider
     */
    template <class FluidState>
    static Scalar enthalpy(const FluidState &fluidState,
                           int phaseIdx)
    {
        assert(0 <= phaseIdx && phaseIdx < numPhases);

        Scalar temperature = fluidState.temperature(phaseIdx);
        Scalar pressure = fluidState.pressure(phaseIdx);
        return Fluid::enthalpy(temperature, pressure);
    }

    using Base::thermalConductivity;
    /*!
     * \brief Thermal conductivity of a fluid phase \f$\mathrm{[W/(m K)]}\f$.
     *
     * Use the conductivity of air and water as a first approximation.
     * Source:
     * http://en.wikipedia.org/wiki/List_of_thermal_conductivities
     *
     * \param fluidState An arbitrary fluid state
     * \param phaseIdx The index of the fluid phase to consider
     */
    template <class FluidState>
    static Scalar thermalConductivity(const FluidState &fluidState,
                                      int phaseIdx)
    {
        assert(0 <= phaseIdx && phaseIdx < numPhases);

        Scalar temperature = fluidState.temperature(phaseIdx);
        Scalar pressure = fluidState.pressure(phaseIdx);
        return Fluid::thermalConductivity(temperature, pressure);
    }

    using Base::heatCapacity;
    /*!
     * \brief Specific isobaric heat capacity of a fluid phase.
     *        \f$\mathrm{[J/(kg*K)]}\f$.
     *
     * \param fluidState An arbitrary fluid state
     * \param phaseIdx The index of the fluid phase to consider
     */
    template <class FluidState>
    static Scalar heatCapacity(const FluidState &fluidState,
                               int phaseIdx)
    {
        assert(0 <= phaseIdx && phaseIdx < numPhases);

        Scalar temperature = fluidState.temperature(phaseIdx);
        Scalar pressure = fluidState.pressure(phaseIdx);
        return Fluid::heatCapacity(temperature, pressure);
    }

};

} // end namespace

#ifdef DUMUX_PROPERTIES_HH
namespace Properties
{
NEW_PROP_TAG(Fluid);
}

/*!
 * \brief A pure single-phase fluid system.
 *
 * This is an adapter to use TwoPImmiscible<TypeTag>, as is
 * done with most other classes in Dumux and all template parameters
 * are usually defined in the property system anyhow.
 */
template<class TypeTag>
class OnePFluidSystem
: public FluidSystems::OneP<typename GET_PROP_TYPE(TypeTag, Scalar),
                            typename GET_PROP_TYPE(TypeTag, Fluid)>
{};
#endif
} // end namespace

#endif
